cd(fileparts(matlab.desktop.editor.getActiveFilename));
addpath(genpath('../src'))

if strcmp(getenv('USER'),'filippob')
    % -----------------------
    % NOTE: we require that those softwares are already installed on the system
    addpath(genpath('~/workSW/gurobi11.0.3_linux64/gurobi1103/linux64/matlab/'))
    addpath(genpath('~/workSW/casadi-3.7.0-linux64-matlab2018b'))
    run('~/workSW/tbxmanager/startup.m')

    % -----------------------
elseif strcmp(getenv("USERNAME"),'Filippo Badalamenti')
    addpath(genpath('C:\gurobi1103\win64\matlab'))
    run(getenv("HOMEPATH")+"\Documents\Research\workSW\tbxmanager\startup.m")
    addpath(genpath('C:\casadi-3.6.7-windows64-matlab2018b'))
else
    addpath(genpath('/Work_stuff/IMT_set_work/softwares/yalmip_latest'))
    addpath(genpath('/Work_stuff/IMT_set_work/Work_set_based/utilities'))
    addpath(genpath('/Work_stuff/IMT_set_work/softwares/mpt_files'))
    addpath(genpath('/Work_stuff/IMT_set_work/softwares/mosek'))
    rmpath(genpath('/Work_stuff/IMT_set_work/softwares/cplex'))
end
import casadi.*

% set variable environment to let CasaDi find Gurobi installation
setenv("GUROBI_VERSION","110");
setenv("GUROBI_HOME","/home/filippob/workSW/gurobi11.0.3_linux64/gurobi1103/linux64")
setenv("PATH",[getenv("GUROBI_HOME") '/bin:' getenv("PATH")])
setenv("LD_LIBRARY_PATH",[getenv("GUROBI_HOME") '/lib:' getenv("LD_LIBRARY_PATH")])

clc


%% System definition
% Dynamics
Ts = 0.1; g = 9.81; kT = 0.91;
[A,B,Bw] = defDiscreteQuadrotorDyn(Ts,g,kT);
nx = size(A,2); nu = size(B,2); nw = size(Bw,2);

% Constraints
% states: [x y z vx vy vz pitch roll wpitch wroll]
HX = [eye(nx); -eye(nx)]; hX = repmat([4,4,2,10,10,5,pi/3,pi/3,pi,pi]',[2,1]);
HU = [eye(nu); -eye(nu)]; hU = [pi/4; pi/4; 2*g-g/kT;  pi/4; pi/4; -(0-g/kT)];
HW = [eye(nw); -eye(nw)]; hW = repmat([0.05; 0.05; 0.1],[2,1]);

% Define LTI system
A_convh = {A}; B_convh = {B};

X = Polyhedron(HX, hX);
U = Polyhedron(HU, hU);
W_dist = Polyhedron(HW, hW);

sys = qLPV(A_convh, B_convh, Bw, X, U, W_dist);

%% CCpolytope
F_tilde = [-eye(nx); ones(1,nx)];
f_MPI = ones(size(F_tilde,1),1);
templateOpts = struct('computeRCI',true,'justFeasibleRCI',true, 'y_sigma', f_MPI);
ccPoly_base = CCPolytope(sys, F_tilde, templateOpts);

% compute cutting planes perpendicular to vertex normal
F_normal = zeros(ccPoly_base.v,nx);
for i=1:ccPoly_base.v
    F_normal(i,:) = mean(ccPoly_base.F(any(ccPoly_base.Vi_s{i},1),:),1);
end
rhs = Polyhedron(ccPoly_base.F, f_MPI).support(F_normal');
P_allcut = Polyhedron([ccPoly_base.F; F_normal],[f_MPI; 0.9*rhs]).minVRep;

% save the new simple CCPolytope
ccPoly = CCPolytope(sys, P_allcut.A, struct('computeRCI',false, 'y_sigma', P_allcut.b));

%% Define simple cost functions
% solve the mRCI set problem
mRCI_costMat = struct("Qy",1,"Qu",0.1);
mRCI_costFun = mRCICostFun(ccPoly, mRCI_costMat); % function signature: @(ys,us)
[y_m,u_m] = precompute_mRCI_set(ccPoly, mRCI_costFun);

% Tracking cost (define for CCTMPC and adapt to HTMPC using Equation(17))
Q_CC = blkdiag(1*eye(ccPoly.f),1*eye(ccPoly.v*sys.nu)); %(y,u(:))
gamma = 0.95;

% aggregate cost functions for CCTMPC and HTMPC
cfm_csd = aggregateCostFun(ccPoly, mRCI_costFun, Q_CC, gamma, y_m,u_m);

%% Build MPC scheme
csd_opts = {'gurobi',struct(),struct('outputflag',0)};%,'BarConvTol',1e-7,'FeasibilityTol',1e-9)};

N_ocp = 3; % Note: y_N = N_ocp+1
cctmpc = CCTMPC(sys, ccPoly, cfm_csd, N_ocp, gamma,csd_opts,false,y_m,u_m);

%% Simulation CCTMPC
N_mpc = 20*N_ocp+1; % simulation steps

% (state space) system dynamics
x_sys = zeros(nx, N_mpc+1);
x_sys(:,1) = (1-0e-6)*cctmpc.findFeasibleX0([4,-4,-2,zeros(1,7)]');%,diag([0.1,0.1,0.1,100*ones(1,7)]));
disp(x_sys(1:3,1)')

u_sys = zeros(nu, N_mpc);
% (parameter space) optimal control dynamics
y0_CL = cell(1,N_mpc);
Lyap_cost = zeros(1,N_mpc);

CC_times = nan(1,N_mpc);
% main loop
for t = 1:N_mpc
    disp("Iter: "+t+"/"+N_mpc)
    tic
    % compute input for the system
    u_sys(:,t) = cctmpc.solve(x_sys(:,t));
    CC_times(t) = toc;

    % propagate system dynamics
    x_sys(:,t+1) = sys.step(x_sys(:,t), u_sys(:,t));

    % save computed data
    y0_CL{t} = cctmpc.y0_sol;
    Lyap_cost(t) = cctmpc.Lyapunov_cost;
end
CC_times = CC_times(~any(isnan(CC_times), 2), :);
disp(min(CC_times,[],"all")*1000)
disp(mean(CC_times,"all")*1000)
disp(max(CC_times,[],"all")*1000)

timeCC = mean(CC_times,"all")*1000;

%% compute Baseline Hausdorff distance
baseHausCC = cctmpc.hausDist;

%% Simulation HTMPC (+ Hausdorff distance)
horIncr = 0:7;
hausH = zeros(1,length(horIncr));
H_s_times = nan(length(horIncr),N_mpc);

% use a simple tracking cost
Q_H = blkdiag(eye(nx),eye(nu),1); %(z,v,a)
cfm_csd = aggregateCostFun_fromH(ccPoly, mRCI_costFun, Q_H, gamma, y_m,u_m);

for n=1:length(horIncr)
    htmpc = HTMPC(sys, ccPoly, cfm_csd, N_ocp+horIncr(n), gamma,csd_opts,false,y_m,u_m);
    hausH(n) = htmpc.hausDist;
    % (state space) system dynamics
    x_sys = zeros(nx, N_mpc+1);
    x_sys(:,1) = (1-0e-6)*htmpc.findFeasibleX0([4,-4,-2,zeros(1,7)]');%,diag([0.1,0.1,0.1,1000*ones(1,7)]));
    disp(x_sys(1:3,1)')

    u_sys = zeros(nu, N_mpc);

    % (parameter space) optimal control dynamics
    y0_CL = cell(1,N_mpc);
    Lyap_cost = zeros(1,N_mpc);

    % main loop
    for t = 1:N_mpc
        disp("Iter: "+t+"/"+N_mpc)
        tic
        % compute input for the system
        u_sys(:,t) = htmpc.solve(x_sys(:,t));
        H_s_times(n,t) = toc;

        % propagate system dynamics
        x_sys(:,t+1) = sys.step(x_sys(:,t), u_sys(:,t));

        % save computed data
        y0_CL{t} = htmpc.y0_sol;
        Lyap_cost(t) = htmpc.Lyapunov_cost;
    end
    % % stop at a certain point to collect data
    % if (hausH(max(n-2,1)) < baseHausCC); break; end
end
H_s_times = H_s_times(~any(isnan(H_s_times), 2), :);
timeH_s = mean(H_s_times,2)*1000;

%% Plot section Case study 1
plot_CC_H_comparison;

%% Define the new, more complex template for HTMPC

F_box = [eye(nx); -eye(nx)];
try % attempt to restore previously computed template
    load("Quadrotor_box.mat")
catch % warning, long computational time
    box_templOpts = struct('computeRCI',true, 'justFeasibleRCI',true, 'y_sigma', ones(2*nx,1));
    box_ccPoly = CCPolytope(sys, F_box, box_templOpts);
    save("Quadrotor_box.mat","box_ccPoly")
end

% % adapt Cost function and mRCI set
mRCI_costMat = struct("Qy",1,"Qu",0.1);
mRCI_costFun = mRCICostFun(box_ccPoly, mRCI_costMat); % function signature: @(ys,us)
[y_box,u_box] = precompute_mRCI_set(box_ccPoly, mRCI_costFun);

% Simple Tracking cost for HTMPC
Q_H = blkdiag(eye(nx),eye(nu),1); %(z,v,a)

gamma = 0.95;

% aggregate cost functions for CCTMPC and HTMPC
cfm_csd = aggregateCostFun_fromH(box_ccPoly, mRCI_costFun, Q_H, gamma, y_box,u_box);

%%
box_N_ocp = 20; % Note: y_N = N_ocp+1
new_htmpc = HTMPC(sys, box_ccPoly, cfm_csd, box_N_ocp, gamma, csd_opts,false, y_box,u_box);

%% Simulation HTMPC
% find a good initial condition
disp(new_htmpc.hausDist);
x0_space = new_htmpc.convh_x0s(1:3,:);

[~,ind] = maxk(vecnorm(x0_space),5);
% or look at the plot
%%
figure; plot(1:box_ccPoly.v+box_ccPoly.f,x0_space); xlim([1,box_ccPoly.v+box_ccPoly.f])
%%
N_mpc = 5*box_N_ocp+1; % simulation steps

% (state space) system dynamics
x_sys = zeros(nx, N_mpc+1);
x_sys(:,1) = (1-1e-5)*new_htmpc.convh_x0s(:,572);
disp(x_sys(1:3,1)')

u_sys = zeros(nu, N_mpc);
% (parameter space) optimal control dynamics
y0_CL = cell(1,N_mpc); z0_CL = cell(1,N_mpc); alpha0_CL = cell(1,N_mpc);
Lyap_cost = zeros(1,N_mpc);

Hbox_times = nan(1,N_mpc);
% main loop
for t = 1:N_mpc
    disp("Iter: "+t+"/"+N_mpc)
    tic
    % compute input for the system
    u_sys(:,t) = new_htmpc.solve(x_sys(:,t));
    Hbox_times(t) = toc;

    % propagate system dynamics
    x_sys(:,t+1) = sys.step(x_sys(:,t), u_sys(:,t));

    % save computed data
    y0_CL{t} = new_htmpc.y0_sol;
    [z0_CL{t},alpha0_CL{t}] = new_htmpc.ocpSol{1:2};
    Lyap_cost(t) = new_htmpc.Lyapunov_cost;
end
Hbox_times = Hbox_times(~any(isnan(Hbox_times), 2), :);
timeHbox = mean(Hbox_times,"all")*1000;

%% Plot section Case study 2

plot_QuadrotorBox;























%% ------------------------- Support functions ----------------------------

function [A_d, B_d, Bw_d] = defDiscreteQuadrotorDyn(Ts,g,kT)

% get Continuous time Quadrotor model
[A_num, B_num, Bw] = getQuadrotorModel(g,kT);

% discretize it using Euler Approximation
% (we don't need a very accurate model, only a test one in discrete time)
A_d = eye(size(A_num,1)) + Ts*A_num;
B_d = Ts*B_num;
Bw_d = Ts*Bw;

end

function [A_num, B_num, Bw] = getQuadrotorModel(g,kT)
d0 = 10; d1 = 8; % coefficients
n0 = 10;

syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 u1 u2 u3 real;

% Define symbolic state and input vectors.
x = [x1; x2; x3; x4; x5; x6; x7; x8; x9; x10];
u = [u1; u2; u3];

% Define the dynamics f(x,u) (set disturbances w=0).
f1 = x4; % + w(1);
f2 = x5; % + w(2);
f3 = x6; % + w(3);
f4 = g * tan(x7);
f5 = g * tan(x8);
f6 = -g;  % note: will be canceled by kT*u3 at equilibrium.
f7 = -d1 * x7 + x9;
f8 = -d1 * x8 + x10;
f9 = -d0 * x3;
f10 = -d0 * x7;

% Combine into f.
f = [f1; f2; f3; f4; f5; f6; f7; f8; f9; f10];

% Define the control influence matrix G symbolically.
G = zeros(10,3);
G(6,3) = kT;
G(9,1) = n0;
G(10,2) = n0;

f = f + G*u;  % Full dynamics

% Compute the Jacobians.
A_sym = jacobian(f, x);
B_sym = jacobian(f, u);

% Define the equilibrium (for instance: all states zero, but u3 = g/kT).
x_eq = zeros(10,1);
u_eq = [0; 0; g/kT];

% Evaluate A_sym:
A_num = double(subs(A_sym, [x; u], [x_eq; u_eq]));
% Similarly for B_sym:
B_num = double(subs(B_sym, [x; u], [x_eq; u_eq]));

Bw = [eye(3); zeros(7,3)];
end

