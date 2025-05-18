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
h = 1/4;

A = [1 h (h^2)/2; 0 1 h; 0 0 1];                nx = size(A,2);
B = [(h^3)/6; (h^2)/2; h];                      nu = size(B,2);
Bw = [h (h^2)/2 (h^3)/6; 0 h (h^2)/2; 0 0 h];   nw = size(Bw,2);

% Constraints
HX = [eye(nx); -eye(nx)]; hX = 5*ones(2*nx,1);
HU = [eye(nu); -eye(nu)]; hU = [3; 3];
HW = [eye(nw); -eye(nw)]; hW = (1/20)*ones(2*nw,1);


% Define LTI system
A_convh = {0.9*A, 1.1*A}; B_convh = {0.9*B, 1.1*B};

X = Polyhedron(HX, hX);
U = Polyhedron(HU, hU);

W_dist = Polyhedron(HW, hW);

sys = qLPV(A_convh, B_convh, Bw, X, U, W_dist);

%% CCpolytope

% % define a template from unit (nx-1)-sphere sampling
% f = 6; %nx+1
% F = Polyhedron(numericalSpherePoints(nx, f)).dual.minHRep().A;
% ccPoly = CCPolytope(sys, F);% struct('computeRCI',false,'y_sigma', ones(f,1)));

% or use the template from CTR_example.m
chosenIter = 10+1;
ccPoly = CCPolytope(sys, save_Triple{chosenIter});

% for graphical purposes (works in 2D and 3D)
% plotVerticesSimplePolytope(ccPoly, save_yM{chosenIter}); 

%%
% Cost Function definition: we can either pass structs of matrices, or
% function handles defined externally. We only require that the cost
% functions are -strictly- convex, and that the function signatures are:
% RCI : @(ys,us)
% Stage/Term_cost : @(y_k/N,u_k/N,ys,us)
Qv = blkdiag(0.1*eye(nx), 0.1*eye(nu));
Qc = blkdiag(1*eye(nx), 1*eye(nu));
mRCI_costMat = struct("Qv",Qv,"Qc",Qc);

% mRCI_costMat = struct("Qy",10,"Qu",0.1);
mRCI_costFun = mRCICostFun(ccPoly, mRCI_costMat); % function signature: @(ys,us)
[y_m,u_m] = precompute_mRCI_set(ccPoly, mRCI_costFun);

% Tracking cost (define for CCTMPC and adapt to HTMPC using Equation(17))
% Q_H = blkdiag(eye(nx),eye(nu),1); %(z,v,a)
Q_CC = blkdiag(1*eye(ccPoly.f), 1*eye(ccPoly.v*sys.nu)); %(y,u(:))

gamma = 0.95;

% aggregate cost functions for CCTMPC and HTMPC
cfm_csd = aggregateCostFun(ccPoly, mRCI_costFun, Q_CC, gamma, y_m,u_m);
% cfm_csd = aggregateCostFun_fromH(ccPoly, mRCI_costFun, Q_H, gamma, y_m,u_m);

%% Build MPC scheme
csd_opts = {'daqp'}; %{name_solver,CasadiOpti_opts,internal_solver_opts}
% csd_opts = {'gurobi',struct(),struct('outputflag',1)};

N_ocp = 3; % Note: y_N = N_ocp+1
cctmpc = CCTMPC(sys, ccPoly, cfm_csd, N_ocp, gamma,csd_opts,false,y_m,u_m);
htmpc = HTMPC(sys, ccPoly, cfm_csd, N_ocp, gamma, csd_opts,false,y_m,u_m);

%% Simulation CCTMPC
N_mpc = 5*N_ocp+1; % simulation steps
x_inits = cartesianProduct(linspace(-5,5,3),3);
x_inits(ismembertol(x_inits, zeros(1,3), 1e-5,'ByRows',true),:) = [];

x_sys_all = cell(1,size(x_inits,1));

CC_times = nan(size(x_inits,1),N_mpc);
for i=1:size(x_inits,1)
    disp("MPC simulation: " + i +"/"+ size(x_inits,1))

    % (state space) system dynamics
    x_sys = zeros(nx, N_mpc+1);
    x_sys(:,1) = (1-1e-5)*cctmpc.findFeasibleX0(x_inits(i,:)');
    u_sys = zeros(nu, N_mpc);

    % (parameter space) optimal control dynamics
    y0_CL = cell(1,N_mpc);
    Lyap_cost = zeros(1,N_mpc);

    % main loop
    for t = 1:N_mpc
        disp("Iter: "+t+"/"+N_mpc)
        tic
        % compute input for the system
        u_sys(:,t) = cctmpc.solve(x_sys(:,t));
        CC_times(i,t) = toc;

        % propagate system dynamics
        x_sys(:,t+1) = sys.step(x_sys(:,t), u_sys(:,t));

        % save computed data
        y0_CL{t} = cctmpc.y0_sol;
        Lyap_cost(t) = cctmpc.Lyapunov_cost;
    end

    x_sys_all{i} = x_sys;
end
%%
CC_times = CC_times(~any(isnan(CC_times), 2), :);
disp(min(CC_times,[],"all")*1000)
disp(mean(CC_times,"all")*1000)
disp(max(CC_times,[],"all")*1000)

%% Simulation HTMPC

H_times = nan(size(x_inits,1),N_mpc);
for i=1:size(x_inits,1)
    disp("MPC simulation: " + i +"/"+ size(x_inits,1))

    % (state space) system dynamics
    x_sys = zeros(nx, N_mpc+1);
    x_sys(:,1) = (1-1e-5)*htmpc.findFeasibleX0(x_inits(i,:)');
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
        H_times(i,t) = toc;

        % propagate system dynamics
        x_sys(:,t+1) = sys.step(x_sys(:,t), u_sys(:,t));

        % save computed data
        y0_CL{t} = htmpc.y0_sol;
        Lyap_cost(t) = htmpc.Lyapunov_cost;
    end
end
%%
H_times = H_times(~any(isnan(H_times), 2), :);
disp(min(H_times,[],"all")*1000)
disp(mean(H_times,"all")*1000)
disp(max(H_times,[],"all")*1000)


%% Compute Regions of Attraction, and compare them with approx. MRCI set
RoA_CCTMPC = cctmpc.computeFeasRegion(200); disp("CC_done.")
RoA_HTMPC = htmpc.computeFeasRegion(200); disp("H_SP_done.")

approxMRCI = ULTISystem('A', A_convh, 'B', B_convh, 'E',Bw...
    ).invariantSet('X',X,'U',U,'D',Polyhedron(HW,hW),'maxIterations',6);
%%
dist_CCTMPC_MRCI = cctmpc.getHausdorffDist(approxMRCI);
dist_HTMPC_SP_MRCI = htmpc.getHausdorffDist(approxMRCI);

disp("Hausdorff distance (X,CCTMPC): " + round(dist_CCTMPC_MRCI,4))
disp("Hausdorff distance (X,HTMPC_SP): " + round(dist_HTMPC_SP_MRCI,4))


%% Plot section
plot_Integrator3D_ss;

%% ------------------------- Support functions ----------------------------

function plotVerticesSimplePolytope(ccPoly,y_sigma) %#ok<DEFNU>
figure; hold on
Polyhedron(ccPoly.F,y_sigma).plot("alpha",0.1,"EdgeColor",[0.9290 0.6940 0.1250],...
    "Color",[0.9290 0.6940 0.1250],"Linewidth",0.5,"LineStyle",':');
ccPoly.sys.X.plot("alpha",0,"EdgeColor",0.3*ones(1,3));

for v=1:ccPoly.v
    x = ccPoly.Vi_s{v}*y_sigma;
    x = num2cell(x);
    text(x{:}, sprintf('$V_{%d} y$', v), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'FontSize', 12,"Interpreter","latex");
end
xlabel("x")
ylabel("y")
end


function combos = cartesianProduct(v, k)
% generate the Cartesian product of the vector v taken k times.

% Create a cell array where each cell contains the vector v.
args = repmat({v}, 1, k);

% Use ndgrid to generate k-dimensional grids.
[G{1:k}] = ndgrid(args{:});

% Determine the total number of combinations.
numPoints = numel(G{1});

combos = zeros(numPoints, k);
% Flatten each grid and put it as a column.
for i = 1:k
    combos(:, i) = G{i}(:);
end
end
