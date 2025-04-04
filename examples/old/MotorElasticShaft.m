addpath(genpath('../src'))

if strcmp(getenv('USER'),'filippob')
    % -----------------------
    % NOTE: we require that those softwares are already installed on the system
    addpath(genpath('~/workSW/gurobi1100/linux64/matlab/'))
    run('~/Documenti/CCTMPC/tbxmanager/startup.m')

    % -----------------------
elseif strcmp(getenv("USERNAME"),'Filippo Badalamenti')
    addpath(genpath('C:\gurobi1201\win64\matlab'))
    run(getenv("HOMEPATH")+"\Documents\Research\workSW\tbxmanager\startup.m")
    addpath('C:\casadi-3.6.7-windows64-matlab2018b')
else
    addpath(genpath('/Work_stuff/IMT_set_work/softwares/yalmip_latest'))
    addpath(genpath('/Work_stuff/IMT_set_work/Work_set_based/utilities'))
    addpath(genpath('/Work_stuff/IMT_set_work/softwares/mpt_files'))
    addpath(genpath('/Work_stuff/IMT_set_work/softwares/mosek'))
    rmpath(genpath('/Work_stuff/IMT_set_work/softwares/cplex'))
end

%% System definition

Ts = 0.1;  % seconds
[A_convh, B_convh, E_convh] = defDiscreteMotorModel(Ts);
Bw = E_convh{1};

nx = size(A_convh{1},2); nu = size(B_convh{1},2); nw = size(Bw,2);

% Output
C = blkdiag(eye(nx));
ny = size(C,1);
D = zeros(ny, nu);

% State and Input Bounds
% State vector: x = [i_a; i_f; omega_m; theta_m; omega_l; theta_l; delta; dot_delta] (8 states)
x_lb = [-15; -10; -200; -0.5; -200; -0.5; -0.2; -50];   % lower bounds [A, A, rad/s, rad, rad/s, rad, rad, rad/s]
x_ub = [15; 10; 200; 0.5; 200; 0.5; 0.2; 50];           % upper bounds

% Input vector: u = [v_a; v_f] (2 inputs)
u_lb = [-24; -24];   % lower bounds [V]
u_ub = [24; 24];     % upper bounds [V]

HX = [eye(nx); -eye(nx)];   hX = [x_ub; -x_lb];              
HU = [eye(nu); -eye(nu)];   hU = [u_ub; -u_lb];  

% Disturbance vector: d = [d_motor; d_load] (2-dimensional)
d_lb = [-0.5; -0.5];   % lower bounds [N*m]
d_ub = [ 0.5;  0.5];   % upper bounds [N*m]

HW = [eye(nw); -eye(nw)]; hW = [d_ub; -d_lb];        

% Constraints
X = Polyhedron(HX, hX);
U = Polyhedron(HU, hU);

W_dist = Bw*Polyhedron(HW, hW);

sys = qLPV(A_convh, B_convh, C, D, X, U, W_dist);
%%
% search for different dimensions a good template
for f=nx+1:nx+10
F = Polyhedron(numericalSpherePoints(nx, f, 1),ones(f,1)).minHRep();
disp("f: "+ size(F.A,1) + ", v: " + size(uniquetol(F.V,1e-10,'ByRows',true),1))
end
%%
f = 14;
% take the polar set, having probability~=1 to have simple polytope (as we compute a simplicial one)
F = Polyhedron(numericalSpherePoints(nx, f,1),ones(f,1)).minHRep();
ccPoly = CCPolytope_test(sys, F.A, true);

%%
% Cost Function definition: we can either pass structs of matrices, or
% function handles defined externally. We only require that the cost
% functions are -strictly- convex, and that the function signatures are:
% RCI : @(ys,us)
% Stage/Term_cost : @(y_k/N,u_k/N,ys,us)
paramTrack = true;

% Qv = blkdiag(0.01*eye(nx), 10*eye(nu));
% Qc = blkdiag(1*eye(nx), 1*eye(nu)); 
% RCI_cost = struct("Qv",Qv,"Qc",Qc);

RCI_cost = struct("Qy",0.1,"Qu",10);
% input cost greatly influences ROA of HTMPC_* methods

gamma = 0.95;
if paramTrack
    Qmat = blkdiag(10*eye(ccPoly.f), 0.1*eye(sys.nu*ccPoly.v));
    Pmat = (1/(1-gamma^2))*Qmat;
else
    [Qcost_han, Pcost_han] = defineCustomTrackCost(sys, ccPoly, Qv, gamma); %#ok<UNRCH>
end

% collect all cost function definition in a single class
costFunMan = CostFunctionManager(ccPoly, Qmat, Pmat, RCI_cost);


% Build MPC scheme
sdp_opts = {'solver','gurobi','verbose',5};
ocp_opts = struct("useBasicOCP", false, "var_convh", false);

N_ocp = 3; % Note: y_N = N_ocp+1
cctmpc = CCTMPC(sys, ccPoly, costFunMan, N_ocp, gamma, sdp_opts, ocp_opts);
htmpc = HTMPC(sys, ccPoly, costFunMan, N_ocp, gamma, sdp_opts, ocp_opts);
htmpc_sp = HTMPC_SinglePolicy(sys, ccPoly, costFunMan, N_ocp, gamma, sdp_opts, ocp_opts);

%%
% first show the vertex V_s{i}*ones(f,1) in the plots:
plotVerticesSimplePolytope(ccPoly); % works in 2D and 3D

%%

% then, specify how many (or which) vertices do you want to fix. 
% H_vertices = 1:ccPoly.v; % maybe also comparing selecting vertices vs random params?
H_vertices = [];
% H_vertices = [1,6];
part_cctmpc = PartialCCTMPC(sys, ccPoly, H_vertices,costFunMan, N_ocp, gamma, sdp_opts, ocp_opts);
%% Simulation

% % for debug reasons, N_mpc = k*T_orbit+1, k \in N
% N_mpc = 2*N_ocp+1; % simulation steps
% 
% % (state space) system dynamics
% x_sys = zeros(nx, N_mpc+1); x_sys(:,1) = [-5; 4]; % initial condition
% u_sys = zeros(nu, N_mpc);
% w_sys = zeros(nx, N_mpc);
% 
% 
% % (parameter space) optimal control dynamics
% OCP_y = cell(1,N_mpc); OCP_u = cell(1,N_mpc);
% OCP_ys = cell(1,N_mpc); OCP_us = cell(1,N_mpc);
% 
% Lyap_cost = zeros(1,N_mpc);
% 
% % main loop
% for t = 1:N_mpc
%     disp("MPC iteration: " + t +"/"+ N_mpc)
% 
%     % compute input for the system
%     u_sys(:,t) = cctmpc.solve(x_sys(:,t));
% 
%     % compute worst case disturbance (with respect to the origin)
%     % w_sys(:,t) = getWorstCaseDist(sys, x_sys(:,t), u_sys(:,t));
% 
%     % propagate system dynamics
%     % x_sys(:,t+1) = sys.step_nominal(x_sys(:,t), u_sys(:,t)) + w_sys(:,t);
%     x_sys(:,t+1) = sys.step(x_sys(:,t), u_sys(:,t));
% 
%     % save computed data
%     [OCP_y{t}, OCP_u{t}, OCP_ys{t}, OCP_us{t}] = cctmpc.ocpSol{1:4};
%     Lyap_cost(t) = cctmpc.Lyapunov_cost;
% end
% 
% Integrator2D_Lyapunov;

%% Compute the feasibility region for all TMPC schemes, and compare them with
% NOTE: what to show? inner (convex hull of x_i) or outer (hyperplane) approximation? 
feasRegionCCTMPC = intersect(cctmpc.computeFeasRegion(10),X); disp("CC_done.")
feasRegionHTMPC = intersect(htmpc.computeFeasRegion(10),X); disp("H_done.")
feasRegionHTMPC_SP = intersect(htmpc_sp.computeFeasRegion(10),X); disp("H_SP_done.")
feasRegionPartCCTMPC = intersect(part_cctmpc.computeFeasRegion(10),X); disp("PartCC_done.")

% NOTE: to get the same result, you can alternatively use:
% MRCI = ULTISystem('A', A_convh, 'B', B_convh).invariantSet('X',X,'U',U,'D',Polyhedron(HW,hW));
MRCI = feasRegionCCTMPC;
% compare set distance of feasible regions from the MRCI
dist_CCTMPC_MRCI = setDistance(feasRegionCCTMPC, MRCI);
dist_HTMPC_MRCI = setDistance(feasRegionHTMPC, MRCI);
dist_HTMPC_SP_MRCI = setDistance(feasRegionHTMPC_SP, MRCI);
dist_PartCCTMPC_MRCI = setDistance(feasRegionPartCCTMPC, MRCI);


disp("Hausdorff distance (MRCI,CCTMPC): " + round(dist_CCTMPC_MRCI,4))
disp("Hausdorff distance (MRCI,HTMPC): " + round(dist_HTMPC_MRCI,4))
disp("Hausdorff distance (MRCI,HTMPC_SP): " + round(dist_HTMPC_SP_MRCI,4))
disp("Hausdorff distance (MRCI,PartCCTMPC): " + round(dist_PartCCTMPC_MRCI,4))
disp(" ")
disp("Is PartCCTMPC equal to CCTMPC: " + mat2str(feasRegionCCTMPC== feasRegionPartCCTMPC))
disp("Is PartCCTMPC equal to HTMPC: " + mat2str(feasRegionHTMPC== feasRegionPartCCTMPC))

disp("Volume coverage: " + round((feasRegionPartCCTMPC.volume/feasRegionCCTMPC.volume)*100,4) + "%")

%%
fig = figure("Renderer","painters","Units","centimeters");
fig.Position(3:4) = 4*[3,2.5]; % [width,height]
hold on

X.plot("Alpha",0,"EdgeColor",0.4*ones(1,3),"Linewidth",1.5)
han_CC = feasRegionCCTMPC.plot("alpha",0.1,"EdgeColor",[0.4660 0.6740 0.1880],"Color",[0.4660 0.6740 0.1880],"Linewidth",1.25);

han_PartCC = feasRegionPartCCTMPC.plot("alpha",0.1,"EdgeColor",[0 0.4470 0.7410],"Color",[0 0.4470 0.7410],"Linewidth",1.25);

han_H = feasRegionHTMPC.plot("Alpha",0.1,"EdgeColor",[0.8500 0.3250 0.0980],"Color",[0.8500 0.3250 0.0980],"Linewidth",1.25);
% 
han_SP = feasRegionHTMPC_SP.plot("Alpha",0.1,"EdgeColor",[0.9290 0.6940 0.1250],"Color",[0.9290 0.6940 0.1250],"Linewidth",1.25);
% 
han_mRPI = Polyhedron(ccPoly.F,htmpc_sp.rpiSol{1}).plot("Alpha",0.8,"Color",0.9*[0.9290 0.6940 0.1250]);
han_mRCI = Polyhedron(ccPoly.F,cctmpc.rciSol{1}).plot("Alpha",0.8,"Color",0.9*[0.4660 0.6740 0.1880]);

% since we have a 3D integrator, we can get the slice at x3=0 and show in
% perspective!

% % bring grid in front of everything
Ax1 = gca;
Ax1.YGrid = 'on'; Ax1.Layer = 'top'; Ax1.GridAlpha = 0.05;
set(Ax1, 'TickLabelInterpreter', 'latex');

han_leg1 = legend(Ax1,[han_CC, han_PartCC,han_H, han_SP,han_mRPI,han_mRCI],... han_ref],...
    {'CCTMPC','PartCCTMPC','HTMPC','HTMPC\_SP','mRPI','mRCI'}, ...
    'Interpreter','latex','Location','northeast');
title("Feasible regions and minimal Robust sets comparison")


% % minimize white borders around plot
set(Ax1,'LooseInset', max(get(Ax1,'TightInset'), 0.01)) % remove border from axis
set(fig,'PaperPositionMode','Auto','PaperUnits',...
    'centimeters','PaperSize',fig.Position(3:4)) % resize pdf page

fig.Renderer = 'painters'; % a way to force saving in vector graphics
% % save the plot as PDF file
saveas(fig, '../figures/Integrator3D_Regions.pdf','pdf')

%% Plot section
% generate plots in separate scripts to improve code readibility
Integrator2D_PhasePlane;
%%
step2plot = 13;
Integrator2D_StepOpenLoop;
%%

%% ------------------------- Support functions ----------------------------

function [Qtrack_han, Ptrack_han] = defineCustomTrackCost(sys,ccPoly,Qv,gamma)
% We provide an instance of a user-specified cost function, passed as an
% handle. In this case, we seek a tradeoff between tracking performance and
% CCPolytope volume. See "Illustrative Example" from Section IV of the paper.

nx = sys.nx; nu = sys.nu;
f = ccPoly.f; v = ccPoly.v;

% RCI cost matrices (in state space)
Qv_x = Qv(1:nx,1:nx);
Qv_u = Qv(nx+1:end,nx+1:end);

V_bar = mean(cat(3, ccPoly.Vi_s{:}), 3); % each Vi_s{i} is a 2D matrix, averaging along 3rd dim.
U_bar = (1/v)*repmat(eye(nu),1,v); % nu x (nu*v)

% cumulating matrix cost (in parameter space)
Qv_yp = zeros(f,f);
for i = 1:v
    Qv_yp = Qv_yp+(ccPoly.Vi_s{i}-V_bar)'*Qv_x*(ccPoly.Vi_s{i}-V_bar);
end

Qv_up = zeros(nu*v, nu*v);
Iu_mats = cell(1,v);
for i = 1:v
    Iu_mats{i} = sparse(nu,nu*v);
    Iu_mats{i}(:,(i-1)*nu+1:i*nu) = eye(nu);
    Qv_up = Qv_up+(Iu_mats{i}-U_bar)'*Qv_u*(Iu_mats{i}-U_bar);
end

Q = blkdiag(Qv_yp,Qv_up)+1e-5*eye(f+v*nu,f+v*nu);
P = (1/(1-gamma^2))*Q;

% define function handles (here quadratic costs)
Qtrack_han = @(y_k,u_k,ys,us) weighted2NormSquared([y_k-ys; u_k(:)-us(:)],Q);
Ptrack_han = @(y_N,u_N,ys,us) weighted2NormSquared([y_N-ys; u_N(:)-us(:)],P);

end



function w_dist = getWorstCaseDist(sys,x_curr,uMPC) %#ok<DEFNU>
% compute worst case disturbance (with respect to the origin) by
% enumerating vertices of disturbance set
W_dist_vert = sys.W_dist.V'; n_d = size(W_dist_vert,2);

x_next = zeros(sys.nx,n_d);
dist_to_orig = zeros(1,n_d);
for j=1:n_d % enumerate disturbance realizations
    x_next(:,j) = sys.A_curr*x_curr + sys.B_curr*uMPC + W_dist_vert(:,j);
    dist_to_orig(j) = norm(x_next(:,j),2)^2;
end
[~,id_worst] = max(dist_to_orig);
w_dist = W_dist_vert(:,id_worst);

end

function dist = setDistance(polytope1, polytope2)
% compute the l-norm Hausdorff distance between polytopes
vert_poly2 = polytope2.V'; % (nx,v)
[nx, v] = size(vert_poly2);

x_pts = sdpvar(nx, v,'full');
unit_ball_pts = sdpvar(nx, v,'full');
eps = sdpvar(1,1,'full');

constr = [vert_poly2 == x_pts + unit_ball_pts;
    polytope1.A*x_pts <= repmat(polytope1.b, 1, v);
    [eye(nx); -eye(nx)]*unit_ball_pts <= repmat(eps*ones(2*nx,1),1,v)];

optimize(constr, eps, sdpsettings('verbose',0));
dist = value(eps);

end

function plotVerticesSimplePolytope(ccPoly)
figure; hold on
ccPoly.ccpolytope.projection(1:2).plot("alpha",0.1,"EdgeColor",[0.9290 0.6940 0.1250],...
    "Color",[0.9290 0.6940 0.1250],"Linewidth",0.5,"LineStyle",':');

for v=1:ccPoly.v
    x = ccPoly.Vi_s{v}*ones(ccPoly.f,1);
    x = num2cell(x);
    text(x{:}, sprintf('$V_{%d} y$', v), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'FontSize', 12,"Interpreter","latex");
end
end

function combinations = allUnitCombs(n)
% Total number of combinations
numComb = 2^n;
% Initialize the matrix to store combinations
combinations = zeros(numComb, n);

% Generate all combinations using binary representation
for i = 0:numComb-1
    % Convert the integer to a binary vector
    % The 'left-msb' option ensures the most significant bit is on the left
    binVec = double(bitget(i, n:-1:1));
    % Map binary digits to -1 and +1 (0 -> -1, 1 -> +1)
    combRow = binVec * 2 - 1;
    % Store the combination in the matrix
    combinations(i+1, :) = combRow;
end

% Display the result
% disp(combinations);
end


function [Ad, Bd, Ed] = eulerDiscretize(A, B, E, Ts)
    % Euler discretization:
    % Ad = I + Ts*A
    % Bd = Ts*B
    % Ed = Ts*E
    Ad = eye(size(A)) + Ts * A;
    Bd = Ts * B;
    Ed = Ts * E;
end

function [A_convexHull, B_convexHull, E_convexHull] = discretizeConvexHull(A0, A1, B, E, Ts)
    % Compute discrete matrices for vertex 1 (lambda = 0): A0
    [Ad1, Bd1, Ed1] = eulerDiscretize(A0, B, E, Ts);
    
    % Compute discrete matrices for vertex 2 (lambda = 1): A0 + A1
    [Ad2, Bd2, Ed2] = eulerDiscretize(A0 + A1, B, E, Ts);
    
    % Pack the results into cell arrays
    A_convexHull = {Ad1, Ad2};
    B_convexHull = {Bd1, Bd2};
    E_convexHull = {Ed1, Ed2};
end

function [A_convh, B_convh, E_convh] = defDiscreteMotorModel(Ts)
%% Parameters (nominal and uncertain vertices)
% Define physical parameters
Jm = 0.01;      % Motor inertia
Jl = 0.02;      % Load inertia
% Other parameters (R_a, L_a etc.) are defined according to your model,
% but for discretization we focus on the state-space matrices.

% Uncertainty in shaft stiffness K_s is modeled via the matrices A0 and A1.
% (Here we assume that the uncertain term enters linearly.)
% For illustration, assume that A0 includes the effect of K_s = 90 N*m/rad,
% while A1 is the incremental effect altering K_s to 110.
% (In practice, you'd derive A0 and A1 from your physical equations.)
%
% An example A0 (8x8) might be:
A0 = [ ...
    0    0    0    0    0     0    0    0;
    0    0    0    0    0     0    0    0;
    0    0   -50  100    0     0    0    0;  % motor dynamics row (example)
    0    0    1    0     0     0    0    0;  % motor position
    0    0    0    0  -30    60    0    0;  % load dynamics row (example)
    0    0    0    0    1     0    0    0;  % load position
    0    0    0    0    0     0    0    1;  % flexible mode: twist derivative
    0    0    0    0    0     0    0    0]; % flexible mode dynamics

% A1 captures the effect of the uncertainty in K_s
A1 = [ ...
    zeros(4,8);  % first four rows not affected
    zeros(2,8);  % next two rows not affected
    zeros(1,8);  % row 7
    zeros(1,8)]; % row 8
% For example, let's assume only the motor and load dynamic rows are affected:
A1(3,3) = -5;      % slight additional damping in motor dynamics
A1(5,3) =  3;      % variation affecting load dynamics

% Define the input matrix B (8x2), assuming the inputs affect the mechanical side:
B = [...
    zeros(2,2);       % electrical dynamics (if any)
    zeros(1,2);
    zeros(1,2);
    zeros(1,2);
    zeros(1,2);
    zeros(1,2);
    zeros(1,2)];  
% Here we provide a simple example: let the voltages control motor & load acceleration
% Replacing some zeros with actual values.
B(3,1) = 1;      % motor acceleration affected by first input
B(5,2) = 1;      % load acceleration affected by second input

% Define the disturbance matrix E (8x2)
E = zeros(8,2);
E(3,1) = 1/Jm;   % disturbance on motor acceleration
E(5,2) = 1/Jl;   % disturbance on load acceleration

%% Discretization of the Convex Hull
[A_convh, B_convh, E_convh] = discretizeConvexHull(A0, A1, B, E, Ts);

% % Display the discrete matrices for each vertex (for verification)
% disp('Discrete A matrices for the convex hull:');
% disp(A_convh{1});
% disp(A_convh{2});
% disp('Discrete B matrices (common for both vertices):');
% disp(B_convh{1});
% disp('Discrete E matrices (disturbance matrices):');
% disp(E_convh{1});
end








