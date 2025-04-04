addpath(genpath('../src'))

if strcmp(getenv('USER'),'filippob')
    % -----------------------
    % NOTE: we require that those softwares are already installed on the system
    addpath(genpath('~/workSW/gurobi1201/linux64/matlab/'))
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
% Dynamics
A = [1 1; 0 1];     nx = size(A,2);
B = [0; 1];         nu = size(B,2);
% Bw = ones(nx,1);    nw =  size(Bw,2);
Bw = eye(nx);    nw =  size(Bw,2);


% Constraints
HX = [eye(nx); -eye(nx)]; hX = [5; 3; 5; 2.5];
HU = [eye(nu); -eye(nu)]; hU = [2;1];
HW = [eye(nw); -eye(nw)]; hW = [0.1; 0.1; 0.1; 0.1];

% Define qLPV system
zeta = 0.15; % mult. uncertainty
A1 = A-zeta; A1(2,1) = 0; A2 = A+zeta; A2(2,1) = 0;
B1 = B-zeta; B1(1,1) = 0; B2 = B+zeta; B2(1,1) = 0;

% A_convh = {A}; B_convh = {B};
A_convh = {A1, A2};
B_convh = {B1, B2};

% Acurr_han = @(~,~,~) A; Bcurr_han = @(~,~,~) B;

X = Polyhedron(HX, hX);
U = Polyhedron(HU, hU);

W_dist = Polyhedron(HW, hW);

sys = qLPV(A_convh, B_convh, Bw, X, U, W_dist);%, Acurr_han, Bcurr_han);

% By using polar coordinates for parameterization, we ensure to get a
% CCPolytope with minimal representation, i.e. f = v
f = 21;
% F = Polyhedron(polyGenEulerAngles(nx, f),ones(f,1)).minHRep();
F = Polyhedron(numericalSpherePoints(nx, f)).normalize.minHRep().dual();
%%
ccPoly = CCPolytope_test(sys, F.A, true,false); % ( sys, F[, computeCC_RCIset] )

%%

%%
plotVerticesSimplePolytope(ccPoly); % works in 2D and 3D


% P_dual = ccPoly.ccpolytope.dual;
% P_dual = Polyhedron(P_dual.A,ccPoly.ccpolytope.support(P_dual.A'));
% hold on
% P_dual.plot("alpha",0.1)
%%
dGap = cat(2,ccPoly.dual_gap{2,:});
disp(vecnorm(dGap,inf))
[normV_s,idx_norm] = sort(vecnorm(dGap),'descend');
disp(idx_norm)
disp(normV_s)

%% select the facet opposite to the vertex
plotVerticesSimplePolytope(ccPoly); % works in 2D and 3D

F_new = ccPoly.F;
fb = ccPoly.ccpolytope.b;

for v=idx_norm(end)

F_tan =  mean(ccPoly.F(any(ccPoly.Vi_s{v},1),:),1);
% ccPoly.F 

% epsi = (1-0.01)*ccPoly.ccpolytope.support(-F_opp');
% epsi = (1-0.1)*ccPoly.ccpolytope.support(F_dual');
epsi = (1-0.2)*ccPoly.ccpolytope.support((F_tan)');


% F_new = [F_new; -F_opp];
% F_new = [F_new; F_dual./epsi];
F_new = [F_new; F_tan./epsi];
fb = [fb; 1];

end

F_cut = Polyhedron(F_new,fb).minHRep().minVRep();
F_cut.plot("alpha",0.2)
F_cut
%%
ccPoly = CCPolytope_test(sys, F_cut.A, true, false); % ( sys, F[, computeCC_RCIset] )
plotVerticesSimplePolytope(ccPoly); % works in 2D and 3D


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

RCI_cost = struct("Qy",10,"Qu",0.1);
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
sdp_opts = sdpsettings('solver','gurobi','verbose',0);

N_ocp = 3; % Note: y_N = N_ocp+1
cctmpc = CCTMPC(sys, ccPoly, costFunMan, N_ocp, gamma, sdp_opts);
htmpc = HTMPC(sys, ccPoly, costFunMan, N_ocp, gamma, sdp_opts, false);
htmpc_sp = HTMPC_SinglePolicy(sys, ccPoly, costFunMan, N_ocp, gamma, sdp_opts);

% H_vertices = 1:ccPoly.v; % maybe also comparing selecting vertices vs random params?
% H_vertices = [];
% H_vertices = [1,6];
% % H_vertices = 1:2:6;
H_vertices = idx_norm(end-2:end); % -5
H_vertices = idx_norm(1:3);

part_cctmpc = PartialCCTMPC(sys, ccPoly, H_vertices,costFunMan, N_ocp, gamma, sdp_opts, ocp_opts);
disp("#:" + length(part_cctmpc.H_params)+ " of "+ ccPoly.f + ", params:" + mat2str(part_cctmpc.H_params'))


%% Simulation
% 
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
feasRegionCCTMPC = cctmpc.computeFeasRegion(100); disp("CC_done.")
feasRegionHTMPC = htmpc.computeFeasRegion(100); disp("H_done.")
feasRegionHTMPC_SP = htmpc_sp.computeFeasRegion(100); disp("H_SP_done.")
feasRegionPartCCTMPC = part_cctmpc.computeFeasRegion(100); disp("PartCC_done.")

% NOTE: to get the same result, you can alternatively use:
MRCI = ULTISystem('A', A_convh, 'B', B_convh, 'E',Bw).invariantSet('X',X,'U',U,'D',Polyhedron(HW,hW));

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
han_CC = feasRegionCCTMPC.plot("Alpha",0.1,"EdgeColor",[0.4660 0.6740 0.1880],"Color",[0.4660 0.6740 0.1880],"Linewidth",1.25);

feasRegionPartCCTMPC.plot("LineStyle",'none',"Color",ones(3,1));
han_PartCC = feasRegionPartCCTMPC.plot("Alpha",0.1,"EdgeColor",[0 0.4470 0.7410],"Color",[0 0.4470 0.7410],"Linewidth",1.25);


feasRegionHTMPC.plot("LineStyle",'none',"Color",ones(3,1));
han_H = feasRegionHTMPC.plot("Alpha",0.1,"EdgeColor",[0.8500 0.3250 0.0980],"Color",[0.8500 0.3250 0.0980],"Linewidth",1.25);

feasRegionHTMPC_SP.plot("LineStyle",'none',"Color",ones(3,1));
han_SP = feasRegionHTMPC_SP.plot("Alpha",0.1,"EdgeColor",[0.9290 0.6940 0.1250],"Color",[0.9290 0.6940 0.1250],"Linewidth",1.25);

han_mRCI = Polyhedron(ccPoly.F,cctmpc.rciSol{1}).plot("Alpha",1,"Color",0.9*[0.4660 0.6740 0.1880]);


% % bring grid in front of everything
Ax1 = gca;
Ax1.YGrid = 'on'; Ax1.Layer = 'top'; Ax1.GridAlpha = 0.05;
set(Ax1, 'TickLabelInterpreter', 'latex');

han_leg1 = legend(Ax1,[han_CC, han_PartCC, han_H, han_SP, han_mRCI],... han_ref],...
    {'CCTMPC','PartCCTMPC','HTMPC','HTMPC\_SP','mRCI'}, ...
    'Interpreter','latex','Location','northeast');
title("Feasible regions and minimal Robust sets comparison")


% % minimize white borders around plot
set(Ax1,'LooseInset', max(get(Ax1,'TightInset'), 0.01)) % remove border from axis
set(fig,'PaperPositionMode','Auto','PaperUnits',...
    'centimeters','PaperSize',fig.Position(3:4)) % resize pdf page

fig.Renderer = 'painters'; % a way to force saving in vector graphics
% % save the plot as PDF file
saveas(fig, '../figures/Integrator2D_Regions_allVerts.pdf','pdf')

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
ccPoly.ccpolytope.plot("alpha",0.1,"EdgeColor",[0.9290 0.6940 0.1250],...
    "Color",[0.9290 0.6940 0.1250],"Linewidth",0.5,"LineStyle",':');
ccPoly.sys.X.plot("alpha",0,"EdgeColor",0.3*ones(1,3));

for v=1:ccPoly.v
    x = ccPoly.Vi_s{v}*ones(ccPoly.f,1);
    x = num2cell(x);
    text(x{:}, sprintf('$V_{%d} y$', v), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'FontSize', 12,"Interpreter","latex");
end
end









