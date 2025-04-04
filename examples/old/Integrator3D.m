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
h = 1/3;

A = [1 h (h^2)/2; 0 1 h; 0 0 1];                nx = size(A,2);
B = [(h^3)/6; (h^2)/2; h];                      nu = size(B,2);
Bw = [h (h^2)/2 (h^3)/6; 0 h (h^2)/2; 0 0 h];   nw = size(Bw,2);

% Constraints
HX = [eye(nx); -eye(nx)]; hX = 5*ones(2*nx,1);%hX = [11.5;6.5;6.1; 4;6;6.5];
HU = [eye(nu); -eye(nu)]; hU = [2; 2];
HW = [eye(nw); -eye(nw)]; hW = (1/15)*ones(2*nw,1);


% Define LTI system
A_convh = {0.95*A, 1.05*A}; B_convh = {0.95*B, 1.05*B};

X = Polyhedron(HX, hX);
U = Polyhedron(HU, hU);

W_dist = Polyhedron(HW, hW);

sys = qLPV(A_convh, B_convh, Bw, X, U, W_dist);

% CCpolytope

% f = 32; % football ball
f = 7; 
% f=10; % dodecahedron

% take the polar set, having probability~=1 to have simple polytope (as we compute a simplicial one)
% F = Polyhedron(polyGenEulerAngles(nx, v),ones(v,1)).minHRep();
F = Polyhedron(numericalSpherePoints(nx, f)).normalize.dual.minHRep();
F.normalize();
% % Optionally perturb the vertices slightly from numericalSpherePoints() to avoid degeneracy:
% epsilon = 1e-6;
% vertices = vertices + epsilon * randn(size(vertices));
% vertices = vertices ./ vecnorm(vertices, 2, 2);
%%
% search for different dimensions a good template
for f=(nx+1):(nx-1):(nx+20)
    F = Polyhedron(numericalSpherePoints(nx, f,1),ones(f,1)).minHRep();
    disp("f: "+ size(F.A,1) + ", v: " + size(uniquetol(F.V,1e-10,'ByRows',true),1))
end
%%
ccPoly = CCPolytope_test(sys, F.A, true,false);


%%
plotVerticesSimplePolytope(ccPoly); % works in 2D and 3D


P_dual = ccPoly.ccpolytope.dual;
P_dual = Polyhedron(P_dual.A,ccPoly.ccpolytope.support(P_dual.A'));
hold on
% P_dual.plot("alpha",0.1)
%%
dGap = cat(2,ccPoly.dual_gap{2,:});
disp(vecnorm(dGap,inf))
[normV_s,idx_norm] = sort(vecnorm(dGap),'descend');
disp(idx_norm)
disp(normV_s)

%%
% Cost Function definition: we can either pass structs of matrices, or
% function handles defined externally. We only require that the cost
% functions are -strictly- convex, and that the function signatures are:
% RCI : @(ys,us)
% Stage/Term_cost : @(y_k/N,u_k/N,ys,us)

% Qv = blkdiag(0.01*eye(nx), 10*eye(nu));
% Qc = blkdiag(1*eye(nx), 1*eye(nu));
% RCI_cost = struct("Qv",Qv,"Qc",Qc);

RCI_cost = struct("Qy",10,"Qu",0.1);
% input cost greatly influences ROA of HTMPC_* methods

gamma = 0.98;
paramTrack = true;
if paramTrack
    Qmat = blkdiag(0.1*eye(ccPoly.f), 0.1*eye(sys.nu*ccPoly.v));
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
% htmpc = HTMPC(sys, ccPoly, costFunMan, N_ocp, gamma, sdp_opts,false);
htmpc_sp = HTMPC_SinglePolicy(sys, ccPoly, costFunMan, N_ocp, gamma, sdp_opts);

% then, specify how many (or which) vertices do you want to fix.
% H_vertices = 1:ccPoly.v; % maybe also comparing selecting vertices vs random params?
% H_vertices = [17,20,8,10];
% H_vertices = [18,16,17,20];
%%
% H_vertices = [9,13,16,6];
H_vertices = idx_norm(end-5:end); % -5
% H_vertices = idx_norm(1:4);

% H_vertices = [];
part_cctmpc = PartialCCTMPC(sys, ccPoly, H_vertices,costFunMan, N_ocp, gamma, sdp_opts);
disp("#:" + length(part_cctmpc.H_params)+ " of "+ ccPoly.f + ", params:" + mat2str(part_cctmpc.H_params'))

%% Simulation

% for debug reasons, N_mpc = k*T_orbit+1, k \in N
N_mpc = 5*N_ocp+1; % simulation steps

% (state space) system dynamics
x_sys = zeros(nx, N_mpc+1);

x_sys(:,1) = htmpc_sp.findFeasibleX0([5,-5,5]');
% x_sys(:,1) = htmpc.findFeasibleX0([5,-5,5]');
disp(x_sys(:,1)')

% x_sys(:,1) = cctmpc.findFeasibleX0([5,-5,5]');
% disp(x_sys(:,1)')
% % 
% x_sys(:,1) = part_cctmpc.findFeasibleX0([5,-5,5]');
% disp(x_sys(:,1)')


% x_sys(:,1) = [0;0;0];
u_sys = zeros(nu, N_mpc);

% (parameter space) optimal control dynamics
y0_CL = cell(1,N_mpc);
Lyap_cost = zeros(1,N_mpc);

%%
% main loop
tic
for t = 1:N_mpc
    disp("MPC iteration: " + t +"/"+ N_mpc)

    % compute input for the system
    u_sys(:,t) = cctmpc.solve(x_sys(:,t));
    % u_sys(:,t) = part_cctmpc.solve(x_sys(:,t));
    % u_sys(:,t) = htmpc.solve(x_sys(:,t));
    u_sys(:,t) = htmpc_sp.solve(x_sys(:,t));



    % propagate system dynamics
    x_sys(:,t+1) = sys.step_noDist(x_sys(:,t), u_sys(:,t));

    % save computed data
    % y0_CL{t} = cctmpc.y0_sol;
    % Lyap_cost(t) = cctmpc.Lyapunov_cost;
    % y0_CL{t} = part_cctmpc.y0_sol;
    % Lyap_cost(t) = part_cctmpc.Lyapunov_cost;
    y0_CL{t} = htmpc_sp.y0_sol; % z+\alpha*ys
    Lyap_cost(t) = htmpc_sp.Lyapunov_cost;

end
toc

%%
scriptName = "Integrator3D";
LyapunovCost;

%% Compute the feasibility region for all TMPC schemes, and compare them with
% NOTE: what to show? inner (convex hull of x_i) or outer (hyperplane) approximation?

feasRegionCCTMPC = intersect(cctmpc.computeFeasRegion(15),X); disp("CC_done.")
feasRegionPartCCTMPC = intersect(part_cctmpc.computeFeasRegion(15),X); disp("PartCC_done.")
feasRegionHTMPC = intersect(htmpc.computeFeasRegion(15),X); disp("H_done.")
feasRegionHTMPC_SP = intersect(htmpc_sp.computeFeasRegion(15),X); disp("H_SP_done.")

%%
% NOTE: to get the same result, you can alternatively use:
% MRCI = ULTISystem('A', A_convh, 'B', B_convh).invariantSet('X',X,'U',U,'D',Polyhedron(HW,hW));
MRCI = feasRegionCCTMPC;
% compare set distance of feasible regions from the MRCI
dist_CCTMPC_MRCI = setDistance(feasRegionCCTMPC, MRCI);
dist_PartCCTMPC_MRCI = setDistance(feasRegionPartCCTMPC, MRCI);
dist_HTMPC_MRCI = setDistance(feasRegionHTMPC, MRCI);
dist_HTMPC_SP_MRCI = setDistance(feasRegionHTMPC_SP, MRCI);


disp("Hausdorff distance (MRCI,CCTMPC): " + round(dist_CCTMPC_MRCI,4))
disp("Hausdorff distance (MRCI,PartCCTMPC): " + round(dist_PartCCTMPC_MRCI,4))
disp("Hausdorff distance (MRCI,HTMPC): " + round(dist_HTMPC_MRCI,4))
disp("Hausdorff distance (MRCI,HTMPC_SP): " + round(dist_HTMPC_SP_MRCI,4))
disp(" ")
disp("Is PartCCTMPC equal to CCTMPC: " + mat2str(feasRegionCCTMPC== feasRegionPartCCTMPC))
disp("Is PartCCTMPC equal to HTMPC: " + mat2str(feasRegionHTMPC== feasRegionPartCCTMPC))

disp("Volume coverage: " + round((feasRegionPartCCTMPC.volume/feasRegionCCTMPC.volume)*100,4) + "%")



%%
fig = figure("Renderer","painters","Units","centimeters");
fig.Position(3:4) = 4*[3,2.5]; % [width,height]
hold on

idx = 1:2;

X.projection(idx).plot("Alpha",0,"EdgeColor",0.4*ones(1,3),"Linewidth",1.5)
han_CC = feasRegionCCTMPC.projection(idx).plot("alpha",0.1,"EdgeColor",[0.4660 0.6740 0.1880],"Color",[0.4660 0.6740 0.1880],"Linewidth",0.75);

han_PartCC = feasRegionPartCCTMPC.projection(idx).plot("alpha",0.1,"EdgeColor",[0 0.4470 0.7410],"Color",[0 0.4470 0.7410],"Linewidth",0.75);

han_H = feasRegionHTMPC.projection(idx).plot("Alpha",0.1,"EdgeColor",[0.8500 0.3250 0.0980],"Color",[0.8500 0.3250 0.0980],"Linewidth",0.75);
%
han_SP = feasRegionHTMPC_SP.projection(idx).plot("Alpha",0.1,"EdgeColor",[0.9290 0.6940 0.1250],"Color",[0.9290 0.6940 0.1250],"Linewidth",0.75);
%
han_mRCI = Polyhedron(ccPoly.F,cctmpc.rciSol{1}).projection(idx).plot("Alpha",0.8,"Color",0.9*[0.4660 0.6740 0.1880]);

% since we have a 3D integrator, we can get the slice at x3=0 and show in
% perspective!

% % bring grid in front of everything
Ax1 = gca;
Ax1.YGrid = 'on'; Ax1.Layer = 'top'; Ax1.GridAlpha = 0.05;
set(Ax1, 'TickLabelInterpreter', 'latex');

% han_leg1 = legend(Ax1,[han_CC, han_PartCC,han_H, han_SP,han_mRPI,han_mRCI],... han_ref],...
han_leg1 = legend(Ax1,[han_CC, han_PartCC,han_H, han_SP,han_mRCI],... han_ref],...
    {'CCTMPC','PartCCTMPC','HTMPC','HTMPC\_SP',"mRCI"}, ...
    'Interpreter','latex','Location','northeast');
% title("Feasible regions and mRCIinimal Robust set comparison")


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
xlabel("x")
ylabel("y")
end


function dir_nx = getDir(nx, subDiv)
verts = polyGenEulerAngles(nx, subDiv);
dir_nx = Polyhedron(verts,ones(size(verts,1),1)).minHRep().A;
end



















