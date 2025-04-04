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
clc
%% System definition
% Dynamics
Ts = 0.1;
[A,B,Bw] = defDiscreteQuadrotorDyn(Ts);
nx = size(A,2); nu = size(B,2); nw = size(Bw,2);

% [~,K,~] = idare(A,B,100*eye(nx),eye(nu));
% A = A-B*K;
global g kT
% Constraints
% states: [x y z vx vy vz pitch roll wpitch wroll]
% HX = [eye(nx); -eye(nx)]; hX = repmat([4,4,2,10,10,2,pi/4,pi/4,pi,pi]',[2,1]);
ny = 3;
HX = [[eye(ny); -eye(ny)], zeros(2*ny,nx-ny) ]; hX = repmat([4,4,2]',[2,1]);

% HX = eye(nx); HX = [HX(4:10,:);-HX(4:10,:)];
% hX = repmat([5,5,3,pi/4,pi/4,pi,pi]',[2,1]);


HU = [eye(nu); -eye(nu)]; hU = [pi/6; pi/6; 2*g-g/kT;  pi/6; pi/6; -(0-g/kT)];
HW = [eye(nw); -eye(nw)]; hW = repmat([0.01; 0.09; 0.2],[2,1]);
clear g kT;

% Define LTI system
A_convh = {A}; B_convh = {B};
% A_convh = {0.9*A, 1.1*A}; B_convh = {0.9*B, 1.1*B};

X = Polyhedron(HX, hX);
U = Polyhedron(HU, hU);

W_dist = Polyhedron(HW, hW);

sys = qLPV(A_convh, B_convh, Bw, X, U, W_dist);

% CCpolytope
F = [-eye(nx); ones(1,nx)];
f_MPI= [zeros(nx,1);1./norm(ones(1,nx))];
F = Polyhedron(F,f_MPI).minHRep();
F = (-0.1*ones(nx,1))+F;
F.normalize();


% for f=(nx+1):nx+3
%     try
%         F = Polyhedron(numericalSpherePoints(nx, f,0,1)).normalize.dual;%.minVRep();
%         F = F.normalize();
%     catch
%         continue
%     end
%     disp("f: "+ size(F.A,1) + ", v: " + size(uniquetol(F.V,1e-12,'ByRows',true),1))
% end

%%
ccPoly = CCPolytope_test(sys, F.A, true, false);

% ccPoly = CCPolytope_test(sys, ccPoly.F, true, false);


%%
% RCI_cost = struct("Qy",10,"Qu",0.1);

costFunMan = CostFunctionManager(ccPoly, Qmat, Pmat, RCI_cost);
[ys,~,~] = compute_mRCI_set(ccPoly, costFunMan);
RCI = Polyhedron(ccPoly.F,ys);
V_s = RCI.V';
RCI_proj = RCI.projection(1:3).minVRep().minHRep();
V_s_proj = RCI_proj.V';
normV_s = vecnorm(V_s_proj,inf);
disp(normV_s)

% [normV_s,idx_norm] = sort(normV_s,'descend');
% disp(normV_s(1:5))
figure; RCI_proj.plot("alpha",0.1); axis equal
RCI_proj.volume



%%
dGap = cat(2,ccPoly.dual_gap{:});
% disp(vecnorm(dGap,inf))
[normV_s,idx_norm] = sort(vecnorm(dGap),'descend');
% [normV_s,idx_norm] = sort(vecnorm(RCI.V(:,1:2)',inf),'descend');
disp(normV_s)

%%
% select the facet opposite to the vertex
F_new = ccPoly.F;
fb = ccPoly.ccpolytope.b;

for v=idx_norm(1)
    % this opposite facet retrieving works only for simplices
    % F_opp = ccPoly.F(~any(ccPoly.Vi_s{v},1),:);
    F_dual = ccPoly.ccpolytope.dual.A(v,:);
    F_tan =  mean(ccPoly.F(any(ccPoly.Vi_s{v},1),:),1);

    % epsi = (1-0.2)*ccPoly.ccpolytope.support(-F_opp');
    epsi = (1-0.001)*ccPoly.ccpolytope.support(F_dual');
    % epsi = (1-0.01)*ccPoly.ccpolytope.support((F_tan)');


    % F_new = [F_new; -F_opp./epsi];
    F_new = [F_new; F_dual./epsi];
    % F_new = [F_new; F_tan./epsi];


    fb = [fb; 1];
end
F_cut = Polyhedron(F_new,fb).minHRep().minVRep();
F_cut.normalize()
%%
ccPoly = CCPolytope_test(sys, F_cut.A, true, false);

%%

%%
% Cost Function definition: we can either pwass structs of matrices, or
% function handles defined externally. We only require that the cost
% functions are -strictly- convex, and that the function signatures are:
% RCI : @(ys,us)
% Stage/Term_cost : @(y_k/N,u_k/N,ys,us)

Qv = blkdiag( diag([100*ones(1,3), 0.1*ones(1,7)]), 1*eye(nu));
Qc = blkdiag(1*eye(nx), 1*eye(nu));
RCI_cost = struct("Qv",Qv,"Qc",Qc);

% RCI_cost = struct("Qy",10,"Qu",0.1);
% input cost greatly influences ROA of HTMPC_* methods

gamma = 0.995;

paramTrack = true;
if paramTrack
    Qmat = blkdiag(0.1*eye(ccPoly.f), 0.1*eye(sys.nu*ccPoly.v));
    Pmat = (1/(1-gamma^2))*Qmat;
else
    [Qmat, Pmat] = defineCustomTrackCost(sys, ccPoly, Qv, gamma);
end

% collect all cost function definition in a single class
costFunMan = CostFunctionManager(ccPoly, Qmat, Pmat, RCI_cost);
%%

% Build MPC scheme
sdp_opts = sdpsettings('solver','gurobi','verbose',0);%,'savesolverinput',1);

N_ocp = 3; % Note: y_N = N_ocp+1
cctmpc = CCTMPC(sys, ccPoly, costFunMan, N_ocp, gamma, sdp_opts);
% htmpc = HTMPC(sys, ccPoly, costFunMan, N_ocp, gamma, sdp_opts, false);
% htmpc_sp = HTMPC_SinglePolicy(sys, ccPoly, costFunMan, N_ocp, gamma, sdp_opts);

%%
mRCI_vs = Polyhedron(ccPoly.F,cctmpc.rciSol{1}).V';
[norm_mRCI_vs,idx_norm_mRCI] = sort(vecnorm(mRCI_vs),'descend');
disp(norm_mRCI_vs(1:5))
%%
H_vertices = [];%[idx_norm_mRCI(end)];
part_cctmpc = PartialCCTMPC(sys, ccPoly, H_vertices,costFunMan, N_ocp, gamma, sdp_opts);
disp("#:" + length(part_cctmpc.H_params) + ", params:" + mat2str(part_cctmpc.H_params'))
%% Simulation

% for debug reasons, N_mpc = k*T_orbit+1, k \in N
N_mpc = 50*N_ocp+1; % simulation steps

% (state space) system dynamics
x_sys = zeros(nx, N_mpc+1);
x_sys(1:3,1) = [-3.3;-2.4;-0.35]; % initial condition
% x_sys(1:3,1) = [-3.6;-2.9;-0.35]; % initial condition
% x_sys(1:3,1) = [-1;-1;0]; % initial condition

x_sys(:,1) = cctmpc.findFeasibleX0([5,5,-5,zeros(1,7)]',diag([0.001*ones(1,3),100*ones(1,7)]));
disp(x_sys(1:3,1)')
% x_sys(:,1) = part_cctmpc.findFeasibleX0([5,5,-5,zeros(1,7)]',diag([0.001*ones(1,3),100*ones(1,7)]));
% disp(x_sys(1:3,1)')
% x_sys(:,1) = htmpc.findFeasibleX0([5,5,-5,zeros(1,7)]',diag([0.001*ones(1,3),100*ones(1,7)]));
% disp(x_sys(1:3,1)')
% x_sys(:,1) = htmpc_sp.findFeasibleX0([5,5,-5,zeros(1,7)]',diag([0.001*ones(1,3),100*ones(1,7)]));
% disp(x_sys(1:3,1)')
%%
x_sys(1:3,1) = [0;0;-0.3]; % initial condition
% x_sys(1:3,1) = [-5;-5;-5]; % initial condition

%%
u_sys = zeros(nu, N_mpc);

% (parameter space) optimal control dynamics
y0_CL = cell(1,N_mpc);

Lyap_cost = zeros(1,N_mpc);

tic
% main loop
for t = 1:N_mpc
    disp("MPC iteration: " + t +"/"+ N_mpc)
    % compute input for the system
    u_sys(:,t) = cctmpc.solve(x_sys(:,t));
    % u_sys(:,t) = part_cctmpc.solve(x_sys(:,t));
    % u_sys(:,t) = htmpc.solve(x_sys(:,t));
    % u_sys(:,t) = htmpc_sp.solve(x_sys(:,t));


    % propagate system dynamics
    x_sys(:,t+1) = sys.step(x_sys(:,t), u_sys(:,t));

    % save computed data
    y0_CL{t} = cctmpc.y0_sol;
    Lyap_cost(t) = cctmpc.Lyapunov_cost;
    % y0_CL{t} = part_cctmpc.y0_sol;
    % Lyap_cost(t) = part_cctmpc.Lyapunov_cost;
    % y0_CL{t} = htmpc.y0_sol;
    % Lyap_cost(t) = htmpc.Lyapunov_cost;
    % y0_CL{t} = htmpc_sp.y0_sol;
    % Lyap_cost(t) = htmpc_sp.Lyapunov_cost;

end
toc

%%
scriptName = "QuadrotorCCTMPC";
LyapunovCost;


%% Compute the feasibility regions for all TMPC schemes
% since nx > 3, we show only relevant components
X_space = X.projection(1:3);
X_vel = X.projection(4:6);
X_angles = X.projection(7:8);
X_angRate = X.projection(9:10);
%%
dir2D = getDir(2, 100); dir3D = getDir(3, 15);

% %%
templ = zeros(size(dir3D,1),nx); templ(:,1:3) = dir3D;
RoA_CCTMPC_space = cctmpc.computeFeasRegion(templ).projection(1:3);

% RoA_HTMPC_space = htmpc.computeFeasRegion(templ).projection(1:3);
%%
idx = 1:3;
% templ = zeros(size(dir3D,1),nx); templ(:,idx) = dir3D;
% RoA_CCTMPC_idx = cctmpc.computeFeasRegion(templ).projection(idx); disp("CC_done.")

fig = figure("Renderer","painters","Units","centimeters");
fig.Position(3:4) = 4*[3,2.5]; % [width,height]
hold on
han_X = X_space.plot("Alpha",0,"EdgeColor",0.4*ones(1,3),"Linewidth",1);

% RoA_CCTMPC_idx.plot("alpha",0.02,"EdgeColor",[0.4660 0.6740 0.1880],"Color",[0.4660 0.6740 0.1880],"Linewidth",1.25);

xlabel('$x(t)$','Interpreter','latex');
ylabel('$y(t)$','Interpreter','latex');
zlabel('$z(t)$','Interpreter','latex');

% project the trajectory dynamics
xt_idx = x_sys(idx,:);
plot3(xt_idx(1,:),xt_idx(2,:),xt_idx(3,:),"Color",[0.8500 0.3250 0.0980],"LineWidth",1)

% view(az,el);

maxInd = find(Lyap_cost <= 5e-4);
% maxInd = N_mpc;
for t=1:floor(maxInd/8):maxInd
    disp(t)
    Poly_y0_space = Polyhedron(ccPoly.F,y0_CL{t}).projection(idx);
    if size(Poly_y0_space.V,1) == 1
    scatter3(xt_idx(1,t),xt_idx(2,t),xt_idx(3,t),100,[0 0.4470 0.7410],"Marker",".")
    elseif Poly_y0_space.volume < 1e-10
    scatter3(xt_idx(1,t),xt_idx(2,t),xt_idx(3,t),100,[0 0.4470 0.7410],"Marker",".")
    else
    han_CC = Poly_y0_space.plot("alpha",0.05,"EdgeColor",[0 0.4470 0.7410],"Color",[0 0.4470 0.7410],"Linewidth",0.5);
    end
end
han_CC = Polyhedron(ccPoly.F,y0_CL{maxInd(1)}).projection(idx).plot("alpha",0.05,"EdgeColor",[0 0.4470 0.7410],"Color",[0 0.4470 0.7410],"Linewidth",0.5);

han_mRCI = Polyhedron(ccPoly.F,cctmpc.rciSol{1}).projection(idx).minVRep().plot("Alpha",0.2,"Color",[0.9290 0.6940 0.1250],"Linewidth",0.5);


Ax1 = gca;
Ax1.YGrid = 'on'; Ax1.Layer = 'top'; Ax1.GridAlpha = 0.05;
set(Ax1, 'TickLabelInterpreter', 'latex');

% han_leg1 = legend(Ax1,[han_CC, han_PartCC,han_H, han_SP,han_mRPI,han_mRCI],... han_ref],...
han_leg1 = legend(Ax1,[han_CC, han_mRCI, han_X],... han_ref],...
    {'$\mathrm{proj}_{[x,y,z]}(X(y_0^*))$','$\mathrm{mRCI}$','$\mathrm{proj}_{[x,y,z]}(\mathcal{X})$'}, ...
    'Interpreter','latex','Location','northeast');


% % minimize white borders around plot
set(Ax1,'LooseInset', max(get(Ax1,'TightInset'), 0.01)) % remove border from axis
set(fig,'PaperPositionMode','Auto','PaperUnits',...
    'centimeters','PaperSize',fig.Position(3:4)) % resize pdf page

fig.Renderer = 'painters'; % a way to force saving in vector graphics
% % save the plot as PDF file
saveas(fig, '../figures/QuadrotorCCTMPC_space.pdf','pdf')

%%
idx = 4:6;
fig = figure("Renderer","painters","Units","centimeters");
fig.Position(3:4) = 4*[3,2.5]; % [width,height]
hold on
X_vel.plot("Alpha",0,"EdgeColor",0.4*ones(1,3),"Linewidth",1.5)

templ = zeros(size(dir3D,1),nx); templ(:,idx) = dir3D;
% RoA_CCTMPC_idx = cctmpc.computeFeasRegion(templ).projection(idx); disp("CC_done.")
% RoA_CCTMPC_idx.plot("alpha",0.02,"EdgeColor",[0.4660 0.6740 0.1880],"Color",[0.4660 0.6740 0.1880],"Linewidth",1.25);

xlabel('X_{dot}'); ylabel('Y_{dot}'); zlabel('Z_{dot}');

% project the trajectory dynamics
xt_idx = x_sys(idx,:);
plot3(xt_idx(1,:),xt_idx(2,:),xt_idx(3,:),"Color","r","Marker",'o')

% for t=1:N_mpc
%     disp(t)
%     Poly_y0_space = Polyhedron(ccPoly.F,OCP_y{t}).projection(idx);
%     Poly_y0_space.plot("alpha",0.02,"EdgeColor",[0 0.4470 0.7410],"Color",[0 0.4470 0.7410],"Linewidth",0.75);
% end

Polyhedron(ccPoly.F,cctmpc.rciSol{1}).projection(idx).plot("Alpha",0.1,"Color",0.9*[0.8500 0.3250 0.0980]);


%%
idx = 7:8;
fig = figure("Renderer","painters","Units","centimeters");
fig.Position(3:4) = 4*[3,2.5]; % [width,height]
hold on
X_angles.plot("Alpha",0,"EdgeColor",0.4*ones(1,3),"Linewidth",1.5)

templ = zeros(size(dir2D,1),nx); templ(:,idx) = dir2D;
RoA_CCTMPC_idx = cctmpc.computeFeasRegion(templ).projection(idx); disp("CC_done.")
RoA_CCTMPC_idx.plot("alpha",0.02,"EdgeColor",[0.4660 0.6740 0.1880],"Color",[0.4660 0.6740 0.1880],"Linewidth",1.25);

xlabel('Pitch'); ylabel('Roll')

% project the trajectory dynamics
xt_idx = x_sys(idx,:);
plot(xt_idx(1,:),xt_idx(2,:),"Color","r","Marker",'o')

for t=1:N_mpc
    disp(t)
    Poly_y0_space = Polyhedron(ccPoly.F,y0_CL{t}).projection(idx);
    Poly_y0_space.plot("alpha",0.02,"EdgeColor",[0 0.4470 0.7410],"Color",[0 0.4470 0.7410],"Linewidth",0.75);
end

Polyhedron(ccPoly.F,cctmpc.rciSol{1}).projection(idx).plot("Alpha",0.1,"Color",0.9*[0.8500 0.3250 0.0980]);


%%
idx = 9:10;
fig = figure("Renderer","painters","Units","centimeters");
fig.Position(3:4) = 4*[3,2.5]; % [width,height]
hold on
X_angRate.plot("Alpha",0,"EdgeColor",0.4*ones(1,3),"Linewidth",1.5)

templ = zeros(size(dir2D,1),nx); templ(:,idx) = dir2D;
RoA_CCTMPC_idx = cctmpc.computeFeasRegion(templ).projection(idx); disp("CC_done.")
RoA_CCTMPC_idx.plot("alpha",0.02,"EdgeColor",[0.4660 0.6740 0.1880],"Color",[0.4660 0.6740 0.1880],"Linewidth",1.25);

xlabel('Pitch_{dot}'); ylabel('Roll_{dot}')

% project the trajectory dynamics
xt_idx = x_sys(idx,:);
plot(xt_idx(1,:),xt_idx(2,:),"Color","r","Marker",'o')

for t=1:N_mpc
    Poly_y0_space = Polyhedron(ccPoly.F,y0_CL{t}).projection(idx);
    Poly_y0_space.plot("alpha",0.02,"EdgeColor",[0 0.4470 0.7410],"Color",[0 0.4470 0.7410],"Linewidth",0.75);
end

Polyhedron(ccPoly.F,cctmpc.rciSol{1}).projection(idx).plot("Alpha",0.1,"Color",0.9*[0.8500 0.3250 0.0980]);




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

% Q = blkdiag(Qv_yp,Qv_up)+1e-2*eye(f+v*nu,f+v*nu);
Q = 1e-1*blkdiag(Qv_yp,Qv_up) + ...
    1*blkdiag(0.1*eye(f), 0.1*eye(nu*v));

P = (1/(1-gamma^2))*Q;

% define function handles (here quadratic costs)
Qtrack_han = @(y_k,u_k,ys,us) weighted2NormSquared([y_k-ys; u_k(:)-us(:)],Q);
Ptrack_han = @(y_N,u_N,ys,us) weighted2NormSquared([y_N-ys; u_N(:)-us(:)],P);

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

for v=1:ccPoly.v
    x = ccPoly.Vi_s{v}*ones(ccPoly.f,1);
    x = num2cell(x);
    text(x{:}, sprintf('$V_{%d} y$', v), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'FontSize', 12,"Interpreter","latex");
end
end



function dir_nx = getDir(nx, subDiv)
verts = polyGenEulerAngles(nx, subDiv);
dir_nx = Polyhedron(verts,ones(size(verts,1),1)).minHRep().A;
end





function [A_d, B_d, Bw_d] = defDiscreteQuadrotorDyn(Ts)

% get Continuous time Quadrotor model
[A_num, B_num, Bw] = getQuadrotorModel();

% discretize it using Euler Approximation
% (we don't need a very accurate model, only a test one in discrete time)
A_d = eye(size(A_num,1)) + Ts*A_num;
B_d = Ts*B_num;
Bw_d = Ts*Bw;

end

function [A_num, B_num, Bw] = getQuadrotorModel()
global g d0 kT n0
g  = 9.81; % gravitational constant
d0 = 10; d1 = 8; % coefficients
kT = 0.91; %
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


% %%
% % for debug reasons, N_mpc = k*T_orbit+1, k \in N
% N_mpc = 10*N_ocp+1; % simulation steps
%
% % (state space) system dynamics
% x_sys = zeros(nx, N_mpc+1); x_sys(1:3,1) = [1;0;0]; % initial condition
% u_sys = zeros(nu, N_mpc);
%
% % (parameter space) optimal control dynamics
% OCP_y = cell(1,N_mpc); OCP_u = cell(1,N_mpc);
%
% Lyap_cost = zeros(1,N_mpc);
%
% % main loop
% for t = 1:N_mpc
%     disp("MPC iteration: " + t +"/"+ N_mpc)
%     % compute input for the system
%     u_sys(:,t) = htmpc.solve(x_sys(:,t));
%
%     % propagate system dynamics
%     x_sys(:,t+1) = sys.step(x_sys(:,t), u_sys(:,t));
%
%     % save computed data
%     [OCP_y{t}, OCP_u{t}] = htmpc.ocpSol{1:2};
%     Lyap_cost(t) = htmpc.Lyapunov_cost;
% end
%
% %%
% scriptName = "Quadrotor_HTMPC";
% LyapunovCost;
%
%
% %%
% % for debug reasons, N_mpc = k*T_orbit+1, k \in N
% N_mpc = 10*N_ocp+1; % simulation steps
%
% % (state space) system dynamics
% x_sys = zeros(nx, N_mpc+1); x_sys(1:3,1) = [1;0;0]; % initial condition
% u_sys = zeros(nu, N_mpc);
%
% % (parameter space) optimal control dynamics
% OCP_y = cell(1,N_mpc); OCP_u = cell(1,N_mpc);
%
% Lyap_cost = zeros(1,N_mpc);
%
% % main loop
% for t = 1:N_mpc
%     disp("MPC iteration: " + t +"/"+ N_mpc)
%     % compute input for the system
%     u_sys(:,t) = part_cctmpc.solve(x_sys(:,t));
%
%     % propagate system dynamics
%     x_sys(:,t+1) = sys.step(x_sys(:,t), u_sys(:,t));
%
%     % save computed data
%     [OCP_y{t}, OCP_u{t}] = part_cctmpc.ocpSol{1:2};
%     Lyap_cost(t) = part_cctmpc.Lyapunov_cost;
% end
%
% %%
% scriptName = "Quadrotor_PartCCTMPC";
% LyapunovCost;
%
function [ys,us,RCI_setcost] = compute_mRCI_set(ccPoly, costFunMan)

% compute the mRCI for the given template
y_m = sdpvar(ccPoly.f,1);
u_m = sdpvar(ccPoly.sys.nu, ccPoly.v,'full');

% cost function
cost = costFunMan.cost_RCI(y_m,u_m);

% constraints
constr = constrSetS(ccPoly, y_m,u_m,y_m);

% solve the OCP
optimize(constr,cost, sdpsettings('solver','gurobi','verbose',5));

ys = value(y_m); us = value(u_m);
RCI_setcost = value(cost);
end

function constr = constrSetS(ccPoly,y,u,yp)
% define the Configuration Constrained RFIT set.
% NOTE: if the system convex hulls are updated, the constraints
% do it accordingly
constr = [];
for j=1:ccPoly.v
    for i=1:ccPoly.sys.nm % #models
        constr = [constr;
            ccPoly.F*(ccPoly.sys.A_convh{i}*ccPoly.Vi_s{j}*y + ccPoly.sys.B_convh{i}*u(:,j) )+ ccPoly.d <= yp];
    end
    constr = [constr;   ccPoly.sys.X.A*ccPoly.Vi_s{j}*y <= ccPoly.sys.X.b;    ccPoly.sys.U.A*u(:,j) <= ccPoly.sys.U.b];

end
constr = [constr; ccPoly.E*y <= 0];
end
