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
    addpath(genpath('C:\casadi-3.6.7-windows64-matlab2018b'))
else
    addpath(genpath('/Work_stuff/IMT_set_work/softwares/yalmip_latest'))
    addpath(genpath('/Work_stuff/IMT_set_work/Work_set_based/utilities'))
    addpath(genpath('/Work_stuff/IMT_set_work/softwares/mpt_files'))
    addpath(genpath('/Work_stuff/IMT_set_work/softwares/mosek'))
    rmpath(genpath('/Work_stuff/IMT_set_work/softwares/cplex'))
end
import casadi.*

% mptopt('zero_tol',1e-14, 'lex_tol',1e-12,'abs_tol', 1e-10)
clc
%% System definition
% Dynamics
h = 1/4;

A = [1 h (h^2)/2; 0 1 h; 0 0 1];                nx = size(A,2);
B = [(h^3)/6; (h^2)/2; h];                      nu = size(B,2);
Bw = [h (h^2)/2 (h^3)/6; 0 h (h^2)/2; 0 0 h];   nw = size(Bw,2);

% Constraints
HX = [eye(nx); -eye(nx)]; hX = 4*ones(2*nx,1);%hX = [11.5;6.5;6.1; 4;6;6.5];
HU = [eye(nu); -eye(nu)]; hU = [2; 2];
HW = [eye(nw); -eye(nw)]; hW = (1/20)*ones(2*nw,1);


% Define LTI system
A_convh = {0.95*A, 1.05*A}; B_convh = {0.95*B, 1.05*B};

X = Polyhedron(HX, hX);
U = Polyhedron(HU, hU);

W_dist = Polyhedron(HW, hW);

sys = qLPV(A_convh, B_convh, Bw, X, U, W_dist);

% CCpolytope
f = 8;%nx+1; 
F = Polyhedron(numericalSpherePoints(nx, f)).normalize.dual.minHRep();
F.normalize();
%%
ccPoly = CCPolytope_test(sys, F.A, true, false);

% ccPoly = CCPolytope_test(sys, S1.A./S1.b, false);

% ccPoly = CCPolytope_test(sys, Ys_best{11}./ys_best{11}, false);
ccPoly
%%
plotVerticesSimplePolytope(ccPoly); % works in 2D and 3D

%%
% Cost Function definition: we can either pass structs of matrices, or
% function handles defined externally. We only require that the cost
% functions are -strictly- convex, and that the function signatures are:
% RCI : @(ys,us)
% Stage/Term_cost : @(y_k/N,u_k/N,ys,us)

Qv = blkdiag(0.1*eye(nx), 0.1*eye(nu));
Qc = blkdiag(1*eye(nx), 1*eye(nu));
RCI_cost = struct("Qv",Qv,"Qc",Qc);

% RCI_cost = struct("Qy",0.1,"Qu",1);
% input cost greatly influences ROA of HTMPC_* methods

gamma = 0.95;

% as place-holders
Q_mat = blkdiag(0.1*eye(ccPoly.f), 0.1*eye(sys.nu*ccPoly.v));
P_mat = (1/(1-gamma^2))*Q_mat;


% collect all cost function definition in a single class
costFunMan = CostFunctionManager(ccPoly, Q_mat, P_mat, RCI_cost);

% equivalent cost between CCTMPC and HTMPC
[y_m,u_m,~] = compute_mRCI_set(ccPoly, costFunMan);

% Q_bar = blkdiag(eye(nx),eye(nu),1);
% P_bar = (1/(1-gamma^2))*Q_bar;
% 
% F_t_cost = [ccPoly.F, zeros(ccPoly.f,nu), y_m; 
%             zeros(ccPoly.v, nx), kron(ones(ccPoly.v,1),eye(nu)), u_m(:)];
% 
% Q_mat = F_t_cost*Q_bar*F_t_cost';
% P_mat = F_t_cost*P_bar*F_t_cost';
% 
% costFunMan = CostFunctionManager(ccPoly, Q_mat, P_mat, RCI_cost);

%%

% Build MPC scheme
sdp_opts = sdpsettings('solver','gurobi','verbose',0);


N_ocp = 3; % Note: y_N = N_ocp+1
cctmpc = CCTMPC(sys, ccPoly, costFunMan, N_ocp, gamma, sdp_opts);
htmpc_sp = HTMPC_SinglePolicy(sys, ccPoly, costFunMan, N_ocp, gamma, sdp_opts);
%% Simulation CCTMPC

% for debug reasons, N_mpc = k*T_orbit+1, k \in N
N_mpc = 4*N_ocp+1; % simulation steps
x_inits = cartesianProduct(linspace(-5,5,3),3);
x_inits(ismembertol(x_inits, zeros(1,3), 1e-5,'ByRows',true),:) = []; 

x_sys_all = cell(1,size(x_inits,1));


%%
CC_times = zeros(size(x_inits,1),N_mpc);
for i=1:size(x_inits,1)
    disp("MPC simulation: " + i +"/"+ size(x_inits,1))
    
    % (state space) system dynamics
    x_sys = zeros(nx, N_mpc+1);
    x_sys(:,1) = (1-2e-3)*cctmpc.findFeasibleX0(x_inits(i,:)');
    u_sys = zeros(nu, N_mpc);

    % (parameter space) optimal control dynamics
    y0_CL = cell(1,N_mpc);
    Lyap_cost = zeros(1,N_mpc);

    
    % main loop
    for t = 1:N_mpc
        % disp("MPC iteration: " + t +"/"+ N_mpc)
        
        tic
        % compute input for the system
        try
        u_sys(:,t) = cctmpc.solve(x_sys(:,t));
        catch
            disp(t + " point out of feasible.")
            break;
        end
        CC_times(i,t) = toc;

        % propagate system dynamics
        x_sys(:,t+1) = sys.step_noDist(x_sys(:,t), u_sys(:,t));

        % save computed data
        y0_CL{t} = cctmpc.y0_sol;
        Lyap_cost(t) = cctmpc.Lyapunov_cost;
    end
    

    x_sys_all{i} = x_sys;
end
%%
CC_times = CC_times(~any(CC_times, 2), :);
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
    x_sys(:,1) = (1-2e-3)*htmpc_sp.findFeasibleX0(x_inits(i,:)');
    u_sys = zeros(nu, N_mpc);

    % (parameter space) optimal control dynamics
    y0_CL = cell(1,N_mpc);
    Lyap_cost = zeros(1,N_mpc);
    

    % main loop
    for t = 1:N_mpc

        % disp("MPC iteration: " + t +"/"+ N_mpc)

        tic
        % compute input for the system
        try
            u_sys(:,t) = htmpc_sp.solve(x_sys(:,t));
        catch
            disp(t + " point out of feasible.")
            break;
            % u_sys(:,t) = htmpc_sp.solve((1-2e-2)*x_sys(:,t));
        end
        H_times(i,t) = toc;

        % propagate system dynamics
        x_sys(:,t+1) = sys.step_noDist(x_sys(:,t), u_sys(:,t));

        % save computed data
        y0_CL{t} = htmpc_sp.y0_sol;
        Lyap_cost(t) = htmpc_sp.Lyapunov_cost;

    end

end

H_times = H_times(~any(isnan(H_times), 2), :);
disp(min(H_times,[],"all")*1000)
disp(mean(H_times,"all")*1000)
disp(max(H_times,[],"all")*1000)

%%
scriptName = "Integrator3D";
LyapunovCost;

%% Compute the feasibility region for all TMPC schemes, and compare them with
% NOTE: what to show? inner (convex hull of x_i) or outer (hyperplane) approximation?

feasRegionCCTMPC = X.intersect(cctmpc.computeFeasRegion(20)).minVRep(); disp("CC_done.")
feasRegionHTMPC_SP = htmpc_sp.computeFeasRegion(20).minVRep(); disp("H_SP_done.")

%%
% NOTE: to get the same result, you can alternatively use:
% compare set distance of feasible regions from the MRCI
MRCI = ULTISystem('A', A_convh, 'B', B_convh, 'E',Bw...
    ).invariantSet('X',X,'U',U,'D',Polyhedron(HW,hW),'maxIterations',5);
%%
dist_CCTMPC_MRCI = setDistance(feasRegionCCTMPC, MRCI);
dist_HTMPC_SP_MRCI = setDistance(feasRegionHTMPC_SP, MRCI);

disp("Hausdorff distance (X,CCTMPC): " + round(dist_CCTMPC_MRCI,4))
disp("Hausdorff distance (X,HTMPC_SP): " + round(dist_HTMPC_SP_MRCI,4))


%%
fig = figure("Renderer","painters","Units","centimeters");
fig.Position(3:4) = 4*[3,2.5]; % [width,height]
hold on

X.plot("wire",true,"linestyle",'--',"edgecolor",0.5*ones(1,3),"edgealpha",0.5)
han_CC = feasRegionCCTMPC.plot("alpha",0.05,"EdgeColor",0.75*[0.4660 0.6740 0.1880],...
    "Color",[0.4660 0.6740 0.1880],"Linewidth",0.5,"edgealpha",0.5);

% verts_H = feasRegionHTMPC_SP.V;
% chull_H = convhull(verts_H,"Simplify",true);
% han_SP = trisurf(chull_H,verts_H(:,1),verts_H(:,1),verts_H(:,1),"FaceAlpha",0.1,"EdgeColor",0.75*[0.9290 0.6940 0.1250],...
%     "FaceColor",[0.9290 0.6940 0.1250],"Linewidth",0.5)
han_SP = feasRegionHTMPC_SP.plot("Alpha",0.05,"EdgeColor",0.75*[0.9290 0.6940 0.1250],...
    "Color",[0.9290 0.6940 0.1250],"Linewidth",0.5,"edgealpha",0.5);

han_mRCI = Polyhedron(ccPoly.F,cctmpc.rciSol{1}).plot("Alpha",0.1,"Color",[0.8500 0.3250 0.0980],...
    "EdgeColor",0.75*[0.8500 0.3250 0.0980],"Linewidth",0.5,"edgealpha",0.5);

for i=1:length(x_sys_all)
    plot3(x_sys_all{i}(1,:),x_sys_all{i}(2,:),x_sys_all{i}(3,:),"LineWidth",0.3,"Linestyle","-","Color",[0 0.4470 0.7410])
    scatter3(x_sys_all{i}(1,1),x_sys_all{i}(2,1),x_sys_all{i}(3,1), 8, 0.2*ones(1,3),"filled")
end

% % bring grid in front of everything
Ax1 = gca;
Ax1.YGrid = 'on'; Ax1.Layer = 'top'; Ax1.GridAlpha = 0.05;
set(Ax1, 'TickLabelInterpreter', 'latex');

han_leg1 = legend(Ax1,[han_CC, han_SP, han_mRCI],...
    {'$\mathcal{O}_{c}(3)$','$\mathcal{O}_h(3)$',"$P(y_m)$"}, ...
    'Interpreter','latex','Location','northeast');
han_leg1.FontSize = 12;

xlabel('$x_1$','Interpreter','latex');
ylabel('$x_2$','Interpreter','latex');
zlabel('$x_3$','Interpreter','latex');

view(-200,35)
%%
% % minimize white borders around plot
set(Ax1,'LooseInset', max(get(Ax1,'TightInset'), 0.01)) % remove border from axis
set(fig,'PaperPositionMode','Auto','PaperUnits',...
    'centimeters','PaperSize',fig.Position(3:4)) % resize pdf page

fig.Renderer = 'painters'; % a way to force saving in vector graphics
% % save the plot as PDF file
saveas(fig, '../figures/Integrator3D_Regions.pdf','pdf')



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
Polyhedron(ccPoly.F,ccPoly.y_sigma).plot("alpha",0.1,"EdgeColor",[0.9290 0.6940 0.1250],...
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



function combos = cartesianProduct(v, k)
% CARTESIANPRODUCT Generates the Cartesian product of the vector v taken k times.
%
%   combos = CARTESIANPRODUCT(v, k) returns a matrix in which each row is one
%   of the n^k possible k-tuples formed from the elements of v.
%
%   Example:
%       v = [10, 20, 30];
%       k = 3;
%       combos = cartesianProduct(v, k)
%       % Output (each row is a combination):
%       %    10    10    10
%       %    20    10    10
%       %    30    10    10
%       %    10    20    10
%       %    20    20    10
%       %    30    20    10
%       %    10    30    10
%       %       ...
%       % (Total of 27 rows since 3^3 = 27.)
%
%   This function uses ndgrid to generate k grids and then stacks the results.
%

    % Create a cell array where each cell contains the vector v.
    args = repmat({v}, 1, k);
    
    % Use ndgrid to generate k-dimensional grids.
    [G{1:k}] = ndgrid(args{:});
    
    % Determine the total number of combinations.
    numPoints = numel(G{1});
    
    % Preallocate for speed.
    combos = zeros(numPoints, k);
    
    % Flatten each grid and put it as a column.
    for i = 1:k
        combos(:, i) = G{i}(:);
    end
end

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




