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
Ts = 0.1;
[A,B,Bw] = defDiscreteQuadrotorDyn(Ts);
nx = size(A,2); nu = size(B,2); nw = size(Bw,2);

global g kT
% Constraints
% states: [x y z vx vy vz pitch roll wpitch wroll]
HX = [eye(nx); -eye(nx)]; hX = repmat([4,4,2,10,10,5,pi/3,pi/3,pi,pi]',[2,1]);

% HX = eye(nx); HX = [HX(4:10,:);-HX(4:10,:)];
% hX = repmat([5,5,3,pi/4,pi/4,pi,pi]',[2,1]);

HU = [eye(nu); -eye(nu)]; hU = [pi/4; pi/4; 2*g-g/kT;  pi/4; pi/4; -(0-g/kT)];
HW = [eye(nw); -eye(nw)]; hW = repmat([0.05; 0.05; 0.1],[2,1]);
clear g kT;

% Define LTI system
A_convh = {A}; B_convh = {B};

X = Polyhedron(HX, hX);
U = Polyhedron(HU, hU);

W_dist = Polyhedron(HW, hW);

sys = qLPV(A_convh, B_convh, Bw, X, U, W_dist);

%% generate a template polytope
F_tilde = [-eye(nx); ones(1,nx)];
f=size(F_tilde,1);
f_MPI = ones(f,1);

% ensure we have a first feasible (F,E,V) triple
ccPoly = CCPolytope(sys, F_tilde);


%% Configuration Triple Refinement (CTR) algorithm
opti = casadi.Opti('conic');
i_max = 3;

save_Triple = cell(1,i_max+1); save_yM = cell(1,i_max+1); save_sigma = inf(1,i_max+1);

% first Configuration Triple
F = ccPoly.F; E = full(ccPoly.E); V = ccPoly.Vi_s;
f = size(F,1); e = size(E,1); v = length(V);

% generate the OCP problem
hausFun = genOCPFun(opti,sys,f,e,v); %{'F','E','V','d','yM_0'}
d = sys.W_dist.support((F*sys.Bw)');

[sigma,yM] = hausFun(F,E,V{:},d, f_MPI);
yM = full(yM); sigma = full(sigma);

disp(sigma)
save_Triple{1} = {F,E,V};  save_yM{1} = yM;  save_sigma(1) = sigma;

y_cut_s = cell(1,i_max);
tic
for i=1:i_max
    % get vertices in state space.
    verts = cellfun(@(v) v*yM, V,"UniformOutput",false);
    v_num = length(verts);

    sigma_s = zeros(1,v_num); yM_s = cell(1,v_num);
    new_F_s = cell(1,v_num); new_E_s = cell(1,v_num); new_V_s = cell(1,v_num);

    new_d_s = cell(1,v_num);
    % generate the OCP problem (we already know the complexity at each iteration)
    f_new = size(F,1)+1; v_new = length(V)+nx-1; e_new = (f_new)*v_new;
    hausFun = genOCPFun(opti,sys,f_new,e_new,v_new);

    for j=1:v_num % parfor can be used here
        disp("Candidate: " + j + "/" + v_num)
        % select a vertex
        vert = verts{j};

        % compute a new facet (e.g. average of facets generating the vertex)
        F_cut = mean(F(any(V{j},1),:),1);

        % get the two highest support functions to find a partial ordering
        [supp_vs,~] = maxk(cellfun(@(v) F_cut*v,verts),2);

        % use \kappa\in(0,1) to ensure the truncation is shallow
        y_cut = mean([supp_vs(1), max(0,supp_vs(2))]); % kappa = 1/2
        y_cut_s{i} = [y_cut_s{i},y_cut];
        % compute the new Triple
        new_F_s{j} = [F; F_cut];
        new_V_s{j} = updateV(V, j, new_F_s{j});
        new_E_s{j} = updateE(new_F_s{j}, new_V_s{j}, E);

        % update the support function
        new_d_s{j} = sys.W_dist.support((new_F_s{j}*sys.Bw)');
        % compute updated CC-RCI set
        [sigma_j,yM_j] = hausFun(new_F_s{j},new_E_s{j},new_V_s{j}{:},new_d_s{j}, [yM;y_cut]);
        yM_s{j} = full(yM_j); sigma_s(j) = full(sigma_j);
    end
    % save best resulting Triple
    [sigma,best_idx] = min(sigma_s);
    disp(sigma)
    F = new_F_s{best_idx}; E = new_E_s{best_idx}; V = new_V_s{best_idx};
    yM = yM_s{best_idx};

    save_Triple{i+1} = {F,E,V}; save_yM{i+1} = yM; save_sigma(i+1) = sigma;
end
toc

%% Plot section
plot_RCICut;
plot_RCIOneIter;
plot_RCIAfterIterations;
plot_sigmaOverIterations;


%% --- support functions ----
function new_V = updateV(V, j, new_F)
F_old = new_F(1:end-1,:); F_cut = new_F(end,:);
nx = size(F_cut,2); f = size(V{1},2);

% Isolate matrix to be removed.
V_tmp = V; V_tmp(j) = [];

% increase the size of each V_j by one empty column
V_tmp = cellfun(@(v_i) [v_i, zeros(nx,1)],V_tmp,"UniformOutput",false);

% get the facets defining the vertex
facets_idx = find(any(V{j},1));

% compute the combination of the facets for each new vertices
combFacets = nchoosek(facets_idx,nx-1);
for i=1:nx %as we're adding an (nx-1)-simplex
    % TODO: remove matrix inversion by using Sherman–Morrison formula

    % take the F_active for the vertex
    F_act = inv(V{j}(:,facets_idx));
    % substitute a single row each time
    F_new = [F_old(combFacets(i,:),:); F_cut];

    % append the new offset-vertex map
    one_map = zeros(nx,f+1); one_map(:,[combFacets(i,:),f+1]) = eye(nx);
    V_tmp{end+1} = F_new \ one_map;
end

new_V = V_tmp;
end


function E = updateE(F, V, old_E)
% TODO: update the conic constraint old_E, instead of recomputing from scratch.
old_E;
E = sparse(0,0);
for j=1:length(V)
    E = [E; F*V{j}-speye(size(F,1))];
end

% % NOTE: not removing redundant rows is often faster on low dimensions.
% E = sparse(E);
% E_min = Polyhedron(E,sparse(size(E,1),1)).minHRep;
% E = sparse(E_min.A);
end


function hausFun = genOCPFun(opti,sys,f,e,v)
opti.subject_to(); % clean the previous formulation

xV_s = sys.X.V'; numVX = size(xV_s,2);

% parameters of our problem
F = opti.parameter(f,sys.nx);
E = opti.parameter(e,f);
V = arrayfun(@(j) opti.parameter(sys.nx,f), 1:v, "UniformOutput", false);

% compute disturbance support for the given template
d = opti.parameter(f,1);

% optimization variables
yM = opti.variable(f,1);
u = opti.variable(sys.nu, v);
z = opti.variable(sys.nx, numVX);

% cost definition
cost = 0;
for k=1:numVX
    cost = cost + squared2norm(z(:,k)-xV_s(:,k));
end
opti.minimize(cost);

constr = {};
% Initial condition constraint
for k =1:numVX
    constr{end+1} = F*z(:,k) <= yM;
end

% Robust Tube set constraint
for j = 1:v   % loop over vertices
    for i = 1:sys.nm   % loop over models
        constr{end+1} = F*(sys.A_convh{i}*V{j}*yM + sys.B_convh{i}*u(:,j)) + d <= yM;
    end
    constr{end+1} = sys.X.A * V{j}*yM <= sys.X.b;
    constr{end+1} = sys.U.A *u(:,j) <= sys.U.b;
end
constr{end+1} = E*yM <= 0;
% origin inclusion
constr{end+1} = yM >= 0;
opti.subject_to(constr);

% csd_opts = {'daqp',struct('error_on_fail',false),struct('eps_prox',1e-4)};
csd_opts = {'gurobi'};%,struct('outputflag',0)};

opti.solver(csd_opts{:});


hausFun = opti.to_function('HausdorffDist', ...
    [{F},{E},V(:)',{d},{yM}], {cost, yM});%,...
%{'F','E','V','d','yM_0'},{'sigma','yM'});
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