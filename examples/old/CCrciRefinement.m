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

clc
%% 3D system definition
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

f = nx+2; % dodecahedron

F = Polyhedron(numericalSpherePoints(nx, f,1)).normalize.dual.minHRep();
F.normalize()

%% 2D system definition
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

% By using polar coordinates for parameterization, we ensure to get a
% CCPolytope with minimal representation, i.e. f = v
f = nx+1;
F = Polyhedron(polyGenEulerAngles(nx, f),ones(f,1)).minHRep();
% F = Polyhedron(numericalSpherePoints(nx, f),ones(f,1)).normalize.dual.minHRep();
F.normalize()

%%
X = Polyhedron(HX, hX);
U = Polyhedron(HU, hU);
W_dist = Polyhedron(HW, hW);

sys = qLPV(A_convh, B_convh, Bw, X, U, W_dist);%, Acurr_han, Bcurr_han);
ccPoly = CCPolytope_test(sys, F.A, true,false);


plotVerticesSimplePolytope(ccPoly); % works in 2D and 3D

hold on
fig_han = gcf;

% ccPoly_best.ccpolytope.plot("alpha",0.1)
% 
% MRCI = ULTISystem('A', A_convh, 'B', B_convh, 'E',Bw).invariantSet('X',X,'U',U,'D',Polyhedron(HW,hW));
% MRCI.plot("alpha",0)

%%
% F = ccPoly.F;
% E = ccPoly.E;
% Vi_s = ccPoly.Vi_s{:};

ccPoly_best = ccPoly;

n_refinements = 10;
cut = 0.1;
for n = 1:n_refinements
    v = ccPoly_best.v;
    ccPoly_candidate = cell(1,v);

    vol_candidates = zeros(1,v);
    y_max = cell(1,v);
    cut = 0.5*cut;
    for i=1:v
        % F_tan =  mean(ccPoly_best.F(any(ccPoly_best.Vi_s{i},1),:),1);
        % F_tan = F_tan./norm(F_tan);
        % epsi = (1-cut)*ccPoly_best.ccpolytope.support((F_tan)');
        % F_cut = [ccPoly_best.F; F_tan./epsi]; 
        % Polyhedron(F_cut, ones(size(F_cut,1),1)).plot("alpha",0.1)
        
        F_dual = ccPoly_best.ccpolytope.dual.A(v,:);
        F_dual = F_dual./norm(F_dual);
        epsi = (1-cut)*ccPoly_best.ccpolytope.support((F_dual)');
        F_cut = [ccPoly_best.F; F_dual./epsi];

        % fb = [fb; 1];

        ccPoly_candidate{i} = CCPolytope_test(sys, F_cut, false);

        % for each of them, maximize the volume (hence min distance to X)
        [vol_candidates(i),y_max{i}] = getMaxVolume(ccPoly_candidate{i});
        % disp(vol_candidates(i))


    end
    disp(vol_candidates)
    [max_vol, idx_max] = min(vol_candidates);
    ccPoly_best = ccPoly_candidate{idx_max};
    ccPoly_best = CCPolytope_test(sys, ccPoly_best.F./y_max{idx_max},false);
    ccPoly_best.ccpolytope.plot("alpha",0.01);
    % Polyhedron(ccPoly_best.F,y_max{idx_max}).plot("alpha",0.01)
    drawnow;
end














%% --------------------------------------------------

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


% function E = computeEMatrix(obj)
%             % conic constraint defining the vertex configuration domain explicitly
%             E = zeros(obj.v*obj.f, obj.f);
%             for i=1:obj.v
%                 E((i-1)*obj.f+1:i*obj.f,:) = obj.F*obj.Vi_s{i}-eye(obj.f);
%             end
%         end
%
% function V = getVi_s(obj)
%             % computing Vi_s such that Vi*y0 is an isolated vertex of the CCPolytope
%             vertices = num2cell(obj.ccpolytope.V',1);
%             V = cell(size(vertices));
%
%             for i=1:length(vertices)
%                 % check for each individual vertex which constraint is active
%                 Vi_mask = abs(obj.F*vertices{i}-obj.ccpolytope.b) <= 1e-10;
%
%                 % assert(nnz(Vi_mask)==obj.sys.nx, "The vertex cannot be defined uniquely")
%
%                 % compute 1_Vi to get Vi = inv(F)*1_Vi
%                 one_mat = zeros(obj.sys.nx,obj.f); % nx x f
%                 one_mat(sub2ind(size(one_mat),1:obj.sys.nx,find(Vi_mask)')) = 1;
%
%                 V{i} = obj.F(Vi_mask,:) \ one_mat;
%             end
%         end


function [vol,ys] = getMaxVolume(ccPoly)

y_m = sdpvar(ccPoly.f,1);
u_m = sdpvar(ccPoly.sys.nu, ccPoly.v,'full');

% cost function

epsilon = sdpvar(ccPoly.f,1);
% 
constr = [];
constr = [constr; epsilon(:) >=0];

% X_verts = ccPoly.sys.X.V';
% 
% for v_x = 1:size(X_verts,2)
%     constr = [constr; ccPoly.F*X_verts(:,v_x) <= ones(ccPoly.f,1) + epsilon];
% end

HX = ccPoly.sys.X.A; hX = ccPoly.sys.X.b;
Lam_X = sdpvar(ccPoly.f,size(HX,1),'full');
constr = [constr; Lam_X(:) >= 0;
    Lam_X*hX <= ones(ccPoly.f,1) + epsilon;
    Lam_X*HX == ccPoly.F];

cost = sum(epsilon);

% constr = [constr; y_m >=0];
% cost = -sum(y_m);


% constraints


constr = [constr; constrSetS(ccPoly,y_m,u_m,y_m)];

% solve the OCP
optimize(constr,cost,sdpsettings('solver','gurobi','verbose',0));

ys = value(y_m);
us = value(u_m);
vol = value(cost);




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















