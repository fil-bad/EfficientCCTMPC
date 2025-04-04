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
%%
clc
clear all
%%
load('model.mat')
nx = 10;
nu = 3;
nw = 3;


% C = [eye(ny), zeros(ny,nx-ny)];
% C = eye(nx);
ny = nx;

g = 9.81;
dt = 0.1;
kT = 0.91;

A_convh = {A}; B_convh = {B};
Bw = [dt*eye(3); zeros(7,3)];
HX = [eye(nx); -eye(nx)]; hX = repmat([4,4,2,10,10,5,pi/2,pi/2,5,5]',[2,1]);
mX = size(HX,1);

HU = [eye(nu); -eye(nu)]; hU = [pi/6; pi/6; 2*g-g/kT;  pi/6; pi/6; -(0-g/kT)];
mU = size(HU,1);
HW = [eye(nw); -eye(nw)]; hW = repmat([0.05; 0.05; 0.1],[2,1]);

% HY = [eye(ny); -eye(ny)]; hY = repmat([4,4,2]',[2,1]); mY = size(HY,1);
% HY = [eye(ny); -eye(ny)]; hY = repmat([4,4,2]',[2,1]); mY = size(HY,1);
C = eye(nx); ny = nx;
HY = HX;
hY = hX;
mY = size(HY,1);



Y_out = Polyhedron(HY,hY);
% Y_out.minHRep().minVRep();

Y_vert = Y_out.V';
nY_vert = size(Y_vert,2);
Y_vert_vec = reshape(Y_vert, ny*nY_vert,1);


F_tilde = [-eye(nx); ones(1,nx)];%./norm(ones(1,nx))];
m=size(F_tilde,1);
f_MPI= [zeros(nx,1);1]-F_tilde*0.5*ones(nx,1)/nx;

% assert(all(f_MPI > 0))
% F_tilde = F_tilde./f_MPI;
f_MPI = ones(m,1);
% F_tilde = Polyhedron(F_tilde,f_MPI).minVRep.A;
% f_MPI = ones(m,1);

% F_tilde = Polyhedron(numericalSpherePoints(nx, nx+1,1)).dual.minHRep().A;
% m=size(F_tilde,1);
% f_MPI=ones(m,1);





%%
nx = 2;
nu = 1;
nw = 2;
A = [1 0.2; -0.2 0.8];
B = [0; 0.2];

A_convh = {A}; B_convh = {B};
Bw = eye(nx);
HX = [eye(nx); -eye(nx)];
hX = [7.7;10;10;5.7];
mX = size(HX,1);
HU = [eye(nu); -eye(nu)];
hU = [10;10];
mU = size(HU,1);
HW = [eye(nx);-eye(nx)];
hW = 0.1*ones(2*nw,1);
% X = Polyhedron(HX,hX);
% X_vert = X.V';
% nX_vert = size(X_vert,2);
% X_vert_vec = reshape(X_vert, nx*nY_vert,1);

C = eye(nx); ny = nx;
HY = HX;
hY = hX;
mY = size(HY,1);

Y_out = Polyhedron(HY,hY);
Y_out.minHRep().minVRep();

Y_vert = Y_out.V';
nY_vert = size(Y_vert,2);
Y_vert_vec = reshape(Y_vert, ny*nY_vert,1);

F_tilde = [-eye(nx); ones(1,nx)./norm(ones(1,nx))];
m=size(F_tilde,1);
f_MPI= [zeros(nx,1);1]-F_tilde*0.5*ones(nx,1)/nx;

assert(all(f_MPI > 0))
F_tilde = F_tilde./f_MPI;
f_MPI = ones(m,1);


%%


h = 1/4;

A = [1 h (h^2)/2; 0 1 h; 0 0 1];                nx = size(A,2);
B = [(h^3)/6; (h^2)/2; h];                      nu = size(B,2);
Bw = [h (h^2)/2 (h^3)/6; 0 h (h^2)/2; 0 0 h];   nw = size(Bw,2);

% Constraints
HX = [eye(nx); -eye(nx)]; hX = 4*ones(2*nx,1); mX = size(HX,1);
HU = [eye(nu); -eye(nu)]; hU = [2; 2]; mU = size(HU,1);
HW = [eye(nw); -eye(nw)]; hW = (1/20)*ones(2*nw,1);


% Define LTI system
A_convh = {0.95*A, 1.05*A}; B_convh = {0.95*B, 1.05*B};

X = Polyhedron(HX, hX);
U = Polyhedron(HU, hU);

W_dist = Polyhedron(HW, hW);

C = eye(nx); ny = nx;
HY = HX;
hY = hX;
mY = size(HY,1);

Y_out = Polyhedron(HY,hY);
Y_out.minHRep().minVRep();

Y_vert = Y_out.V';
nY_vert = size(Y_vert,2);
Y_vert_vec = reshape(Y_vert, ny*nY_vert,1);


% F_tilde = Polyhedron(numericalSpherePoints(nx, nx+1,1),ones(nx+1,1)).minVRep.A;
% m=size(F_tilde,1);
% f_MPI=ones(m,1);


F_tilde = [-eye(nx); ones(1,nx)];%./norm(ones(1,nx))];
m=size(F_tilde,1);
f_MPI= [zeros(nx,1);1]-F_tilde*0.5*ones(nx,1)/nx;

% assert(all(f_MPI > 0))
% F_tilde = F_tilde./f_MPI;
% F_tilde = Polyhedron(F_tilde).minVRep.A;
f_MPI = ones(m,1);


% Polyhedron(F_tilde,f_MPI)

 
%%
%Step 0
justFeasible = false;

WH_P=Polyhedron(F_tilde,f_MPI);
template_vert=WH_P.V';
import casadi.*
opti = Opti();
W    = opti.variable(nx, nx);
Winv = opti.variable(nx, nx);
opti.subject_to( W * Winv == eye(nx) );
numTemplates = size(template_vert,2);
u_vert = opti.variable(nu,numTemplates);
lambdas=opti.variable(m*numTemplates,2*nw);

for i = 1:numTemplates
    lambdas_loc = lambdas((i-1)*m+1:i*m,:);
    opti.subject_to(lambdas_loc(:)>=0)
    opti.subject_to(lambdas_loc*HW ==F_tilde*Winv*Bw)
    for nm = 1:length(A_convh)
        opti.subject_to( ...
            F_tilde * Winv * (A_convh{nm}*W*template_vert(:,i)+B_convh{nm}*u_vert(:,i))+lambdas_loc*hW<= f_MPI);
    end
    opti.subject_to( HY * C * W * template_vert(:,i) <= hY );
    opti.subject_to(HU * u_vert(:,i) <= hU);
    vector = W*template_vert(:,i);
end
W_init = eye(nx);
opti.set_initial(W, W_init);
opti.set_initial(Winv, inv(W_init));

if justFeasible
    cost = 0;
else

    if ny == nx
        if all(C == eye(nx))
        epsilon = opti.variable(m);
        lambda_2 = opti.variable(m,mX);
        opti.subject_to(lambda_2(:)>=0);
        opti.subject_to(lambda_2*hX<=f_MPI+epsilon);
        opti.subject_to(lambda_2*HX==F_tilde*Winv);
        end
    else
        epsilon = opti.variable(mY);
        lambda_y = opti.variable(mY,m);
        opti.subject_to(lambda_y(:)>=0);
        opti.subject_to(hY <= lambda_y*f_MPI + epsilon);
        opti.subject_to(lambda_y*F_tilde == HY*C*W);
    end

    opti.subject_to(epsilon >=0)
    cost = sum(epsilon);%'*epsilon;
end

opti.minimize(cost);
p_opts = struct('expand', true);
s_opts = struct('tol', 1e-9, 'max_iter', 2000, 'print_level', 5);
opti.solver('ipopt', p_opts, s_opts);
sol = opti.solve();
W = sol.value(W);
Winv = sol.value(Winv);
u_vert=sol.value(u_vert);
Y_invar = F_tilde*Winv;
S1 = Polyhedron(Y_invar,f_MPI).normalize.minHRep();
v_RPI = [];
for i=1:numTemplates
    v_RPI(:,i)=W*template_vert(:,i);
end
%%
v_RPI = [];
for i=1:numTemplates
    v_RPI(:,i)=W*template_vert(:,i);
end

% Y_invar = ccPoly.F;
% f_MPI = ccPoly.y_sigma;
% v_RPI = Polyhedron(Y_invar,f_MPI).minVRep.V';

polyopts.alpha = 0.1;
polyopts2.color = 'green';
polyopts2.alpha = 0.001;
regularization = 0;%1e-6;
figure
num_refine =10+1;
Ys = {};
ys = {};
y_RPI = {};
Ys{1,1} = Y_invar;
ys{1,1} = f_MPI;
y_RPI{1,1} = f_MPI;
norm_dual_v = {};
max_norm_dual = [];
cost_compare = {};
id_interest = [];

Ys_best = {}; ys_best = {};

vert_indices = {};

costOverTime = [];

for idxx = 1:num_refine
    disp("Iteration: " + idxx)
    num_candidates = length(Ys(idxx,:));
    V=cell(num_candidates,1);

    for c = 1:num_candidates
        disp(char(9) + "Candidate: " + c + "/"+num_candidates)
        Y = Ys{idxx,c};
        y0 = ys{idxx,c};

        if idxx == 1
            WH_vert = v_RPI;
        else
            Y_old = Y(1:end-1,:); Y_cut = Y(end,:);
            y_old = y0(1:end-1,:); y_cut = y0(end,:);

            [~,idx_active] = sort( abs(Y_old*v_RPI(:,vert_idx(c)) - y_old) );
            vCut_mask = false(m,1);
            vCut_mask(idx_active(1:nx)) = true(nx,1);

            combFacets = nchoosek(find(vCut_mask),nx-1);

            new_verts = [];
            for v = 1:nx
                F_vert = [Y_old(combFacets(v,:),:); Y_cut];
                new_verts = [new_verts, F_vert \ [y_old(combFacets(v,:),:); y_cut] ];
            end
            WH_vert  = [v_RPI(:,[1:vert_idx(c)-1, vert_idx(c)+1:end]), new_verts];
            WH_vert;
        end

        m_bar = size(WH_vert,2);
        m = size(Y,1);

        %Get index sets
        for i=1:m_bar

            [~,idx_active] = sort(abs(Y*WH_vert(:,i)-y0));
            Vi_mask = false(m,1);
            Vi_mask(idx_active(1:nx)) = true(nx,1);
            % compute 1_Vi to get Vi = inv(F)*1_Vi
            one_mat = zeros(nx,m); % nx x f
            one_mat(sub2ind(size(one_mat),1:nx,find(Vi_mask)')) = 1;

            V{c,i} = Y(Vi_mask,:) \ one_mat;
        end

        E=[];
        for i=1:size(WH_vert,2)
            E=[E; Y*V{c,i}-eye(m)];
        end
        E=sparse(E);
        l=size(E,1);
        E_min = minHRep(Polyhedron(E,sparse(l,1)));
        E = sparse(E_min.A);

        U_mats = {};
        for i=1:m_bar
            U_mats{i}=zeros(nu,nu*m_bar);
            U_mats{i}(:,(i-1)*nu+1:i*nu)=eye(nu);
        end

        W_dist = Polyhedron(HW,hW);
        d  = W_dist.support((Y*Bw)');

        %Maximal parameterized RPI set  ------------------------
        A_ineq = [];
        b_ineq = [];
        for i=1:m_bar
            A_ineq = [A_ineq; Y*A*V{c,i}-eye(m) Y*B*U_mats{i} sparse(m,nx*nY_vert)];
            b_ineq = [b_ineq; -d];
        end
        for i=1:m_bar
            A_ineq = [A_ineq; HY*C*V{c,i} sparse(mY,nu*m_bar+nx*nY_vert)];
            b_ineq = [b_ineq; hY];
        end
        for i=1:m_bar
            A_ineq = [A_ineq; sparse(mU,m) HU*U_mats{i} sparse(mU,nx*nY_vert)];
            b_ineq = [b_ineq; hU];
        end
        A_ineq = [A_ineq; E sparse(size(E,1),nu*m_bar+nx*nY_vert)];
        b_ineq = [b_ineq; zeros(size(E,1),1)];

        A_ineq = [A_ineq; repmat(-speye(m),nY_vert,1) sparse(m*nY_vert,nu*m_bar) kron(speye(nY_vert),Y)];
        b_ineq = [b_ineq; zeros(m*nY_vert,1)];

        % % % to enforce all the facets to be positive! (hence zero included by definition)
        % A_ineq = [A_ineq; -speye(m) sparse(m,nu*m_bar+nx*nY_vert)];
        % b_ineq = [b_ineq; -0.001*ones(m,1)];

        selector = [sparse(ny*nY_vert,m+nu*m_bar) kron(speye(nY_vert),C)];
        Q_curr = selector'*selector+regularization*speye(m+nu*m_bar+nx*nY_vert);
        c_curr = -2*selector'*Y_vert_vec;
        n_vars = m + nu*m_bar + nx*nY_vert;

        model = [];
        model.Q    = Q_curr;
        model.obj  = c_curr;
        model.A    = -A_ineq;
        model.rhs  = -b_ineq;
        model.sense= repmat('>', size(b_ineq,1), 1);
        model.lb   = -inf(n_vars,1);
        model.ub   = inf(n_vars,1);
        model.vtype= repmat('C', n_vars, 1);
        params.outputflag = 0;
        params.FeasibilityTol = 1e-9;
        params.OptimalityTol  = 1e-9;
        result = gurobi(model, params);
        if strcmp(result.status, 'OPTIMAL')
            soln_QP = result.x;
            y_RPI{idxx,c} = soln_QP(1:m);
            u_RPI = reshape(soln_QP(m+1:m+nu*m_bar),nu,m_bar);

            error = selector*soln_QP-Y_vert_vec;
            cost_compare{idxx,c} = error'*error;
        elseif strcmp(result.status, 'SUBOPTIMAL')
            str = "Gurobi did not find an optimal solution. Status: " +  result.status;
            warning(str);
            soln_QP = result.x;
            error = selector*soln_QP-Y_vert_vec;
            cost_compare{idxx,c} = error'*error;
        else
            str = "Gurobi did not find an optimal solution. Status: " +  result.status;
            warning(str);
            cost_compare{idxx,c} = inf;
        end
    end

    % select the best improvement (minimal cost)
    [min_cost, best_idx] = min([cost_compare{idxx,:}]);
    if isinf(min_cost)
        warning("It wasn't possible to improve the template with this choice of hyperplanes.")
    end
    Ys_best{idxx} = Ys{idxx,best_idx};
    ys_best{idxx} = y_RPI{idxx,best_idx};

    v_RPI = [];
    for i=1:m_bar
        v_RPI(:,i) = V{best_idx,i}*y_RPI{idxx,best_idx};
    end

    % perpendicular normals
    Y_normal = [];
    for v=1:m_bar
        Y_normal(v,:) =  mean(Ys{idxx,best_idx}(any(V{best_idx,v},1),:),1);
        Y_normal(v,:) = Y_normal(v,:);%./norm(Y_normal(v,:));
    end

    V;

    Ys_poly=Polyhedron(Ys{idxx,best_idx},y_RPI{idxx,best_idx});
    rhs = Ys_poly.support(Y_normal');

    % prepare next batch to be tested
    vert_idx = zeros(1,m_bar);
    for v=1:m_bar
        supp_vs = max(Y_normal(v,:)*v_RPI,[],1); % value at each vertex in that direction
        [supp_vs,supp_idx] = sort(supp_vs,'descend');
        vert_idx(v) = supp_idx(1);
        zero_supp = Y_normal(v,:)*zeros(nx,1);

        % cut_value = 0.01*supp_vs(1) + 0.99*max(zero_supp, supp_vs(find( (supp_vs(1)-supp_vs)/supp_vs(1) > 1e-6, 1)));
        cut_value = 0.1*supp_vs(1) + 0.9*max(zero_supp, supp_vs(2));

        Ys{idxx+1,v} = [Ys{idxx,best_idx}; Y_normal(v,:)];
        ys{idxx+1,v} = [y_RPI{idxx,best_idx}; cut_value ];
    end


    % if nx==2
    %     subplot(1,2,1)
    %     Polyhedron(HX,hX).plot("alpha",0)
    %     hold on
    %     plot(polytope(Y_normal,rhs),polyopts2)
    %     hold on
    %     plot(polytope(Ys_best{idxx},ys_best{idxx}),polyopts)
    %     hold on
    %     % TODO: maintain the color if the number of facet increased,
    %     % highlight the cost in red for which the number of facet remained
    %     % constant!
    %     scatter(v_RPI(1,best_idx),v_RPI(2,best_idx),'bo')
    %     hold on
    % else
    %     dims = 1:3;
    %     subplot(1,2,1)
    %     plot(polytope(Y_normal,rhs).projection(dims),polyopts2)
    %     hold on
    %     plot(polytope(Ys_best{idxx},ys_best{idxx}).projection(dims),polyopts)
    % end
    subplot(1,2,2)
    costOverTime = [costOverTime; cost_compare{idxx,best_idx}];
    scatter(idxx, cost_compare{idxx,best_idx},'bo');
    hold on
    pause(0.001)
end

%%

fig = figure("Renderer","painters","Units","centimeters");
fig.Position(3:4) = 3*[2,2.5]; % [width,height]
hold on


if nx == 2
    MRCI = ULTISystem('A', {A}, 'B', {B}, 'E',Bw).minVRep.invariantSet('X',X,'U',Polyhedron(HU,hU),'D',Polyhedron(HW,hW));
    MRCI.plot("alpha",0)
end

Polyhedron(HX,hX).plot("wire",true,"linestyle",'--',"edgecolor",0.5*ones(1,3),"edgealpha",0.5)

limit = min(20,length(Ys_best)-1);
han_best = Polyhedron(Ys_best{limit}, ys_best{limit}).minVRep.plot(...
    "alpha",0.1,"Color",[0.4660 0.6740 0.1880],"EdgeColor",0.75*[0.4660 0.6740 0.1880],"Linewidth",0.75);
han_start = Polyhedron(Ys_best{1}, ys_best{1}).minVRep.plot(...
    "alpha",0.2,"Color",[0.8500 0.3250 0.0980],"EdgeColor",0.75*[0.8500 0.3250 0.0980],"Linewidth",0.75);


Ax1 = gca;
Ax1.YGrid = 'on'; Ax1.Layer = 'top'; Ax1.GridAlpha = 0.05;
Ax1.TickLabelInterpreter = 'latex';

% Ax1.CameraPosition = Ax1.CameraPosition + [35, 145, 0];
view(-170,15)

% str_best = "$\mathcal{P}(F_{"+ limit + "},y_{"+limit+"}^*)$";
% han_leg = legend(Ax1,[han_start, han_best],...
%     {'$\mathcal{P}(F_0,y_0^*)$',str_best}, ...
%     'Interpreter','latex','Location','northeast');
han_leg = legend(Ax1,[han_start, han_best],...
    {'Initial Set',"Final Set"}, ...
    'Interpreter','latex','Location','northeast');
han_leg.FontSize = 10;

xlabel('$x_1$','Interpreter','latex');
ylabel('$x_2$','Interpreter','latex');
zlabel('$x_3$','Interpreter','latex');

%%
% % minimize white borders around plot
Ax1.LooseInset = max(get(Ax1,'TightInset'), 0.01);
% set(Ax1,'LooseInset', max(get(Ax1,'TightInset'), 0.01)) % remove border from axis
fig.PaperPositionMode = "auto";
fig.PaperUnits = "centimeters";
fig.PaperSize = fig.Position(3:4);
% set(fig,'PaperPositionMode','Auto','PaperUnits',...
%     'centimeters','PaperSize',fig.Position(3:4)) % resize pdf page

fig.Renderer = 'painters'; % a way to force saving in vector graphics
% % save the plot as PDF file

% % save the plot as PDF file
saveas(fig, '../figures/RCI_afterIters.pdf','pdf')




%%
fig = figure("Renderer","painters","Units","centimeters");
fig.Position(3:4) = 3*[2,2.5]; % [width,height]
hold on

Polyhedron(HX,hX).plot("wire",true,"linestyle",'--',"edgecolor",0.5*ones(1,3),"edgealpha",0.5)

polyChoice = 1;
Polyhedron(Ys_best{polyChoice}, ys_best{polyChoice}).plot("wire",true,"linestyle",':',"edgealpha",0.5);

[~ ,idx_best_cut] = min([cost_compare{polyChoice+1,:}]);
Polyhedron(Ys{polyChoice+1,idx_best_cut}, ys{polyChoice+1,idx_best_cut}).plot(...
    "alpha",0.2,"Color",[0.8500 0.3250 0.0980],"EdgeColor",0.5*[0.8500 0.3250 0.0980],"Linewidth",0.75);

h_facets = Polyhedron( Ys{polyChoice+1,idx_best_cut}, ys{polyChoice+1,idx_best_cut}).minHRep.getFacet;
han_facet = h_facets(end).minVRep().plot("alpha",0);
hatchfill2(han_facet,"HatchColor",0.75*[0.8500 0.3250 0.0980],"HatchAngle",0)

vert_cut = Polyhedron(Ys_best{polyChoice}, ys_best{polyChoice}).V(end,:)';
scatter3(vert_cut(1),vert_cut(2),vert_cut(3),25,"MarkerEdgeColor",[0.4660 0.6740 0.1880],"LineWidth",1.25)


Ax1 = gca;
Ax1.YGrid = 'on'; Ax1.Layer = 'top'; Ax1.GridAlpha = 0.05;
Ax1.TickLabelInterpreter = 'latex';

% Ax1.CameraPosition = Ax1.CameraPosition + [35, 145, 0];
view(-170,15)

xlabel('$x_1$','Interpreter','latex');
ylabel('$x_2$','Interpreter','latex');
zlabel('$x_3$','Interpreter','latex');

%%
% % minimize white borders around plot
Ax1.LooseInset = max(get(Ax1,'TightInset'), 0.01);
% set(Ax1,'LooseInset', max(get(Ax1,'TightInset'), 0.01)) % remove border from axis
fig.PaperPositionMode = "auto";
fig.PaperUnits = "centimeters";
fig.PaperSize = fig.Position(3:4);
% set(fig,'PaperPositionMode','Auto','PaperUnits',...
%     'centimeters','PaperSize',fig.Position(3:4)) % resize pdf page

fig.Renderer = 'painters'; % a way to force saving in vector graphics

% % save the plot as PDF file
saveas(fig, '../figures/RCI_Cut.pdf','pdf')

%%

fig = figure("Renderer","painters","Units","centimeters");
fig.Position(3:4) = 3*[2,2.5]; % [width,height]
hold on

Polyhedron(HX,hX).plot("wire",true,"linestyle",'--',"edgecolor",0.5*ones(1,3),"edgealpha",0.5)

polyChoice = 1;
Polyhedron(Ys_best{polyChoice}, ys_best{polyChoice}).plot(...
    "alpha",0.2,"Color",[0.8500 0.3250 0.0980],"EdgeColor",0.5*[0.8500 0.3250 0.0980],"Linewidth",0.75);

vert_cut = Polyhedron(Ys_best{polyChoice}, ys_best{polyChoice}).V(end,:)';
scatter3(vert_cut(1),vert_cut(2),vert_cut(3),25,"MarkerEdgeColor",[0.4660 0.6740 0.1880],"LineWidth",1.25)

Polyhedron(Ys_best{polyChoice+1}, ys_best{polyChoice+1}).plot("alpha",0.1,...
    "alpha",0.1,"Color",[0 0.4470 0.7410],"EdgeColor",0.75*[0 0.4470 0.7410],"Linewidth",0.75);


Ax1 = gca;
Ax1.YGrid = 'on'; Ax1.Layer = 'top'; Ax1.GridAlpha = 0.05;
Ax1.TickLabelInterpreter = 'latex';

% Ax1.CameraPosition = Ax1.CameraPosition + [35, 145, 0];
view(-170,15)

xlabel('$x_1$','Interpreter','latex');
ylabel('$x_2$','Interpreter','latex');
zlabel('$x_3$','Interpreter','latex');

%%
% % minimize white borders around plot
Ax1.LooseInset = max(get(Ax1,'TightInset'), 0.01);
% set(Ax1,'LooseInset', max(get(Ax1,'TightInset'), 0.01)) % remove border from axis
fig.PaperPositionMode = "auto";
fig.PaperUnits = "centimeters";
fig.PaperSize = fig.Position(3:4);
% set(fig,'PaperPositionMode','Auto','PaperUnits',...
%     'centimeters','PaperSize',fig.Position(3:4)) % resize pdf page

fig.Renderer = 'painters'; % a way to force saving in vector graphics

saveas(fig, '../figures/RCI_OneIter.pdf','pdf')


%%
fig = figure("Renderer","painters","Units","centimeters");
fig.Position(3:4) = 3*[2,2.5]; % [width,height]

scatter(0:length(costOverTime)-1, costOverTime,'filled');

xlabel('Iterations (i)','Interpreter','latex');
ylabel('$\sigma^i$','Interpreter','latex');

ylim([180,224])

Ax1 = gca;
Ax1.YGrid = 'on'; Ax1.Layer = 'top'; Ax1.GridAlpha = 0.05;
Ax1.TickLabelInterpreter = 'latex';

xlabel('$i$','Interpreter','latex');
ylabel('$\sigma^i$','Interpreter','latex');


% % minimize white borders around plot
Ax1.LooseInset = max(get(Ax1,'TightInset'), 0.01);
% set(Ax1,'LooseInset', max(get(Ax1,'TightInset'), 0.01)) % remove border from axis
fig.PaperPositionMode = "auto";
fig.PaperUnits = "centimeters";
fig.PaperSize = fig.Position(3:4);
% set(fig,'PaperPositionMode','Auto','PaperUnits',...
%     'centimeters','PaperSize',fig.Position(3:4)) % resize pdf page

fig.Renderer = 'painters'; % a way to force saving in vector graphics

saveas(fig, '../figures/Cost_iterations.pdf','pdf')



% ~,~


%% Quadrotor 

fig = figure("Renderer","painters","Units","centimeters");
fig.Position(3:4) = 4*[2,2.5]; % [width,height]
hold on

idx = 1:3;

polyChoice = 1;
Polyhedron(Ys_best{polyChoice}, ys_best{polyChoice}).projection(1:3).plot("alpha",0.2,"Color",[0.8500 0.3250 0.0980]);

Polyhedron(Ys_best{polyChoice+1}, ys_best{polyChoice+1}).projection(1:3).plot("alpha",0.05,...
    "Color",[0.3010 0.7450 0.9330]);%,"edgealpha",0.75);



