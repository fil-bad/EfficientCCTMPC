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
clear all
load('model.mat')

%%
nx = 10;
nu = 3;
nw = 3;
g = 9.81;
dt = 0.1;
kT = 0.91;
Bw = [dt*eye(3); zeros(7,3)];
HX = [eye(nx); -eye(nx)]; hX = repmat([4,4,2,5,5,3,pi/4,pi/4,pi,pi]',[2,1]);
HU = [eye(nu); -eye(nu)]; hU = [pi/6; pi/6; 2*g-g/kT;  pi/6; pi/6; -(0-g/kT)];
HW = [eye(nw); -eye(nw)]; hW = repmat([0.01; 0.09; 0.2],[2,1]);
epsw = 0.01*ones(nw,1);
X = Polyhedron(HX,hX);
X_vert = X.V';
nX_vert = size(X_vert,2)
%%
% take the polar set, having probability~=1 to have simple polytope (as we compute a simplicial one)
% F = Polyhedron(polyGenEulerAngles(nx, v),ones(v,1)).minHRep();
randSeed = 1;
for f=11%nx+1:nx+10
    v_s = numericalSpherePoints(nx, f, randSeed);
    F = Polyhedron(v_s,ones(f,1)).minHRep();
    disp("f: "+ f + ", v: " + size(uniquetol(F.V,1e-16,'ByRows',true),1))
    % F_v = Polyhedron(uniquetol(F.V,1e-16,'ByRows',true));
end

%%
% % For the simplex -----
Y = [-eye(nx); ones(1,nx)];
m=size(Y,1);
y_MPI= [zeros(nx,1);1]-Y*0.5*ones(nx,1)/nx;
F = Polyhedron(Y,y_MPI).minHRep();
% % ---------------------

Y = F.A./F.b;
m=size(Y,1);
y_MPI = ones(m,1);
WH_P=Polyhedron(Y,y_MPI);
template_vert=WH_P.V';

import casadi.*
opti = Opti();
W    = opti.variable(nx, nx);
Winv = opti.variable(nx, nx);

opti.subject_to( W * Winv == eye(nx) );
numTemplates = size(template_vert,2);
u_vert = opti.variable(nu,numTemplates);

lambda=opti.variable(m,2*nw);
w_hat = opti.variable(m,1);

opti.subject_to(lambda(:)>=0)
opti.subject_to(w_hat>=0)

opti.subject_to(lambda*hW <= w_hat)
opti.subject_to(lambda*HW==Y*Winv*Bw)

cost = 0;
for i = 1:numTemplates
    opti.subject_to( Y * Winv * (A*W*template_vert(:,i)+B*u_vert(:,i))+ w_hat <= y_MPI );
    opti.subject_to( HX * W * template_vert(:,i) <= hX );
    opti.subject_to( HU * u_vert(:,i) <= hU );
end


justFeasible = true;
if justFeasible
    cost = 0; %feasible solution
else
    epsilon = opti.variable(m,1);
    opti.subject_to(epsilon >= 0);

    cost = sum(epsilon); % 1-norm

    % using strong duality to reduce the number of constraints
    Lam_X = opti.variable(m,size(HX,1));
    opti.subject_to( Lam_X(:) >= 0);

    opti.subject_to( Lam_X*hX <= ones(m,1) + epsilon);
    opti.subject_to( Lam_X*HX == Y * Winv);
end

W_init = eye(nx);
opti.set_initial(W, W_init);
opti.set_initial(Winv, inv(W_init));

opti.minimize(cost);
p_opts = struct('expand', true);
s_opts = struct('tol', 1e-6, 'max_iter', 2000, 'print_level', 5);
opti.solver('ipopt', p_opts, s_opts);
sol = opti.solve();
W = sol.value(W);
Winv = sol.value(Winv);
u_vert=sol.value(u_vert);
Y_invar = Y*Winv;
S1 = minHRep(Polyhedron(Y_invar,y_MPI));
v_RPI = [];
for i=1:numTemplates
    v_RPI(:,i)=W*template_vert(:,i);
end
S1V = minVRep(S1);

%%
Y = S1.A./S1.b;
m = size(Y,1);
y = ones(m,1);
w_vert = Polyhedron(HW,hW).V';
nW = size(w_vert,2);
error=[];
x_nexts = [];
for i=1:size(v_RPI,2)
    for j=1:nW
        x_next = A*v_RPI(:,i)+B*u_vert(:,i)+Bw*w_vert(:,j);
        x_nexts = [x_nexts x_next];
        error = [error Y*x_next-y];
    end
end
max(max(error))


polyopts.alpha = 0.1;
figure
subplot(2,2,1)
idx = [1 2 3];
plot(polytope(projection(S1,idx)),polyopts)
hold on
scatter3(x_nexts(idx(1),:),x_nexts(idx(2),:),x_nexts(idx(3),:))
X.projection(idx).plot("alpha",0)

subplot(2,2,2)
idx = [4 5 6];
plot(polytope(projection(S1,idx)),polyopts)
hold on
scatter3(x_nexts(idx(1),:),x_nexts(idx(2),:),x_nexts(idx(3),:))
X.projection(idx).plot("alpha",0)

subplot(2,2,3)
idx = [7 8];
plot(polytope(projection(S1,idx)),polyopts)
hold on
scatter(x_nexts(idx(1),:),x_nexts(idx(2),:))
X.projection(idx).plot("alpha",0)

subplot(2,2,4)
idx = [9 10];
plot(polytope(projection(S1,idx)),polyopts)
hold on
scatter(x_nexts(idx(1),:),x_nexts(idx(2),:))
X.projection(idx).plot("alpha",0)

%%
% % For the simplex -----
% % Y = [-eye(nx); ones(1,nx)];
% % m=size(Y,1);
% % y_MPI= [zeros(nx,1);1]-Y*0.5*ones(nx,1)/nx;
% % ---------------------
% Y = F.A./F.b;
% m=size(Y,1);
% y_MPI = ones(m,1);
% WH_P=Polyhedron(Y,y_MPI);
% template_vert=WH_P.V';
% vert_rand = unique(randi(nX_vert,200,1));
% import casadi.*
% opti = Opti();
% W    = opti.variable(nx, nx);
% Winv = opti.variable(nx, nx);
% opti.subject_to( W * Winv == eye(nx) );
% numTemplates = size(template_vert,2);
% u_vert = opti.variable(nu,numTemplates);
% 
% lambda=opti.variable(m,2*nw);
% w_hat = opti.variable(m,1);
% 
% opti.subject_to(lambda(:)>=0)
% opti.subject_to(w_hat>=0)
% 
% opti.subject_to(lambda*hW <= w_hat)
% opti.subject_to(lambda*HW==Y*Winv*Bw)
% 
% cost = 0;
% for i = 1:numTemplates
%     opti.subject_to( Y * Winv * (A*W*template_vert(:,i)+B*u_vert(:,i))+ w_hat <= y_MPI );
%     opti.subject_to( HX * W * template_vert(:,i) <= hX );
%     opti.subject_to( HU * u_vert(:,i) <= hU );
%     vector = W*template_vert(:,i);
%     % cost = cost+vector'*vector;
% end
% 
% W_init = eye(nx);
% opti.set_initial(W, W_init);
% opti.set_initial(Winv, inv(W_init));
% b1_pts = opti.variable(nx, length(vert_rand));
% % for i = 1:length(vert_rand)
% %     i
% %     id_inter = vert_rand(i);
% %     opti.subject_to( Y * Winv * b1_pts(:,i) <= y_MPI );
% %     error= X_vert(:,id_inter)-b1_pts(:,i);
% %     cost = cost+error'*error;
% % end
% opti.minimize(cost);
% p_opts = struct('expand', true);
% s_opts = struct('tol', 1e-6, 'max_iter', 2000, 'print_level', 5);
% opti.solver('ipopt', p_opts, s_opts);
% sol = opti.solve();
% W = sol.value(W);
% Winv = sol.value(Winv);
% u_vert=sol.value(u_vert);
% Y_invar = Y*Winv;
% S1 = minHRep(Polyhedron(Y_invar,y_MPI));
% v_RPI = [];
% for i=1:numTemplates
%     v_RPI(:,i)=W*template_vert(:,i);
% end
% S1V = minVRep(S1);

