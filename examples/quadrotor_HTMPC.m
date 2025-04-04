clc
clear all

 load('model.mat')
 nx = 10;
 nu = 3;
 nw = 3;
g = 9.81;
dt = 0.1;
kT = 0.91;

Bw = [dt*eye(3); zeros(7,3)];

 HX = [eye(nx); -eye(nx)]; hX = repmat([4,4,2,10,10,5,pi/3,pi/3,pi,pi]',[2,1]);
 mX = size(HX,1);
HU = [eye(nu); -eye(nu)]; hU = [pi/4; pi/4; 2*g-g/kT;  pi/6; pi/6; -(0-g/kT)];
HW = [eye(nw); -eye(nw)]; hW = repmat([0.05; 0.05; 0.1],[2,1]);
epsw = 0.01*ones(nw,1);
X = Polyhedron(HX,hX);
X_vert = X.V';
nX_vert = size(X_vert,2);

%%
Y = [eye(nx);-eye(nx)];
y_MPI = ones(2*nx,1);
m = size(Y,1);

WH_P=Polyhedron(Y,y_MPI);
template_vert=WH_P.V';

import casadi.*

opti = Opti();
W    = opti.variable(nx, nx);
Winv = opti.variable(nx, nx);
opti.subject_to( W * Winv == eye(nx) );
numTemplates = size(template_vert,2);
u_vert = opti.variable(nu,numTemplates);
lambdas=opti.variable(m*numTemplates,2*nw);
con = {};
con_X={};
con_U = {};
for i = 1:numTemplates
    i
    lambdas_loc = lambdas((i-1)*m+1:i*m,:);
    opti.subject_to(lambdas_loc(:)>=0)
    opti.subject_to(lambdas_loc*HW==Y*Winv*Bw)
    con{i} = Y * Winv * (A*W*template_vert(:,i)+B*u_vert(:,i))+lambdas_loc*hW<= y_MPI;
    opti.subject_to( con{i});
    con_X{i} = HX * W * template_vert(:,i) <= hX;
    opti.subject_to( con_X{i});
    con_U{i} = HU * u_vert(:,i) <= hU;
    opti.subject_to(con_U{i});
    vector = W*template_vert(:,i);
end
W_init = eye(nx);
opti.set_initial(W, W_init);
opti.set_initial(Winv, inv(W_init));
% cost = sum(sum(W'*W));
cost = 0;

opti.minimize(cost);
p_opts = struct('expand', true);
s_opts = struct('tol', 1e-6, 'max_iter', 2000, 'print_level', 5);
opti.solver('ipopt', p_opts, s_opts);
sol = opti.solve();
W = sol.value(W);
Winv = sol.value(Winv);
u_vert=sol.value(u_vert);
Y_invar = Y*Winv;
save('quad_box.mat')

v_RPI = [];
for i=1:numTemplates
    v_RPI(:,i)=W*template_vert(:,i);
end

d = [];
for i=1:m
    ww = linprog(-(Y_invar(i,:)*Bw)',HW,hW);
    d(i,1) = Y_invar(i,:)*Bw*ww;
end

%%
Y = Y_invar./y_MPI;
m = size(Y,1);
y = ones(m,1);
index = {};
V = {};
x_recon = [];
residual = [];
m_bar = numTemplates;
for i=1:m_bar
    residual(:,i) = Y*v_RPI(:,i)-y;
    [~, sortedIdx] = sort(abs(residual(:, i)), 'ascend');
    index{i} = sortedIdx(1:nx);
    Y_rows = Y(index{i},:);
    V{i} = zeros(nx,m);
    V{i}(:,index{i}) = inv(Y_rows);
    x_recon(:,i) = V{i}*y;
end
E = [];
for i=1:m_bar
    E = [E; Y*V{i}-eye(m)];
end
cone = minHRep(Polyhedron(E,zeros(size(E,1),1)));
E = cone.H(:,1:end-1);
U = {};
for i=1:m_bar
    UU = zeros(nu,nu*m_bar);
    UU(:,(i-1)*nu+1:i*nu)=eye(nu);
    U{i} = UU;
end

d = [];
for i=1:m
    ww = linprog(-(Y(i,:)*Bw)',HW,hW);
    d(i,1) = Y(i,:)*Bw*ww;
end

%%
ym = sdpvar(m,1);
um = sdpvar(nu*m_bar,1);
con=[E*ym<=0];
for i=1:m_bar
    con=[con; Y*(A*V{i}*ym+B*U{i}*um)+d<=ym];
    con=[con; HX*V{i}*ym<=hX];
    con=[con; HU*U{i}*um<=hU];
end
cost = ym'*ym+um'*um;
optimize(con,cost,sdpsettings('solver','gurobi'));
ym = value(ym);
um = value(um);
mRCI = projection(Polyhedron(Y,ym),1:3);


%%
h_xm = [];
h_um = [];
for i=1:m_bar
    h_xm(:,i) = HX*V{i}*ym;
    h_um(:,i) = HU*U{i}*um;
end
h_xm = max(h_xm')';
h_um = max(h_um')';
H = [blkdiag(HX,HU) [h_xm;h_um]];
h = [hX; hU];

N = 20;
beta = 0.95;
z = sdpvar(nx,N+1,'full');
v = sdpvar(nu,N+1,'full');
alpha = sdpvar(1,N+1,'full');
x0_des = sdpvar(nx,1);
x0_proj = sdpvar(nx,1);
error = x0_des-x0_proj;
cost = error'*error;
con=[Y*x0_proj<=alpha(1)*ym+Y*z(:,1)];
con = [con; alpha>=0];
for t=1:N
    con=[con; Y*(A*z(:,t)+B*v(:,t))+(1-alpha(t))*d+alpha(t)*ym<=alpha(t+1)*ym+Y*z(:,t+1)];
    con=[con; H*[z(:,t);v(:,t);alpha(:,t)]<=h];
end
t=N+1;
con=[con; Y*(A*z(:,t)+B*v(:,t))+(1-alpha(t))*d+alpha(t)*ym<=(beta*alpha(t)+1-beta)*ym+beta*Y*z(:,t)];
con=[con; H*[z(:,t);v(:,t);alpha(:,t)]<=h];
x0_projector = optimizer(con,cost,sdpsettings('solver','quadprog','verbose',0),{x0_des},{x0_proj});


%%
vX = size(X_vert,2);
x0_projected = [];
distance = [];
parfor i=1:vX
    i
    x0_projected(:,i)=x0_projector(X_vert(:,i));
    distance(i) = norm(X_vert(:,i)-x0_projected(:,i),2);
end
hD = max(distance)
plot(x0_projected(1:3,:)')

%%
Q = eye(nx+nu+1);
P = Q/(1-beta^2);
z = sdpvar(nx,N+1,'full');
v = sdpvar(nu,N+1,'full');
alpha = sdpvar(1,N+1,'full');
x0 = sdpvar(nx,1);
cost = 0;
con=[Y*x0<=alpha(1)*ym+Y*z(:,1)];
con = [con; alpha>=0];
for t=1:N
    con=[con; Y*(A*z(:,t)+B*v(:,t))+(1-alpha(t))*d+alpha(t)*ym<=alpha(t+1)*ym+Y*z(:,t+1)];
    con=[con; H*[z(:,t);v(:,t);alpha(:,t)]<=h];
    vector = [z(:,t);v(:,t);alpha(:,t)-1];
    cost = cost+vector'*Q*vector;
end
t=N+1;
con=[con; Y*(A*z(:,t)+B*v(:,t))+(1-alpha(t))*d+alpha(t)*ym<=(beta*alpha(t)+1-beta)*ym+beta*Y*z(:,t)];
con=[con; H*[z(:,t);v(:,t);alpha(:,t)]<=h];
vector = [z(:,t);v(:,t);alpha(:,t)-1];
cost = cost+vector'*P*vector;
MPC = optimizer(con,cost,sdpsettings('solver','quadprog','verbose',0),{x0},{z,alpha,cost});

x_curr = sdpvar(nx,1);
u_curr = sdpvar(nu,1);
y_next = sdpvar(m,1);
con=[HU*u_curr<=hU];
con=[con; Y*(A*x_curr+B*u_curr)+d<=y_next];
cost = u_curr'*u_curr;
u_compute = optimizer(con,cost,sdpsettings('solver','quadprog','verbose',0),{x_curr,y_next},{u_curr});

%%
N_sim = 100;
W_vert = Polyhedron(HW,hW).V';
mW = size(W_vert,2);
index = 798;
x = x0_projected(:,index);
u = [];
w = [];
y = [];
zs = {};
alphas = {};
cost = [];
time = [];
for t=1:N_sim
    t
    tic
    soln = MPC{x(:,t)};
    time(t) = toc;
    zs{t} = soln{1};
    alphas{t} = soln{2};
    cost(t) = soln{3};
    y(:,t) = alphas{t}(1)*ym+Y*zs{t}(:,1);
    y_next = alphas{t}(2)*ym+Y*zs{t}(:,2);
    u(:,t) = u_compute{x(:,t),y_next};
    w(:,t) = W_vert(:,randi(mW));
    x(:,t+1)=A*x(:,t)+B*u(:,t)+Bw*w(:,t);
end


%%
sets_space = {};
parfor t=1:N_sim
    t
    if cost(t)>=3e-3
        sets_space{t} = projection(Polyhedron(Y,y(:,t)),1:3);
    end
end
polyopts1.alpha = 0.4;
polyopts1.color = 'red';
polyopts2.alpha = 0.2;
polyopts2.color = 'black';
polyopts3.alpha = 0.1;
polyopts3.color = 'white';
X_set = projection(Polyhedron(HX,hX),1:3);

%%
stabRegion = convhull(x0_projected(1:3,:)');

%%
fig = figure("Renderer","painters","Units","centimeters");
fig.Position(3:4) = 4*[3,2.5]; % [width,height]
hold on

X_set.plot("wire",true,"linestyle",'-',"edgecolor",0.5*ones(1,3),"edgealpha",0.5)
hold on


% plot(polytope(mRCI),polyopts1)
hold on
han_stab = trisurf(stabRegion,x0_projected(1,:),x0_projected(2,:),x0_projected(3,:),...
    "FaceColor",[0.4660 0.6740 0.1880],"FaceAlpha",0.05,"EdgeColor",[0.4660 0.6740 0.1880],...
    "LineWidth",0.05,"LineStyle","-");

% han_stab = Polyhedron(x0_projected(1:3,:)').plot();

% for t=1:N_mpc
%     if rem(t,5)==0 && cost(t)>=3e-3
%         plot(polytope(sets{t}),polyopts2)
%         hold on
%     end
%     pause(0.0001)
% end
idZeroSet = find(cost <= 3e-3);
for t=[1:floor((idZeroSet-1)/3),floor((idZeroSet-1)/3)+1:3:floor((idZeroSet-1)/2)+1,floor((idZeroSet-1)/2)+1:6:idZeroSet-1]
    % plot(polytope(sets{t}),polyopts2)
    han_H = sets_space{t}.plot("alpha",0.05,"EdgeColor",[0 0.4470 0.7410],"Color",[0 0.4470 0.7410],"Linewidth",0.5,"edgealpha",0.5);
    hold on

    pause(0.0001)
end

han_mRCI = mRCI.plot("Alpha",0.8,"Color",[0.8500 0.3250 0.0980],...
    "EdgeColor",0.0*[0.8500 0.3250 0.0980],"Linewidth",0.5,"edgealpha",0.5);

plot3(x(1,:),x(2,:),x(3,:),'black','LineWidth',1.5)
hold on
scatter3(x(1,2:end),x(2,2:end),x(3,2:end),15,'ko','filled')
hold on
scatter3(x(1,1),x(2,1),x(3,1),15,'ro','filled')
view(-170,5)

% text(x(1,1)+0.2,x(2,1)+0.1,x(3,1)-0.1,'$x_0$','Interpreter','latex')

Ax1 = gca;
Ax1.YGrid = 'on'; Ax1.Layer = 'top'; Ax1.GridAlpha = 0.05;
Ax1.TickLabelInterpreter = 'latex';

% han_leg1 = legend(Ax1,[han_CC, han_PartCC,han_H, han_SP,han_mRPI,han_mRCI],... han_ref],...
han_leg1 = legend(Ax1,[han_H, han_mRCI, han_stab],... han_ref],...
    {'$P(y_0^*)$','$P(y_m)$','$\mathcal{O}_h(20)$'}, ...
    'Interpreter','latex','Location','northeast');
han_leg1.FontSize = 12;

xlabel('$x_1$','Interpreter','latex');
ylabel('$x_2$','Interpreter','latex');
zlabel('$x_3$','Interpreter','latex');


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

saveas(fig, '../figures/QuadrotorBox_space.pdf','pdf')




% %% Velocities
% sets = {};
% parfor t=1:N_sim
%     t
%     if cost(t)>=3e-3
%         sets{t} = projection(Polyhedron(Y,y(:,t)),4:6);
%     end
% end
% polyopts1.alpha = 0.4;
% polyopts1.color = 'red';
% polyopts2.alpha = 0.2;
% polyopts2.color = 'black';
% polyopts3.alpha = 0.1;
% polyopts3.color = 'white';
% X_set = projection(Polyhedron(HX,hX),4:6);
% 
% %%
% stabRegion = convhull(x0_projected(4:6,:)');
% 
% %%
% fig = figure("Renderer","painters","Units","centimeters");
% fig.Position(3:4) = 4*[3,2.5]; % [width,height]
% hold on
% 
% X_set.plot("wire",true,"linestyle",'-',"edgecolor",0.5*ones(1,3),"edgealpha",0.5)
% hold on
% 
% 
% % plot(polytope(mRCI),polyopts1)
% hold on
% han_stab = trisurf(stabRegion,x0_projected(4,:),x0_projected(5,:),x0_projected(6,:),...
%     "FaceColor",[0.4660 0.6740 0.1880],"FaceAlpha",0.05,"EdgeColor",[0.4660 0.6740 0.1880],...
%     "LineWidth",0.1,"LineStyle",":");
% 
% % han_stab = Polyhedron(x0_projected(1:3,:)').plot();
% 
% % for t=1:N_mpc
% %     if rem(t,5)==0 && cost(t)>=3e-3
% %         plot(polytope(sets{t}),polyopts2)
% %         hold on
% %     end
% %     pause(0.0001)
% % end
% idZeroSet = find(cost <= 3e-3);
% for t=[2:floor((idZeroSet-1)/3),floor((idZeroSet-1)/3)+1:3:floor((idZeroSet-1)/2)+1,floor((idZeroSet-1)/2)+1:6:idZeroSet-1]
%     % plot(polytope(sets{t}),polyopts2)
%     disp(t)
%     han_H = sets{t}.plot("alpha",0.05,"EdgeColor",[0 0.4470 0.7410],"Color",[0 0.4470 0.7410],"Linewidth",0.5,"edgealpha",0.5);
%     hold on
% 
%     pause(0.0001)
% end
% 
% han_mRCI = projection(Polyhedron(Y,ym),4:6).plot("Alpha",0.8,"Color",[0.8500 0.3250 0.0980],...
%     "EdgeColor",0.0*[0.8500 0.3250 0.0980],"Linewidth",0.5,"edgealpha",0.5);
% 
% plot3(x(4,:),x(5,:),x(6,:),'black','LineWidth',1.5)
% hold on
% scatter3(x(4,:),x(5,:),x(6,:),20,'ko','filled')
% hold on
% 
% view(-170,5)
% 
% Ax1 = gca;
% Ax1.YGrid = 'on'; Ax1.Layer = 'top'; Ax1.GridAlpha = 0.05;
% Ax1.TickLabelInterpreter = 'latex';
% 
% % % han_leg1 = legend(Ax1,[han_CC, han_PartCC,han_H, han_SP,han_mRPI,han_mRCI],... han_ref],...
% % han_leg1 = legend(Ax1,[han_H, han_mRCI, han_stab],... han_ref],...
% %     {'$P(y_0^*)$','$P(y_m)$','$\mathcal{O}_h(20)$'}, ...
% %     'Interpreter','latex','Location','northeast');
% % han_leg1.FontSize = 12;
% 
% xlabel('$x_4$','Interpreter','latex');
% ylabel('$x_5$','Interpreter','latex');
% zlabel('$x_6$','Interpreter','latex');
% 
% 
% % % minimize white borders around plot
% Ax1.LooseInset = max(get(Ax1,'TightInset'), 0.01);
% % set(Ax1,'LooseInset', max(get(Ax1,'TightInset'), 0.01)) % remove border from axis
% fig.PaperPositionMode = "auto";
% fig.PaperUnits = "centimeters";
% fig.PaperSize = fig.Position(3:4);
% % set(fig,'PaperPositionMode','Auto','PaperUnits',...
% %     'centimeters','PaperSize',fig.Position(3:4)) % resize pdf page
% 
% fig.Renderer = 'painters'; % a way to force saving in vector graphics
% % % save the plot as PDF file
% 
% saveas(fig, '../figures/QuadrotorBox_vel.pdf','pdf')








