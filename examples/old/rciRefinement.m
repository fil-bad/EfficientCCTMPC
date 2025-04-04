clc
clear all
%%
 load('model.mat')
 nx = 10;
 nu = 3;
 nw = 3;
g = 9.81;
dt = 0.1;
kT = 0.91;
Bw = [dt*eye(3); zeros(7,3)];
 HX = [eye(nx); -eye(nx)]; hX = repmat([4,10,10,5,5,5,pi/4,pi/4,pi/4,pi/4]',[2,1]);
 mX = size(HX,1);
HU = [eye(nu); -eye(nu)]; hU = [pi/6; pi/6; 2*g-g/kT;  pi/6; pi/6; -(0-g/kT)];
HW = [eye(nw); -eye(nw)]; hW = repmat([0.01; 0.01; 0.01],[2,1]);
epsw = 0.01*ones(nw,1);
X = Polyhedron(HX,hX);
X_vert = X.V';
nX_vert = size(X_vert,2);
Y = [-eye(nx); ones(1,nx)];
m=size(Y,1);
y_MPI= [zeros(nx,1);1]-Y*0.5*ones(nx,1)/nx;
%%
%Step 0
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
con_casadi = {};
con_X={};
con_U = {};
for i = 1:numTemplates
    i
    lambdas_loc = lambdas((i-1)*m+1:i*m,:);
    opti.subject_to(lambdas_loc(:)>=0)
    opti.subject_to(lambdas_loc*HW==Y*Winv*Bw)
    con_casadi{i} = Y * Winv * (A*W*template_vert(:,i)+B*u_vert(:,i))+lambdas_loc*hW<= y_MPI;
    opti.subject_to( con_casadi{i});
    con_X{i} = HX * W * template_vert(:,i) <= hX;
    opti.subject_to( con_X{i});
    con_U{i} = HU * u_vert(:,i) <= hU;
    opti.subject_to(con_U{i});
    vector = W*template_vert(:,i);
end
W_init = eye(nx);
opti.set_initial(W, W_init);
opti.set_initial(Winv, inv(W_init));
epsilon = opti.variable(m);
lambda_2 = opti.variable(m,mX);
opti.subject_to(lambda_2(:)>=0);
opti.subject_to(lambda_2*hX<=y_MPI+epsilon);
opti.subject_to(lambda_2*HX==Y*Winv);
cost = epsilon'*epsilon;
opti.minimize(cost);
p_opts = struct('expand', true);
s_opts = struct('tol', 1e-8, 'max_iter', 2000, 'print_level', 5);
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
polyopts.alpha = 0.1;
num_refine = 100;
Ys = {};
ys = {};
y_RPI = {};
Ys{1} = Y_invar;
ys{1} = y_MPI;
y_RPI{1} = y_MPI;
norm_dual_v = {};
max_norm_dual = [];
cost_compare = [];
id_interest = [];
for idxx = 1:num_refine
    %Configuration constraints ------------------------
    Y = Ys{idxx};
    y0 = ys{idxx};
    poly_refined = minHRep(Polyhedron(Y,y0));
    Y = poly_refined.A;
    y0 = poly_refined.b;
    Y = Y./y0;
    y0 = 0*y0+1;
    m = size(Y,1);
    WH_P=Polyhedron(Y,y0);
    WH_vert=WH_P.V';
    
    WC=nchoosek(1:size(Y,1),nx);
    V={};
    index_mat={};
    one_Vi={};
    for i=1:size(WH_vert,2)
        for j=1:size(WC,1)
            V_test=Y(WC(j,:),:);
            h_test=y0(WC(j,:));
            if norm(V_test*WH_vert(:,i)-h_test,'inf')<=1e-9
                index_mat{i}=WC(j,:)';
                ones_mat=zeros(nx,size(Y,1));
                for k=1:nx
                    ones_mat(k,WC(j,k))=1;
                end
                one_Vi{i}=ones_mat;
                V{i}=inv(V_test)*one_Vi{i};
            end
        end
    end
    E=[];
    for i=1:size(WH_vert,2)
        E=[E; Y*V{i}-eye(m)];
    end
    E=sparse(E);
    m_bar=length(V); %Number of vertices
    l=size(E,1);
    U_mats = {};
    for i=1:m_bar
        U_mats{i}=zeros(nu,nu*m_bar);
        U_mats{i}(:,(i-1)*nu+1:i*nu)=eye(nu);
    end
     d = [];
    for i=1:m
        ww = linprog(-(Y(i,:)*Bw)',HW,hW);
        d(i,1) = Y(i,:)*Bw*ww;
    end
    %Maximal parameterized RPI set  ------------------------
    import casadi.*
    opti = Opti();
    y_tilde = opti.variable(m);
    u_vert = opti.variable(nu,m_bar);
    con_casadi = {E*y_tilde<=0};
    con_X={};
    con_U = {};
    for i = 1:m_bar
        con_casadi{i} = Y*(A*V{i}*y_tilde+B*u_vert(:,i))+d<= y_tilde;
        opti.subject_to( con_casadi{i});
        con_X{i} = HX*V{i} *y_tilde<=hX;
        opti.subject_to( con_X{i});
        con_U{i} = HU*u_vert(:,i)<=hU;
        opti.subject_to(con_U{i});
    end
    b_pts = opti.variable(nx,nX_vert);
    cost = 0.;
    for i=1:nX_vert
        opti.subject_to(Y*b_pts(:,i)<=y_tilde);
        vector = b_pts(:,i)-X_vert(:,i);
        cost = cost+vector'*vector;
    end
    opti.minimize(cost);
    p_opts = struct('expand', true);
    s_opts = struct('tol', 1e-9, 'max_iter', 2000, 'print_level', 0);
    opti.solver('ipopt', p_opts, s_opts);
    sol = opti.solve();
    y_tilde = sol.value(y_tilde);
    cost_compare(:,idxx) = sol.value(cost)
    y_RPI{idxx} = y_tilde;
    v_RPI = [];
    for i=1:m_bar
        v_RPI(:,i) = V{i}*y_tilde;
    end
    %Get i_star ------------------------
    norm_duals = [];
    duals = [];
    for i=1:size(v_RPI,2)
        duals(:,i) = sol.value(opti.dual(con_casadi{i}));
        norm_duals(i) = norm(duals(:,i),2);
    end
    norm_dual_v{idxx} = norm_duals;
    max_norm_dual(idxx) = max(norm_duals);
    id_interest(idxx) = find(norm_duals==max(norm_duals));
    %Build new template ------------------------
    Y_dual = v_RPI';
    rhs = [];
    for i=1:size(Y_dual,1)
        xx =  linprog(-Y_dual(i,:)',Ys{idxx},ys{idxx});
        rhs(i,1) = Y_dual(i,:)*xx;
    end
    Ys{idxx+1} = [Y; Y_dual(id_interest(idxx),:)];
    ys{idxx+1} = [y_tilde; 0.99*rhs(id_interest(idxx))];
end