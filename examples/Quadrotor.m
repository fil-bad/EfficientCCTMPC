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
Ts = 0.1;
[A,B,Bw] = defDiscreteQuadrotorDyn(Ts);
nx = size(A,2); nu = size(B,2); nw = size(Bw,2);

global g kT
% Constraints
% states: [x y z vx vy vz pitch roll wpitch wroll]
HX = [eye(nx); -eye(nx)]; hX = repmat([4,4,2,10,10,5,pi/3,pi/3,pi,pi]',[2,1]);
% HX = [[eye(ny); -eye(ny)], zeros(2*ny,nx-ny) ]; hX = repmat([4,4,2]',[2,1]);

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

% CCpolytope
F_tilde = [-eye(nx); ones(1,nx)];%./norm(ones(1,nx))];
m=size(F_tilde,1);
% f_MPI= [zeros(nx,1);1]-F_tilde*0.5*ones(nx,1)/nx;
% F_tilde = F_tilde./f_MPI;
f_MPI = ones(m,1);




%%
ccPoly_base = CCPolytope_test(sys, F_tilde, true, true);
%%
ccPoly = ccPoly_base;

%%
% perpendicular normals
Y_normal = [];
for i=1:ccPoly_base.v
    Y_normal(i,:) = mean(ccPoly_base.F(any(ccPoly_base.Vi_s{i},1),:),1);
    Y_normal(i,:) = Y_normal(i,:);%./norm(Y_normal(i,:));
end

Ys_poly=Polyhedron(ccPoly_base.F, f_MPI);
rhs = Ys_poly.support(Y_normal');
P_allcut = Polyhedron([ccPoly_base.F; Y_normal],[f_MPI; 0.9*rhs]).minVRep;
P_allcut
%%
ccPoly = CCPolytope_test(sys, P_allcut.A./P_allcut.b, false, true);

% ccPoly = CCPolytope_test(sys, Ys_best{end}./ys_best{end}, false);


% ccPoly = CCPolytope_test(sys, normalizePolytope(Ys_best{1}, ys_best{1}), false);
% ccPoly = CCPolytope_test(sys, S1.A./S1.b, false);



%%
% Cost Function definition: we can either pwass structs of matrices, or
% function handles defined externally. We only require that the cost
% functions are -strictly- convex, and that the function signatures are:
% RCI : @(ys,us)
% Stage/Term_cost : @(y_k/N,u_k/N,ys,us)

Qv = blkdiag( diag([1*ones(1,3), 1*ones(1,7)]), 0.01*eye(nu));
Qc = blkdiag(0.1*eye(nx), 0.1*eye(nu));
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
%%
% Q_bar = blkdiag(eye(nx),eye(nu),1);
% P_bar = (1/(1-gamma^2))*Q_bar;
%
% F_t_cost = [ccPoly.F, zeros(ccPoly.f,nu), y_m;
%             zeros(ccPoly.v*nu, nx), kron(ones(ccPoly.v,1),eye(nu)), u_m(:)];
%
% Q_mat = F_t_cost*Q_bar*F_t_cost';
% P_mat = F_t_cost*P_bar*F_t_cost';
%
% costFunMan = CostFunctionManager(ccPoly, Q_mat, P_mat, RCI_cost);
%%
% Build MPC scheme
sdp_opts = sdpsettings('solver','gurobi','verbose',0);%,'savesolverinput',1);

N_ocp = 3; % Note: y_N = N_ocp+1
cctmpc = CCTMPC(sys, ccPoly, costFunMan, N_ocp, gamma, sdp_opts); disp("CCTMPC done");

baseHausCC = cctmpc.hausDist;
%%
% htmpc_sp = HTMPC_SinglePolicy(sys, ccPoly, costFunMan, N_ocp, gamma, sdp_opts); disp("HTMPC done");


%% Simulation CCTMPC

N_mpc = 20*N_ocp+1; % simulation steps

% (state space) system dynamics
x_sys = zeros(nx, N_mpc+1);
% x_sys(:,1) = (1-0e-4)*cctmpc.findFeasibleX0([5,-5,-5,zeros(1,7)]',diag([0.1,0.1,0.1,1000*ones(1,7)]));
disp(x_sys(1:3,1)')
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

    % propagate system dynamics
    x_sys(:,t+1) = sys.step(x_sys(:,t), u_sys(:,t));

    % save computed data
    y0_CL{t} = cctmpc.y0_sol;
    Lyap_cost(t) = cctmpc.Lyapunov_cost;
end
timeCC = toc;
timeCC = timeCC/N_mpc*1000;
disp(timeCC)

%% Simulation CCTMPC

hausH = zeros(1,10);
timeH_s = zeros(1,10);
horIncr = 0:4;
for n=1:length(horIncr)

    htmpc_sp = HTMPC_SinglePolicy(sys, ccPoly, costFunMan, N_ocp+horIncr(n), gamma, sdp_opts); disp("HTMPC done");
    hausH(n) = htmpc_sp.hausDist;

    % (state space) system dynamics
    x_sys = zeros(nx, N_mpc+1);
    % x_sys(:,1) = (1-0e-4)*htmpc_sp.findFeasibleX0([5,-5,-5,zeros(1,7)]',diag([0.1,0.1,0.1,1000*ones(1,7)]));
    % x_sys(1:3,1) = [-2,-2,0]';
    disp(x_sys(1:3,1)')

    u_sys = zeros(nu, N_mpc);

    % (parameter space) optimal control dynamics
    y0_CL = cell(1,N_mpc);

    Lyap_cost = zeros(1,N_mpc);

    tic
    % main loop
    for t = 1:N_mpc
        disp("MPC iteration: " + t +"/"+ N_mpc)
        % compute input for the system
        % u_sys(:,t) = htmpc_sp.solve(x_sys(:,t));
        try
            u_sys(:,t) = htmpc_sp.solve(x_sys(:,t));
        catch
            disp(t + " point out of feasible.")
            % timeH_s(n) = toc;
            % timeH_s(n) = timeH_s(n)/t*1000;
            break;
            % u_sys(:,t) = htmpc_sp.solve((1-2e-2)*x_sys(:,t));
        end

        % propagate system dynamics
        x_sys(:,t+1) = sys.step(x_sys(:,t), u_sys(:,t));

        % save computed data
        y0_CL{t} = htmpc_sp.y0_sol;
        Lyap_cost(t) = htmpc_sp.Lyapunov_cost;
    end
    timeH_s(n) = toc;
    timeH_s(n) = timeH_s(n)/t*1000; % to catch infeasible iterations
    disp(timeH_s(n))
end

%%

fig = figure("Renderer","painters","Units","centimeters");
fig.Position(3:4) = 4*[3,1.5]; % [width,height]
hold on

% --- Plot data for left y-axis ---
yyaxis left                   % Activate left y-axis
% Here we use a light blue dashed line to connect the points,
% while the markers remain a full blue color.
plot(N_ocp+horIncr, hausH(1:length(horIncr)), ...
    'LineStyle', ':', ...                         % Dashed connecting line
    'LineWidth', 0.2, ...                           % Very thin line
    'Marker', 'o', ...                              % Marker style: circle
    'MarkerFaceColor', [0 0.4470 0.7410], ...                     % Filled markers in blue
    'MarkerSize', 5,...
    'Color', 0.75*[0 0.4470 0.7410]);                          % Light blue connecting line
ylabel("$d(\mathcal{X},\mathcal{O}_h(N))$",'Interpreter','latex');          % Left y-axis label
% Draw horizontal reference line for left y-axis using the same light blue
yline(baseHausCC, '-', "$d(\mathcal{X},\mathcal{O}_c(3))$", 'Color', 0.75*[0 0.4470 0.7410],'Interpreter','latex');
set(gca, 'YColor',  0.75*[0 0.4470 0.7410],'TickLabelInterpreter', 'latex');                           % Left axis tick labels in blue

% --- Plot data for right y-axis ---
yyaxis right                  % Switch to right y-axis
% set(gca, 'YScale', 'log')
% Similar approach: light red dashed line with fully-colored red markers.
plot(N_ocp+horIncr, timeH_s(1:length(horIncr)), ...
    'LineStyle', ':', ...                         % Dashed connecting line
    'LineWidth', 0.2, ...                           % Very thin line
    'Marker', 'square', ...                              % Marker style: circle
    'MarkerFaceColor', [0.8500 0.3250 0.0980], ...                     % Filled markers in blue
    'MarkerSize', 5,...
    'Color', 0.75*[0.8500 0.3250 0.0980]);                          % Light blue connecting line
% Draw horizontal reference line for right y-axis using the same light red,
% and position its label to the left instead of at the end.
yline(timeCC, '-', "CCTMPC(3)", 'Color', 0.75*[0.8500 0.3250 0.0980], 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'bottom','Interpreter','latex');
ylabel('HTMPC($N$) [ms]','Interpreter','latex');           % Right y-axis label
set(gca, 'YColor', 0.75*[0.8500 0.3250 0.0980],'TickLabelInterpreter', 'latex');                           % Right axis tick labels in red
ylim([120,650])
% Common settings for the entire plot
xlabel('Prediction Horizon HTMPC($N$)','Interpreter','latex');
xlim([3-0.04,7+0.04])
% --- Set x-axis ticks to display each integer ---
xVals = N_ocp + horIncr;
xMin = floor(min(xVals));
xMax = ceil(max(xVals));
xticks(xMin:1:xMax);  % Create an integer tick for every number in the range

Ax1 = gca;
% % minimize white borders around plot
set(Ax1,'LooseInset', max(get(Ax1,'TightInset'), 0.01)) % remove border from axis
set(fig,'PaperPositionMode','Auto','PaperUnits',...
    'centimeters','PaperSize',fig.Position(3:4)) % resize pdf page

fig.Renderer = 'painters'; % a way to force saving in vector graphics
% % save the plot as PDF file
saveas(fig, '../figures/HausTimeComparison.pdf','pdf')




%% Compute now relevant trajectories

% (state space) system dynamics
x_sys = zeros(nx, N_mpc+1);
x_sys(:,1) = (1-0e-4)*htmpc_sp.findFeasibleX0([-5,5,-5,zeros(1,7)]',diag([0.1,0.1,0,1000*ones(1,7)]));
disp(x_sys(1:3,1)')

u_sys = zeros(nu, N_mpc);

% (parameter space) optimal control dynamics
y0_CL = cell(1,N_mpc);
Lyap_cost = zeros(1,N_mpc);

% main loop
for t = 1:N_mpc
    disp("MPC iteration: " + t +"/"+ N_mpc)
    % compute input for the system
    u_sys(:,t) = htmpc_sp.solve(x_sys(:,t));

    % propagate system dynamics
    x_sys(:,t+1) = sys.step(x_sys(:,t), u_sys(:,t));

    % save computed data
    y0_CL{t} = htmpc_sp.y0_sol;
    Lyap_cost(t) = htmpc_sp.Lyapunov_cost;
end


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
RoA_CCTMPC_space = intersect(cctmpc.computeFeasRegion(templ).projection(1:3),X_space).minVRep();

% RoA_HTMPC_space = htmpc.computeFeasRegion(templ).projection(1:3);
%%
idx = 1:3;
% templ = zeros(size(dir3D,1),nx); templ(:,idx) = dir3D;
% RoA_CCTMPC_idx = cctmpc.computeFeasRegion(templ).projection(idx); disp("CC_done.")

fig = figure("Renderer","painters","Units","centimeters");
fig.Position(3:4) = 4*[3,2.5]; % [width,height]
hold on
han_X = X_space.plot("wire",true,"linestyle",'--',"edgecolor",0.5*ones(1,3),"edgealpha",0.5);

% RoA_CCTMPC_space.plot("alpha",0.02,"EdgeColor",[0.4660 0.6740 0.1880],"Color",[0.4660 0.6740 0.1880],"Linewidth",1.25);

xlabel('$x(t)$','Interpreter','latex');
ylabel('$y(t)$','Interpreter','latex');
zlabel('$z(t)$','Interpreter','latex');

% project the trajectory dynamics
xt_idx = x_sys(idx,:);

% view(az,el);

maxInd = find(Lyap_cost <= 5e-4);
maxInd = N_mpc;

idxxx = [1 3 6 12 16 24 ];
for i=1:length(idxxx)%floor(maxInd/4):maxInd
    t = idxxx(i);
    disp(t)
    Poly_y0_proj = Polyhedron(ccPoly.F,y0_CL{t}).projection(idx);
    if size(Poly_y0_proj.V,1) == 1
        scatter3(xt_idx(1,t),xt_idx(2,t),xt_idx(3,t),100,[0 0.4470 0.7410],"Marker",".")
    elseif Poly_y0_proj.volume < 1e-10
        scatter3(xt_idx(1,t),xt_idx(2,t),xt_idx(3,t),100,[0 0.4470 0.7410],"Marker",".")
    else
        han_CC = Poly_y0_proj.plot("alpha",0.05,"EdgeColor",[0 0.4470 0.7410],"Color",[0 0.4470 0.7410],"Linewidth",0.5);
    end
end
han_CC = Polyhedron(ccPoly.F,y0_CL{maxInd(1)}).projection(idx).plot("alpha",0.05,"EdgeColor",[0 0.4470 0.7410],"Color",[0 0.4470 0.7410],"Linewidth",0.5);

han_mRCI = Polyhedron(ccPoly.F,cctmpc.rciSol{1}).projection(idx).minVRep().plot("Alpha",0.2,"Color",[0.9290 0.6940 0.1250],"Linewidth",0.5);

plot3(xt_idx(1,:),xt_idx(2,:),xt_idx(3,:),"Color",[0.8500 0.3250 0.0980],"LineWidth",1.25)


view(-50,15)


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
saveas(fig, '../figures/QuadrotorHTMPC_space.pdf','pdf')

%%
idx = 4:6;
fig = figure("Renderer","painters","Units","centimeters");
fig.Position(3:4) = 3*[2.5,2.5]; % [width,height]
hold on


xlabel('$v_x(t)$','Interpreter','latex');
ylabel('$v_y(t)$','Interpreter','latex');
zlabel('$v_z(t)$','Interpreter','latex');

X_vel.plot("wire",true,"linestyle",'--',"edgecolor",0.5*ones(1,3),"edgealpha",0.5);
% project the trajectory dynamics
xt_idx = x_sys(idx,:);
plot3(xt_idx(1,:),xt_idx(2,:),xt_idx(3,:),"Color",[0.8500 0.3250 0.0980],"LineWidth",1.25)

idxxx = [1 3 6 12 16 24 ];
for i=1:length(idxxx)%floor(maxInd/4):maxInd
    t = idxxx(i);
    disp(t)
    Poly_y0_proj = Polyhedron(ccPoly.F,y0_CL{t}).projection(idx);
    if size(Poly_y0_proj.V,1) == 1
        scatter3(xt_idx(1,t),xt_idx(2,t),xt_idx(3,t),100,[0 0.4470 0.7410],"Marker",".")
    elseif Poly_y0_proj.volume < 1e-10
        scatter3(xt_idx(1,t),xt_idx(2,t),xt_idx(3,t),100,[0 0.4470 0.7410],"Marker",".")
    else
        han_CC = Poly_y0_proj.plot("alpha",0.05,"EdgeColor",[0 0.4470 0.7410],"Color",[0 0.4470 0.7410],"Linewidth",0.5);
    end
end


Polyhedron(ccPoly.F,cctmpc.rciSol{1}).projection(idx).plot("Alpha",0.2,"Color",[0.9290 0.6940 0.1250],"Linewidth",0.5);


view(-50,15)
%%
Ax1 = gca;
Ax1.YGrid = 'on'; Ax1.Layer = 'top'; Ax1.GridAlpha = 0.05;
set(Ax1, 'TickLabelInterpreter', 'latex');


% % minimize white borders around plot
set(Ax1,'LooseInset', max(get(Ax1,'TightInset'), 0.01)) % remove border from axis
set(fig,'PaperPositionMode','Auto','PaperUnits',...
    'centimeters','PaperSize',fig.Position(3:4)) % resize pdf page

fig.Renderer = 'painters'; % a way to force saving in vector graphics
% % save the plot as PDF file
saveas(fig, '../figures/QuadrotorHTMPC_vel.pdf','pdf')


%%
idx = 7:8;

fig = figure("Renderer","painters","Units","centimeters");
fig.Position(3:4) = 3*[2.5,2.5]; % [width,height]
hold on


xlabel('$\phi(t)$','Interpreter','latex');
ylabel('$\theta(t)$','Interpreter','latex');

X_angles.plot("wire",true,"linestyle",'--',"edgecolor",0.5*ones(1,3),"edgealpha",0.5);
% project the trajectory dynamics
xt_idx = x_sys(idx,:);
plot(xt_idx(1,:),xt_idx(2,:),"Color",[0.8500 0.3250 0.0980])

idxxx = [1 3 6 12 16 24 ];
for i=1:length(idxxx)%floor(maxInd/4):maxInd
    t = idxxx(i);
    disp(t)
    Poly_y0_proj = Polyhedron(ccPoly.F,y0_CL{t}).projection(idx);
    if size(Poly_y0_proj.V,1) == 1
        scatter(xt_idx(1,t),xt_idx(2,t),100,[0 0.4470 0.7410],"Marker",".")
    elseif Poly_y0_proj.volume < 1e-10
        scatter(xt_idx(1,t),xt_idx(2,t),100,[0 0.4470 0.7410],"Marker",".")
    else
        han_CC = Poly_y0_proj.plot("alpha",0.05,"EdgeColor",[0 0.4470 0.7410],"Color",[0 0.4470 0.7410],"Linewidth",0.5);
    end
end

Polyhedron(ccPoly.F,cctmpc.rciSol{1}).projection(idx).plot("Alpha",0.2,"Color",[0.9290 0.6940 0.1250],"Linewidth",0.5);
%%
Ax1 = gca;
Ax1.YGrid = 'on'; Ax1.Layer = 'top'; Ax1.GridAlpha = 0.05;
set(Ax1, 'TickLabelInterpreter', 'latex');


% % minimize white borders around plot
set(Ax1,'LooseInset', max(get(Ax1,'TightInset'), 0.01)) % remove border from axis
set(fig,'PaperPositionMode','Auto','PaperUnits',...
    'centimeters','PaperSize',fig.Position(3:4)) % resize pdf page

fig.Renderer = 'painters'; % a way to force saving in vector graphics
% % save the plot as PDF file
saveas(fig, '../figures/QuadrotorHTMPC_angles.pdf','pdf')


%%
idx = 9:10;
fig = figure("Renderer","painters","Units","centimeters");
fig.Position(3:4) = 3*[2.5,2.5]; % [width,height]
hold on
% X_angRate.plot("Alpha",0,"EdgeColor",0.4*ones(1,3),"Linewidth",1.5)

xlabel('$\dot{\phi}(t)$','Interpreter','latex');
ylabel('$\dot{\theta}(t)$','Interpreter','latex');

X_angRate.plot("wire",true,"linestyle",'--',"edgecolor",0.5*ones(1,3),"edgealpha",0.5);
% project the trajectory dynamics
xt_idx = x_sys(idx,:);
plot(xt_idx(1,:),xt_idx(2,:),"Color",[0.8500 0.3250 0.0980])

idxxx = [1 3 6 12 16 24 ];
for i=1:length(idxxx)%floor(maxInd/4):maxInd
    t = idxxx(i);
    disp(t)
    Poly_y0_proj = Polyhedron(ccPoly.F,y0_CL{t}).projection(idx);
    if size(Poly_y0_proj.V,1) == 1
        scatter(xt_idx(1,t),xt_idx(2,t),100,[0 0.4470 0.7410],"Marker",".")
    elseif Poly_y0_proj.volume < 1e-10
        scatter(xt_idx(1,t),xt_idx(2,t),100,[0 0.4470 0.7410],"Marker",".")
    else
        han_CC = Poly_y0_proj.plot("alpha",0.05,"EdgeColor",[0 0.4470 0.7410],"Color",[0 0.4470 0.7410],"Linewidth",0.5);
    end
end

Polyhedron(ccPoly.F,cctmpc.rciSol{1}).projection(idx).plot("Alpha",0.2,"Color",[0.9290 0.6940 0.1250],"Linewidth",0.5);
%%
Ax1 = gca;
Ax1.YGrid = 'on'; Ax1.Layer = 'top'; Ax1.GridAlpha = 0.05;
set(Ax1, 'TickLabelInterpreter', 'latex');


% % minimize white borders around plot
set(Ax1,'LooseInset', max(get(Ax1,'TightInset'), 0.01)) % remove border from axis
set(fig,'PaperPositionMode','Auto','PaperUnits',...
    'centimeters','PaperSize',fig.Position(3:4)) % resize pdf page

fig.Renderer = 'painters'; % a way to force saving in vector graphics
% % save the plot as PDF file
saveas(fig, '../figures/QuadrotorHTMPC_angleRates.pdf','pdf')




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




function [A_new, flipIndicator] = normalizePolytope(A, b)

tol = 1e-8;
[m, n] = size(A);
A_new = zeros(m, n);
flipIndicator = false(m,1);

for i = 1:m
    if abs(b(i)) < tol
        error('b(%d) is zero (or nearly zero), cannot normalize.', i);
    elseif b(i) > 0
        A_new(i,:) = A(i,:) / b(i);  % Normal case.
    else  % b(i) < 0
        warning('Inequality %d has b<0. Normalization will flip the inequality sign.', i);
        % Flipping: the original inequality
        %   A(i,:)*x <= b(i)
        % is replaced by
        %   (-A(i,:))*x <= -b(i)
        % and then normalized:
        A_new(i,:) = (-A(i,:)) / (-b(i));  % same as: -A(i,:)/abs(b(i))
        flipIndicator(i) = true;
    end
end
end









