% NOTE: this script must be called from the main one. this separation has
% been done in order to clean up the main code (as the plotting part is
% less relevant compared to the core code. For this reason, a simple check
% is done to avoid running this script directly. Obviously we're not
% responsible for a different use (that it's up to you.)
% db_check = dbstack;
% if length(db_check) < 3 && ...
%         ~any(strcmp(db_check(end).name,{'Ball_Stick','evaluateCode'}))
%     error("Script was called directly. Execute the example one first.")
% end
% close all
% % create figure
fig = figure("Renderer","painters","Units","centimeters");
fig.Position(3:4) = 3*[4,2]; % [width,height]
hold on; grid on

timesteps = 0:N_mpc-1;

% draw reference
plot(timesteps, r_sys, "Color",[0.85 0.325 0.098],"LineWidth",0.5)

% draw X projection for position
xb_interval = X.projection(1).V;
yline(min(xb_interval),"LineWidth",2);
yline(max(xb_interval),"LineWidth",2);

ylim(1.2*[min(xb_interval),max(xb_interval)])

% % draw closed-loop xb projection tubes
% y_MPC_0 (xb_y)
p_yMPC = fill([timesteps,timesteps(end:-1:1)], [min(projVerts.xb_y),max(projVerts.xb_y(:,end:-1:1))],...
    [0 0.447 0.741],'FaceAlpha',0.075,'EdgeColor',[0 0.447 0.741],"LineStyle",':');
% % ys_MPC (e_y)
p_ysMPC = fill([timesteps,timesteps(end:-1:1)], [min(projVerts.xb_ys),max(projVerts.xb_ys(:,end:-1:1))],...
    [0.466 0.674 0.188],"FaceAlpha",0.075,"EdgeColor",[0.466 0.674 0.188],"LineStyle",':');
% % y_rci (e_y)
p_rci = fill([timesteps,timesteps(end:-1:1)], [min(projVerts.xb_rci),max(projVerts.xb_rci(:,end:-1:1))],...
    [0.85 0.325 0.098],"FaceAlpha",0.125,"EdgeColor",[0.85 0.325 0.098],"LineStyle",':');

% draw x_b
plot(timesteps, x_sys(1,1:end-1),"Color",[0 0.447 0.741], "LineWidth",1)

xlabel("timesteps",'Interpreter','latex'); ylabel("$x_b\ [m]$",'Interpreter','latex')

% % bring grid in front of everything
Ax1 = gca;
Ax1.YGrid = 'on'; Ax1.Layer = 'top'; Ax1.GridAlpha = 0.05;
set(Ax1, 'TickLabelInterpreter', 'latex');

lgd_han = legend([p_yMPC,p_ysMPC,p_rci],...
    '$\mathrm{proj}_{x_b}(X(y^*_0(p_t)))$',...
    '$\mathrm{proj}_{x_b}(X(z^*_0(p_t)))$',...
    '$\mathrm{proj}_{x_b}(X(z^{\mathrm{o}}(p_t)))$',...
    'Interpreter','latex','Location','southeast','FontSize',10,'LineWidth',0.02);

% % minimize white borders around plot
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.01)) % remove border from axis
set(fig,'PaperPositionMode','Auto','PaperUnits',...
    'centimeters','PaperSize',fig.Position(3:4)) % resize pdf page

fig.Renderer = 'painters'; % a way to force saving in vector graphics
% % save the plot as PDF file
saveas(fig, '../figures/BallStick_TubeEvol_xB.pdf','pdf')


%% SHOW \Theta_P evolution
fig = figure("Renderer","painters","Units","centimeters");
fig.Position(3:4) = 3*[4,2]; % [width,height]
hold on; grid on

timesteps = 0:N_mpc-1;

% we don't have an output reference for this quantity
% plot(timesteps, r_sys, "Color",[0.85 0.325 0.098],"LineWidth",0.5)

% draw X projection for position
thP_interval = X.projection(3).V;
yline(min(thP_interval),"LineWidth",2);
yline(max(thP_interval),"LineWidth",2);

ylim(1.2*[min(thP_interval),max(thP_interval)])

% % draw closed-loop theta_p projection tubes
% y_MPC_0 (xb_y)
p_yMPC = fill([timesteps,timesteps(end:-1:1)], [min(projVerts.thP_y),max(projVerts.thP_y(:,end:-1:1))],...
    [0 0.447 0.741],'FaceAlpha',0.075,'EdgeColor',[0 0.447 0.741],"LineStyle",':');
% % ys_MPC (e_y)
p_ysMPC = fill([timesteps,timesteps(end:-1:1)], [min(projVerts.thP_ys),max(projVerts.thP_ys(:,end:-1:1))],...
    [0.466 0.674 0.188],"FaceAlpha",0.075,"EdgeColor",[0.466 0.674 0.188],"LineStyle",':');
% % y_rci (e_y)
p_rci = fill([timesteps,timesteps(end:-1:1)], [min(projVerts.thP_rci),max(projVerts.thP_rci(:,end:-1:1))],...
    [0.85 0.325 0.098],"FaceAlpha",0.125,"EdgeColor",[0.85 0.325 0.098],"LineStyle",':');

% draw theta_p
plot(timesteps, x_sys(3,1:end-1),"Color",[0 0.447 0.741], "LineWidth",1)

xlabel("timesteps",'Interpreter','latex'); ylabel("$\theta_p\ [rad]$",'Interpreter','latex')

% % bring grid in front of everything
Ax1 = gca;
Ax1.YGrid = 'on'; Ax1.Layer = 'top'; Ax1.GridAlpha = 0.05;
set(Ax1, 'TickLabelInterpreter', 'latex');

lgd_han = legend([p_yMPC,p_ysMPC,p_rci],...
    '$\mathrm{proj}_{\theta_p}(X(y^*_0(p_t)))$',...
    '$\mathrm{proj}_{\theta_p}(X(z^*_0(p_t)))$',...
    '$\mathrm{proj}_{\theta_p}(X(z^{\mathrm{o}}(p_t)))$',...
    'Interpreter','latex','Location','southeast','FontSize',10,'LineWidth',0.02);

% % minimize white borders around plot
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.01)) % remove border from axis
set(fig,'PaperPositionMode','Auto','PaperUnits',...
    'centimeters','PaperSize',fig.Position(3:4)) % resize pdf page

fig.Renderer = 'painters'; % a way to force saving in vector graphics
% % save the plot as PDF file
saveas(fig, '../figures/BallStick_TubeEvol_thP.pdf','pdf')

return

%%

% draw a two-lane road
center = round(max(X.V(:,1))+min(X.V(:,1)),5)/2;
roadwidth = 3;
yline(-roadwidth+center,"LineWidth",2);
hold on
yline(roadwidth+center,"LineWidth",2);
yline(center,"LineWidth",1.25,"LineStyle","--");

% draw reference
stairs(longDyn.pos, r_sys, "Color",[0.85 0.325 0.098],"LineWidth",1)
% draw e_y
plot(longDyn.pos, x_sys(1,1:end-1),"Color",[0 0.447 0.741], "LineWidth",1)

% % draw closed-loop e_y projection tubes
% y_MPC_0 (e_y)
p_yMPC = fill([longDyn.pos,longDyn.pos(end:-1:1)], [min(projVerts.ey_y),max(projVerts.ey_y(:,end:-1:1))],...
    [0 0.447 0.741],'FaceAlpha',0.075,'EdgeColor',[0 0.447 0.741],"LineStyle",':');
% ys_MPC (e_y)
p_ysMPC = fill([longDyn.pos,longDyn.pos(end:-1:1)], [min(projVerts.ey_ys),max(projVerts.ey_ys(:,end:-1:1))],...
    [0.466 0.674 0.188],"FaceAlpha",0.075,"EdgeColor",[0.466 0.674 0.188],"LineStyle",':');
% y_rci (e_y)
p_rci = fill([longDyn.pos,longDyn.pos(end:-1:1)], [min(projVerts.ey_rci),max(projVerts.ey_rci(:,end:-1:1))],...
    [0.85 0.325 0.098],"FaceAlpha",0.125,"EdgeColor",[0.85 0.325 0.098],"LineStyle",':');

% % draw admissible and current facing e_psi
epsi_init = 1; epsi_inter = 15;
% get position
x_quiv = longDyn.pos(epsi_init:epsi_inter:end-1);
y_quiv = x_sys(1,epsi_init:epsi_inter:end-2);

% % e_psi ccPoly projection
% up
u_quiv = cos(max(projVerts.epsi_y(:,epsi_init:epsi_inter:end-1),[],1));
v_quiv = sin(max(projVerts.epsi_y(:,epsi_init:epsi_inter:end-1),[],1));
h_up = quiver(x_quiv,y_quiv,u_quiv,v_quiv, 'Color',[0.466 0.674 0.188]);
h_up.Head.LineStyle = 'solid';
set(h_up,'AutoScale','on', 'AutoScaleFactor',0.25, 'LineWidth',0.75)
% down
u_quiv = cos(min(projVerts.epsi_y(:,epsi_init:epsi_inter:end-1),[],1));
v_quiv = sin(min(projVerts.epsi_y(:,epsi_init:epsi_inter:end-1),[],1));
h_down = quiver(x_quiv,y_quiv,u_quiv,v_quiv, 'Color',[0.466 0.674 0.188]);
h_down.Head.LineStyle = 'solid';
set(h_down,'AutoScale','on', 'AutoScaleFactor',0.25, 'LineWidth',0.75)
% current e_psi
th_quiv = x_sys(3,epsi_init:epsi_inter:end-2);
u_quiv = cos(th_quiv); v_quiv = sin(th_quiv);
h_curr = quiver(x_quiv,y_quiv,u_quiv,v_quiv, 'Color',[0.929 0.694 0.125]);
h_curr.Head.LineStyle = 'solid';
set(h_curr,'AutoScale','on', 'AutoScaleFactor',0.25, 'LineWidth',0.75)

% % add labels and axes limit
xlabel("$x\ [m]$",'Interpreter','latex'); ylabel("$y\ [m]$",'Interpreter','latex')
xlim([0, max(longDyn.pos)]); ylim([-roadwidth+center, roadwidth+center])

% % bring grid in front of everything
Ax1 = gca;
Ax1.YGrid = 'on'; Ax1.Layer = 'top'; Ax1.GridAlpha = 0.05;
set(Ax1, 'TickLabelInterpreter', 'latex');

% % add legend inside the plot
lgd_han = legend([p_yMPC,p_ysMPC,p_rci,h_up,h_curr],...
    '$\mathrm{proj}_{e_y}(X(y_0^*(x_t,r_t))$',...
    '$\mathrm{proj}_{e_y}(X(y_{\mathrm{s}}^*(x_t,r_t))$',...
    '$\mathrm{proj}_{e_y}(X(y_{\mathrm{o}}(r_t)))$',...
    '$\mathrm{proj}_{e_\psi}(X(y_0^*(x_t,r_t))$','$e_\psi(t)$',...
    'Interpreter','latex','Location','southeast','FontSize',10,'LineWidth',0.02);
lgd_han.Position(1:2) = [0.6,0.1925];

% % minimize white borders around plot
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.01)) % remove border from axis
set(fig,'PaperPositionMode','Auto','PaperUnits',...
    'centimeters','PaperSize',fig.Position(3:4)) % resize pdf page

fig.Renderer = 'painters'; % a way to force saving in vector graphics
% % save the plot as PDF file
saveas(fig, '../figures/LateralDynamics.pdf','pdf')