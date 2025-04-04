% NOTE: this script must be called from the main one. this separation has
% been done in order to clean up the main code (as the plotting part is
% less relevant compared to the core code. For this reason, a simple check
% is done to avoid running this script directly. Obviously we're not
% responsible for a different use (that it's up to you.)

db_check = dbstack;
% if length(db_check) < 3 && ...
%         ~any(strcmp(db_check(end).name,{'LTI_2D','evaluateCode'}))
%     error("Script was called directly. Execute the example one first.")
% end

% % create figure
fig = figure("Renderer","painters","Units","centimeters");
fig.Position(3:4) = 4*[3,2.5]; % [width,height]

hold on

% % plot state constraints
rectangle('Position',[min(X.V(:,1))-2, max(X.V(:,2)),15,2],'FaceColor',0.8*ones(3,1),'LineStyle','none')
rectangle('Position',[max(X.V(:,1)), min(X.V(:,2))-2,2,15],'FaceColor',0.8*ones(3,1),'LineStyle','none')
rectangle('Position',[min(X.V(:,1))-2, min(X.V(:,2))-2,15,2],'FaceColor',0.8*ones(3,1),'LineStyle','none');
rectangle('Position',[min(X.V(:,1))-2, min(X.V(:,2))-2,2,15],'FaceColor',0.8*ones(3,1),'LineStyle','none');
X.plot('alpha',0,'LineWidth',1.5,'EdgeColor',[0.5,0.5,0.5],'xdimension',2)


% plot optimal DRTO_set (assumed constant for plot simplicity)
for j = 1:T_orbit
    han_yper = Polyhedron(ccPoly.F, DRTO_ys{1}(:,j)).plot('alpha',0.5,'Linewidth',0.1,'color',[0.85 0.325 0.098]);%[0.466 0.674 0.188]);
end
pause;

% % plot closed-loop CCPolytopes
for t = 1:N_mpc
    % Polyhedron(ccPoly.F, ys_tot(:,j)).plot('alpha',1,'Linewidth',0.1,'color',[1,1,1]); % white background
    han_ys = Polyhedron(ccPoly.F, OCP_ys{t}(:,1)).plot('alpha',0.5,'Linewidth',0.1,'color',[0.466 0.674 0.188]);%[0.85 0.325 0.098]);
   
    % Polyhedron(ccPoly.F, ys_tot(:,j)).plot('alpha',1,'Linewidth',0.1,'color',[1,1,1]); % white background
    han_y = Polyhedron(ccPoly.F, OCP_y{t}(:,1)).plot('alpha',0.5,'Linewidth',0.1,'color',[0 0.447 0.741]);%[0.929 0.694 0.125]);
    
    pause(0.25)
end


% reference to be tracked
% if nx==ny
% han_ref = plot([r_sys(1,:)], [r_sys(2,:)],'-o',...
%     'Color','k','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor',[0.6350 0.0780 0.1840]);
% end
% % add labels and axes limit
xlabel("$x_1$",'Interpreter','latex'); ylabel("$x_2$",'Interpreter','latex')
% xlim([min(X.V(:,1))-0.2, max(X.V(:,1)+0.2)]); 
xlim([min(X.V(:,1))-0.2, max(X.V(:,1)+0.2)]); 

ylim([min(X.V(:,2))-0.2, max(X.V(:,2))+0.2])

% % bring grid in front of everything
Ax1 = gca;
Ax1.YGrid = 'on'; Ax1.Layer = 'top'; Ax1.GridAlpha = 0.05;
set(Ax1, 'TickLabelInterpreter', 'latex');

% % add first legend inside the plot for CCPolytopes
han_leg1 = legend(Ax1,[han_y, han_ys, han_yper],... han_ref],...
    {'$X(y^*_0(p_t))$','$X(z^*_0(p_t))$','$X(z^{\mathrm{o}}(p_t))$','$\ \ \ \ p_t$'}, ...
    'Interpreter','latex','Location','southwest');

% % minimize white borders around plot
set(Ax1,'LooseInset', max(get(Ax1,'TightInset'), 0.01)) % remove border from axis
set(fig,'PaperPositionMode','Auto','PaperUnits',...
    'centimeters','PaperSize',fig.Position(3:4)) % resize pdf page

fig.Renderer = 'painters'; % a way to force saving in vector graphics
% % save the plot as PDF file
saveas(fig, '../figures/Integrator2D_Periodic.pdf','pdf')