
close all
%%
% % create figure
fig = figure("Renderer","painters","Units","centimeters");
fig.Position(3:4) = 3*[3,3]; % [width,height]
hold on; grid on




vx_s = linspace(vx_interval(1),vx_interval(2),100);

n = 2;
P_curr_s = cell(1,n);
P_curr_s_han = cell(1,n);
% min P0 = 40/3.6 % max 60/3.6;
% min_vx = [52.5,57.5]/3.6;
min_vx = [14.5, 16.5]; % [m/s]


for i=1:n
P_curr_s{i} = Polyhedron([P_hat.A;-1,0],[P_hat.b;-min(min_vx(i),max(P_hat.V(:,1))-5/3.6 )]);
P_curr_s_han{i} = P_curr_s{i}.plot('alpha',0.1*i,"color",[0 0.4470 0.7410],'LineWidth',1.25,'EdgeColor',[0.4,0.4,0.4],'xdimension',2);
end
han_P0 = P_hat.plot('alpha',0.0,"color",[0 0.4470 0.7410],'LineWidth',1.25,'EdgeColor',[0.1,0.1,0.1],'xdimension',2);

% xl = xline(min_vx(2),'--');



han_trueModel = plot(vx_s,1./vx_s,LineWidth=1.25, Color="#D95319");

xlabel("$v^x\ [m/s]$",'Interpreter','latex'); ylabel("$(v^x)^{-1}\ [m/s]$",'Interpreter','latex')



% % bring grid in front of everything
Ax1 = gca;
Ax1.YGrid = 'on'; Ax1.Layer = 'top'; Ax1.GridAlpha = 0.05;
set(Ax1, 'TickLabelInterpreter', 'latex');



% % add legend inside the plot
lgd_han = legend([han_trueModel,han_P0,P_curr_s_han{1},P_curr_s_han{2}],...
    'True model',...
    '$\hat{\mathcal{P}}$',...
    '$\mathcal{P}(14.5)$',...
    '$\mathcal{P}(16.5)$',...
    'Interpreter','latex','Location','best','FontSize',10,'LineWidth',0.02);

% % minimize white borders around plot
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.01)) % remove border from axis
set(fig,'PaperPositionMode','Auto','PaperUnits',...
    'centimeters','PaperSize',fig.Position(3:4)) % resize pdf page

fig.Renderer = 'painters'; % a way to force saving in vector graphics
% % save the plot as PDF file
saveas(fig, '../figures/LateralDynamicsParamSpace.pdf','pdf')


