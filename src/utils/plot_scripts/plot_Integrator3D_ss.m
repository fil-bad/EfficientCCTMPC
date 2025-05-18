fig = figure("Renderer","painters","Units","centimeters");
fig.Position(3:4) = 4*[3,2.5]; % [width,height]
hold on

X.plot("wire",true,"linestyle",'--',"edgecolor",0.5*ones(1,3),"edgealpha",0.5)
han_CC = RoA_CCTMPC.plot("alpha",0.05,"EdgeColor",0.75*[0.4660 0.6740 0.1880],...
    "Color",[0.4660 0.6740 0.1880],"Linewidth",0.5,"edgealpha",0.5);

han_H = RoA_HTMPC.plot("Alpha",0.05,"EdgeColor",0.75*[0.9290 0.6940 0.1250],...
    "Color",[0.9290 0.6940 0.1250],"Linewidth",0.5,"edgealpha",0.5);

han_mRCI = Polyhedron(ccPoly.F,cctmpc.rciSol{1}).plot("Alpha",0.1,"Color",[0.8500 0.3250 0.0980],...
    "EdgeColor",0.75*[0.8500 0.3250 0.0980],"Linewidth",0.5,"edgealpha",0.5);

for i=1:length(x_sys_all)
    plot3(x_sys_all{i}(1,:),x_sys_all{i}(2,:),x_sys_all{i}(3,:),"LineWidth",0.3,"Linestyle","-","Color",[0 0.4470 0.7410])
    scatter3(x_sys_all{i}(1,1),x_sys_all{i}(2,1),x_sys_all{i}(3,1), 8, 0.2*ones(1,3),"filled")
end

% % bring grid in front of everything
Ax1 = gca;
Ax1.YGrid = 'on'; Ax1.Layer = 'top'; Ax1.GridAlpha = 0.05;
set(Ax1, 'TickLabelInterpreter', 'latex');

han_leg1 = legend(Ax1,[han_CC, han_H, han_mRCI],...
    {'$\mathcal{O}_{c}(3)$','$\mathcal{O}_h(3)$',"$P(y_m)$"}, ...
    'Interpreter','latex','Location','northeast');
han_leg1.FontSize = 12;

% xlabel('$x_1$','Interpreter','latex');
% ylabel('$x_2$','Interpreter','latex');
% zlabel('$x_3$','Interpreter','latex');

view(-200,35)
%%
% % minimize white borders around plot
% set(Ax1,'LooseInset', max(get(Ax1,'TightInset'), 0.01)) % remove border from axis
set(fig,'PaperPositionMode','Auto','PaperUnits',...
    'centimeters','PaperSize',fig.Position(3:4)) % resize pdf page

fig.Renderer = 'painters'; % a way to force saving in vector graphics
% % save the plot as PDF file
saveas(fig, '../figures/Integrator3D_Regions.pdf','pdf')
