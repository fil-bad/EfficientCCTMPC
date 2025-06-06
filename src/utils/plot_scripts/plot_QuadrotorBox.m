%%
fig = figure("Renderer","painters","Units","centimeters");
fig.Position(3:4) = 4*[3,2.5]; % [width,height]
hold on

X.projection(1:3).plot("wire",true,"linestyle",'-',"edgecolor",0.5*ones(1,3),"edgealpha",0.5)
hold on

projFeas = Polyhedron(Polyhedron(new_htmpc.convh_x0s').projection(1:3).V);
han_feas = projFeas.plot("alpha",0.05,"EdgeColor",[0.4660 0.6740 0.1880],...
    "Color",[0.4660 0.6740 0.1880],"Linewidth",0.1,"edgealpha",0.5);


%%
idZeroSet = find(Lyap_cost <= 1e-6); idZeroSet = idZeroSet(1);
sets_space = {};
% for i=1:idZeroSet
%     disp(i)
%     sets_space{end+1} = Polyhedron(Polyhedron(box_ccPoly.F,y0_CL{i}).projection(1:3).V);
% end

% since it's homothetic, it's enough computing ,mRCI set do z+a*mRCI
mRCI = Polyhedron(Polyhedron(box_ccPoly.F,y_box).projection(1:3).V);
for i=1:idZeroSet
    disp(i)
    sets_space{end+1} = z0_CL{i}(1:3,1)+(alpha0_CL{i}(1)*mRCI);
    % sets_space{end+1} = Polyhedron(Polyhedron(box_ccPoly.F,y0_CL{i}).projection(1:3).V);
end


%%
idxs_plot = [1:round((idZeroSet-1)/4),...
    round((idZeroSet-1)/4)+1:2:round((idZeroSet-1)/3)+1,...
    round((idZeroSet-1)/3)+1:3:round((idZeroSet-1)/2)+1,...
    round((idZeroSet-1)/2)+1:4:idZeroSet-1];
idxs_plot = [1:10,12:2:20,23:3:32,36:4:40];
for t = idxs_plot
    % plot(polytope(sets{t}),polyopts2)
    han_H = sets_space{t}.plot("alpha",0.05,"EdgeColor",[0 0.4470 0.7410],"Color",[0 0.4470 0.7410],"Linewidth",0.5,"edgealpha",0.5);
    hold on

    % pause(0.0001)
end


han_mRCI = mRCI.plot("Alpha",0.8,"Color",[0.8500 0.3250 0.0980],...
    "EdgeColor",0.0*[0.8500 0.3250 0.0980],"Linewidth",0.5,"edgealpha",0.5);



%%
x_t_pos = x_sys(1:3,:);
plot3(x_t_pos(1,:),x_t_pos(2,:),x_t_pos(3,:),'black','LineWidth',1.2)
hold on
scatter3(x_t_pos(1,2:end),x_t_pos(2,2:end),x_t_pos(3,2:end),10,'ko','filled')
hold on
scatter3(x_t_pos(1,1),x_t_pos(2,1),x_t_pos(3,1),20,'ro','filled')
view(-4,8)

% text(x(1,1)+0.2,x(2,1)+0.1,x(3,1)-0.1,'$x_0$','Interpreter','latex')

Ax1 = gca;
Ax1.YGrid = 'on'; Ax1.Layer = 'top'; Ax1.GridAlpha = 0.05;
Ax1.TickLabelInterpreter = 'latex';

han_leg1 = legend(Ax1,[han_H, han_mRCI, han_feas],...
    {'$P(y_0^*)$','$P(y_m)$','$\mathcal{O}_h(20)$'}, ...
    'Interpreter','latex','Location','northeast');
han_leg1.FontSize = 12;

% xlabel('$x_1$','Interpreter','latex');
% ylabel('$x_2$','Interpreter','latex');
% zlabel('$x_3$','Interpreter','latex');


% % minimize white borders around plot
% Ax1.LooseInset = max(get(Ax1,'TightInset'), 0.01);
% set(Ax1,'LooseInset', max(get(Ax1,'TightInset'), 0.01)) % remove border from axis
fig.PaperPositionMode = "auto";
fig.PaperUnits = "centimeters";
fig.PaperSize = fig.Position(3:4);
% set(fig,'PaperPositionMode','Auto','PaperUnits',...
%     'centimeters','PaperSize',fig.Position(3:4)) % resize pdf page

fig.Renderer = 'painters'; % a way to force saving in vector graphics
% % save the plot as PDF file

saveas(fig, '../figures/QuadrotorBox_space.pdf','pdf')

