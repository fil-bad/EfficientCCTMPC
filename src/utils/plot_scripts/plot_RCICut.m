fig = figure("Renderer","painters","Units","centimeters");
fig.Position(3:4) = 3*[2,2.5]; % [width,height]
hold on

X.plot("wire",true,"linestyle",'--',"edgecolor",0.5*ones(1,3),"edgealpha",0.5)

polyChoice = 1;
P_iter0 = Polyhedron(save_Triple{polyChoice}{1}, save_yM{polyChoice});
P_iter0.plot("wire",true,"linestyle",':',"edgealpha",0.5);

P_iter1_feas = Polyhedron(save_Triple{polyChoice+1}{1}, [save_yM{polyChoice};y_cut_s{polyChoice}(1)]);
P_iter1_feas.plot("alpha",0.2,"Color",[0.8500 0.3250 0.0980],"EdgeColor",0.5*[0.8500 0.3250 0.0980],"Linewidth",0.75);

h_facets = P_iter1_feas.minHRep.getFacet;
han_facet = h_facets(end).minVRep().plot("alpha",0);
hatchfill2(han_facet,"HatchColor",0.75*[0.8500 0.3250 0.0980],"HatchAngle",0)

vert_cut = P_iter0.V(2,:)';
scatter3(vert_cut(1),vert_cut(2),vert_cut(3),25,"MarkerEdgeColor",[0.4660 0.6740 0.1880],"LineWidth",1.5)


Ax1 = gca;
Ax1.YGrid = 'on'; Ax1.Layer = 'top'; Ax1.GridAlpha = 0.05;
Ax1.TickLabelInterpreter = 'latex';

view(-200,20)

% xlabel('$x_1$','Interpreter','latex');
% ylabel('$x_2$','Interpreter','latex');
% zlabel('$x_3$','Interpreter','latex');

%%
% % minimize white borders around plot
fig.PaperPositionMode = "auto";
fig.PaperUnits = "centimeters";
fig.PaperSize = fig.Position(3:4);

fig.Renderer = 'painters'; % a way to force saving in vector graphics

% % save the plot as PDF file
saveas(fig, '../figures/RCI_Cut.pdf','pdf')