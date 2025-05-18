fig = figure("Renderer","painters","Units","centimeters");
fig.Position(3:4) = 3*[2,2.5]; % [width,height]
hold on

X.plot("wire",true,"linestyle",'--',"edgecolor",0.5*ones(1,3),"edgealpha",0.5)

han_best = Polyhedron(save_Triple{end}{1}, save_yM{end}).plot(...
    "alpha",0.1,"Color",[0.4660 0.6740 0.1880],"EdgeColor",0.75*[0.4660 0.6740 0.1880],"Linewidth",0.75);
han_start = Polyhedron(save_Triple{1}{1}, save_yM{1}).plot(...
    "alpha",0.2,"Color",[0.8500 0.3250 0.0980],"EdgeColor",0.75*[0.8500 0.3250 0.0980],"Linewidth",0.75);


Ax1 = gca;
Ax1.YGrid = 'on'; Ax1.Layer = 'top'; Ax1.GridAlpha = 0.05;
Ax1.TickLabelInterpreter = 'latex';

view(-200,20)

han_leg = legend(Ax1,[han_start, han_best],...
    {'Initial Set','Final Set'}, ...
    'Interpreter','latex','Location','northeast');
han_leg.FontSize = 10;

% % minimize white borders around plot
fig.PaperPositionMode = "auto";
fig.PaperUnits = "centimeters";
fig.PaperSize = fig.Position(3:4);

fig.Renderer = 'painters'; % a way to force saving in vector graphics

% % save the plot as PDF file
saveas(fig, "../figures/RCI_afterIters.pdf",'pdf')