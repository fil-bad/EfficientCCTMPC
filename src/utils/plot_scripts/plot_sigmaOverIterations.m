fig = figure("Renderer","painters","Units","centimeters");
fig.Position(3:4) = 3*[2,2.5]; % [width,height]

han_sigmas = scatter(0:i_max, save_sigma,'filled');

% xlabel('Iterations (i)','Interpreter','latex');
% ylabel('$\sigma^i$','Interpreter','latex');

% ylim([180,224])

Ax1 = gca;
Ax1.YGrid = 'on'; Ax1.Layer = 'top'; Ax1.GridAlpha = 0.05;
Ax1.TickLabelInterpreter = 'latex';

han_leg = legend(Ax1,han_sigmas,{'$\sigma_i$'},'Interpreter','latex','Location','northeast');
han_leg.FontSize = 14;

% % minimize white borders around plot
fig.PaperPositionMode = "auto";
fig.PaperUnits = "centimeters";
fig.PaperSize = fig.Position(3:4);

fig.Renderer = 'painters'; % a way to force saving in vector graphics

saveas(fig, '../figures/sigma_overIters.pdf','pdf')
