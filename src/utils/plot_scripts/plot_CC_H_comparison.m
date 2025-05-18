
fig = figure("Renderer","painters","Units","centimeters");
fig.Position(3:4) = 3*[3,1.5]; % [width,height]
hold on

% --- Plot data for left y-axis ---
yyaxis left % Activate left y-axis
% Here we use a light blue dashed line to connect the points,
% while the markers remain a full blue color.
haus_han = plot(N_ocp+horIncr, hausH(1:length(horIncr)), ...
    'LineStyle', ':', ...                         % Dashed connecting line
    'LineWidth', 0.2, ...                         % Very thin line
    'Marker', 'o', ...                            % Marker style: circle
    'MarkerFaceColor', [0 0.4470 0.7410], ...     % Filled markers in blue
    'MarkerSize', 5,...
    'Color', 0.75*[0 0.4470 0.7410]);             % Light blue connecting line
% Draw horizontal reference line for left y-axis using the same light blue
yline(baseHausCC, '-', "$d_\mathcal{X}(\mathcal{O}_c(3))$",...
    'Color', 0.75*[0 0.4470 0.7410],'Interpreter','latex');
set(gca, 'YColor',  0.75*[0 0.4470 0.7410],'TickLabelInterpreter', 'latex'); % Left axis tick labels in blue
ylim([13.3,14.65])

% --- Plot data for right y-axis ---
yyaxis right                  % Switch to right y-axis
set(gca, 'YScale', 'log')
% Similar approach: light red dashed line with fully-colored red markers.
time_han = plot(N_ocp+horIncr, timeH_s(1:length(horIncr)), ...
    'LineStyle', ':', ...                         % Dashed connecting line
    'LineWidth', 0.2, ...                           % Very thin line
    'Marker', 'square', ...                              % Marker style: circle
    'MarkerFaceColor', [0.8500 0.3250 0.0980], ...                     % Filled markers in blue
    'MarkerSize', 5,...
    'Color', 0.75*[0.8500 0.3250 0.0980]);                          % Light blue connecting line
% Draw horizontal reference line for right y-axis using the same light red,
% and position its label to the left instead of at the end.
yline(timeCC, '-', "CCTMPC$(3)$", 'Color', 0.75*[0.8500 0.3250 0.0980], 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'top','Interpreter','latex');
set(gca, 'YColor', 0.75*[0.8500 0.3250 0.0980],'TickLabelInterpreter', 'latex');                           % Right axis tick labels in red
ylim([5,500])

% Common settings for the entire plot
% xlabel('Prediction Horizon HTMPC($N$)','Interpreter','latex');

% --- Set x-axis ticks to display each integer ---
xVals = N_ocp + horIncr;
xMin = floor(min(xVals));
xMax = ceil(max(xVals));
xticks(xMin:1:xMax);  % Create an integer tick for every number in the range

xlim([3-0.1,10+0.1])

leg_han = legend([haus_han,time_han],{'$d_\mathcal{X}(\mathcal{O}_h(N))$','HTMPC$(N)$ [ms]'}, ...
    'Interpreter','latex','Location','southwest');

% leg_han.FontSize = 10;
leg_han.Position(1:2) = [0.15,0.35];
leg_han.Box = "off";

Ax1 = gca;
% % minimize white borders around plot
% set(Ax1,'LooseInset', max(get(Ax1,'TightInset'), 0.01)) % remove border from axis
set(fig,'PaperPositionMode','Auto','PaperUnits',...
    'centimeters','PaperSize',1.1*fig.Position(3:4)) % resize pdf page

fig.Renderer = 'painters'; % a way to force saving in vector graphics
% % save the plot as PDF file
saveas(fig, '../figures/HausTimeComparison.pdf','pdf')


