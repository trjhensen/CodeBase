function fig = createCognitionForestPlot(filtered_res_G, filtered_res_G_plasma, saveDir)

xTitle = "Regression coefficient"; 
plotLabelSize = 14;

% clc; close all; 
fig = figure('Position',[680,504,1037,374*0.8]);
t = tiledlayout(1,2);

nexttile;
plotTitle = "Predicted blood fluxes"; hideLegend = false;
createForestPlot(filtered_res_G.estimate, filtered_res_G{:,["low","high"]}, filtered_res_G.Reaction, filtered_res_G.pValue, plotTitle, xTitle, hideLegend)

% Change title weight
T = get(gca,'Title');
T.FontWeight = 'normal';
set(gca,'TitleFontSizeMultiplier',1.3)

% Annotate plot with plot names
set(gca,'TitleHorizontalAlignment','center')
text(min(xlim), get(gca,'Title').Position(2),'l','FontSize',plotLabelSize,'VerticalAlignment','bottom','HorizontalAlignment','center')

nexttile;
plotTitle = "Metabolomic plasma abundances"; hideLegend = true;
createForestPlot(filtered_res_G_plasma.estimate, filtered_res_G_plasma{:,["low","high"]}, filtered_res_G_plasma.Reaction, filtered_res_G_plasma.pValue, plotTitle, xTitle, hideLegend)
% Change title weight
T = get(gca,'Title');
T.FontWeight = 'normal';
set(gca,'TitleFontSizeMultiplier',1.3)

% Annotate plot with plot names
set(gca,'TitleHorizontalAlignment','center')
text(min(xlim), get(gca,'Title').Position(2),'m','FontSize',plotLabelSize,'VerticalAlignment','bottom','HorizontalAlignment','center')


tl = title(t, {'','',''},'FontSize',6);

% Get position of the first tile
ax1 = nexttile(1);  % First tile
pos = get(ax1, 'Position');  % Position of first tile in normalized units [left bottom width height]

% Add custom title using 'text' in normalized figure units
annotationX = pos(1);         % x-position aligned with first tile
annotationY = pos(2) + pos(4) + 0.15;  % y-position slightly above the first tile

% Use annotation or uicontrol, or best: use axes('Position', ...) + text
axes('Position', [0 0 1 1], 'Units', 'normalized', 'Visible', 'off'); % invisible overlay axes

fSize = 14;
figTitle = 'Flux and metabolomic associations with the global cognition scores';
text(annotationX, annotationY, figTitle, ...
     'Units', 'normalized', ...
     'HorizontalAlignment', 'left', ...
     'FontSize', fSize, ...
     'FontWeight', 'bold');

% Save figure 
% saveDir = paths.figures;
exportgraphics(fig,fullfile(saveDir,'global_cognition_FP.png'),'Resolution',450)
end