function fig = createDementiaBoxplots(rxnsOfInterest, preparedInputTable, preparedMetadata, resultsTable, saveDir)
% Visualise dementia regression resultsTable for the reactions of interest
% input:
% rxnsOfInterest
% preparedInputTable
% preparedMetadata
% resultsTable

% Test if the input data contains metabolite names or microbe names

% Generate boxplot
fig = figure('Position',[147,432,1512,445*1.1]);
t = tiledlayout(2,numel(rxnsOfInterest),'TileSpacing','tight','Padding','loose');
type = "blood fluxes";
metNames = renameAdrcVmhToMetName(rxnsOfInterest); % Get metabolite names

for i=1:numel(metNames)
    nexttile;
    metName = metNames(i);
    dementiaBoxPlots(preparedMetadata,preparedInputTable, type, rxnsOfInterest, resultsTable, metName, i);
end

fSize = 14;
%ylabel(t,"Normalised log2 flux in mmol/day/person",'FontSize',fSize)
title(t,["Healthy controls (CN) vs Mild cognitive impairment (MCI) vs Dementia patients (DEM)",""],'FontSize',fSize,'HorizontalAlignment','center')

% Save figure 
exportgraphics(fig,fullfile(saveDir,'flux_DM_boxplots.png'),'Resolution',300)
end

function plt = dementiaBoxPlots(preparedMetadata,preparedInputTable, type, rxnsOfInterest, resultsTable, metName, i)

% Define figure inputs
group = preparedMetadata.NACCUDSD;%categorical(preparedMetadata.NACCUDSD,{'normal cognition','MCI','Dementia'});
values = preparedInputTable.(rxnsOfInterest(i));
metRes = resultsTable(matches(resultsTable.Reaction,rxnsOfInterest(i)),{'subgroups','pValue'}); % Get p-values for each regression


% Create boxplot
%plt = boxplot(values,group);
plt = boxchart(group,values,'GroupByColor',group,'BoxWidth',0.8);
ylabel("Flux in mmol/day/person")
hold on

% Increase vertical space by x%
rIncr = 0.2; % Percentage increase
yRange = max(ylim) - min(ylim);
ax = gca;
ax.YLim(2) = max(ylim) + (yRange*rIncr);

% Add black horizontal lines to figure

% Define horizontal coordinates for the black lines
horlen = 2/3; 
midPoint = 2;
endPoint = midPoint + (midPoint-horlen);

hl1 = [0+horlen,midPoint]; % Left line
hl2 = [midPoint,endPoint]; % Right line
hl3 = [horlen, endPoint]; % Top line

% Define vertical coords
newYrange = max(ylim) - min(ylim);

vl1 = max(ylim) - (newYrange*0.15); % 85% from the top
vl2 = max(ylim) - (newYrange*0.1); % 90% from the top
vl3 = max(ylim) - (newYrange*0.05); % 95% from the top

% Add black lines
plot(hl1,[vl1,vl1], '-k', 'LineWidth',1) % Left line
plot(hl2,[vl2,vl2], '-k', 'LineWidth',1) % Right line
plot(hl3,[vl3,vl3], '-k', 'LineWidth',1) % Top line

% Add p-values

% Find the p-values 
pcn_mci = sprintf('%0.2E',metRes.pValue(matches(metRes.subgroups,"normal cognition-MCI")) );
pmic_dm = sprintf('%0.2E',metRes.pValue(matches(metRes.subgroups,"MCI-Dementia")) );
pcn_dm = sprintf('%0.2E',metRes.pValue(matches(metRes.subgroups,"normal cognition-Dementia")) );

% Find the p-value x coordinates
p1Xcoord = mean(hl1);
p2Xcoord = mean(hl2);
p3Xcoord = mean(hl3);

% Get the p-value y coordinates
offset = 0.02;
p1Ycoord = vl1 + (newYrange*offset);
p2Ycoord = vl2 + (newYrange*offset);
p3Ycoord = vl3 + (newYrange*offset);

% Annotate plot with p-value
fSize = 10;
text(p1Xcoord,p1Ycoord,pcn_mci,'HorizontalAlignment','center','FontName','Arial','FontSize',fSize,'Interpreter','none');
text(p2Xcoord,p2Ycoord,pmic_dm,'HorizontalAlignment','center','FontName','Arial','FontSize',fSize,'Interpreter','none');
text(p3Xcoord,p3Ycoord,pcn_dm,'HorizontalAlignment','center','FontName','Arial','FontSize',fSize,'Interpreter','none');

% Realine xticks
% set(gca,'XTick',[0,horlen,midPoint,endPoint])

% Add group names to xtick labels
% ax.XTickLabel(1) = {'DM'}; % WHY DOES THIS WORK?
% ax.XTickLabel(3) = {'MCI'}; % WHY DOES THIS WORK?
% ax.XTickLabel(2) = {'CN'}; % WHY DOES THIS WORK?

% Add title
title({strcat("Predicted ", metName), type},'FontWeight','normal')
set(gca,'FontSize',12)
set(gca,'TitleFontSizeMultiplier',1.2)

% Annotate plot with plot names
subPlotNames = 'abcdefghij';
titleProperties = get(gca,'Title');
text(0.0, titleProperties.Position(2),subPlotNames(i),'FontSize',12*1.5,'VerticalAlignment','bottom')
end