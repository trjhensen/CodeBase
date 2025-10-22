% performADRCanalysis
% Goal: To investigate associations between the predicted and cognitive
% decline. 
% Is there a difference between AD and controls?
% Are these results unique to AD status or could they be explained by dementia status?
% Is there a difference between dementia and controls?

% Clean up workspace
clearvars -except paths; clc; close all;

% Prepare flux and metadata for further analysis
[preparedInputTable, preparedMetadata] = prepareDataForStatsADRC(paths.fluxPath, paths.metadata);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Is there a difference between AD and controls? %

response = 'AD'; % State response variable of interest
predictor = 'Flux'; % State predictor variable of interest
confounders = {... % State covariates to control for
    'Sex',... % Male/Female    
    'age_at_collection',... % Age in years
    'mapped_species_reads',...
    'lane',... % Sequencing run lane
    'Ethanol_added',...
    'HYPERT',... % Hypertension present
    'DAILY_ALCOHOL'...% Daily alcohol consumption
    };

% Define regression formula
regFormula = strcat(response,'~',predictor,'+',strjoin(confounders,'+')); 

% Perform regressions
[results,regressions] = performRegressions(preparedInputTable,preparedMetadata,regFormula); 

% Find metabolites to plot
rxnsToPlot = results.Flux.Reaction(results.Flux.FDR<0.05);


% Visualise results
fig = generateAdBoxPlots(regressions,rxnsToPlot);

% Save figure
exportgraphics(fig,fullfile(paths.figures,'flux_AD_boxplots.png'),'Resolution',300)

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Is there a difference between Dementia, MCI, and
% and controls? %

% Associate the predicted fluxes with changes in fluxes from normal cognition to mild cognitively impaired to full dementia. 
preparedInputTable = preparedInputTable(:,["ID",rxnsToPlot']);

% Remove samples with non-AD MCI and non-AD Dementia
toRm = find(matches(preparedMetadata.NACCUDSD,'MCI') & matches(preparedMetadata.AD,'no_AD'));
toRm = [toRm; find(matches(preparedMetadata.NACCUDSD,'Dementia') & matches(preparedMetadata.AD,'no_AD'))];
preparedInputTable(toRm,:) = [];
preparedMetadata(toRm,:) = [];

[results, regressions] = performCogDeclineAnalysis(preparedInputTable, preparedMetadata);


% Generate boxplot
fig = figure('Position',[147,432,1512,445]);
t = tiledlayout(1,4,'TileSpacing','tight','Padding','loose');
type = " fluxes in blood";
metNames = renameAdrcVmhToMetName(rxnsToPlot); % Get metabolite names

for i=1:numel(metNames)
    nexttile;
    metName = metNames(i);
    dementiaBoxPlots(preparedMetadata,preparedInputTable, type, rxnsToPlot, results, metName, i);
end

fSize = 14;
ylabel(t,"Normalised log2 flux in mmol/day/person",'FontSize',fSize)
xlabel(t,"Healthy controls (CN) vs AD-Mild cognitive impairment (AD-MCI) vs AD-Dementia patients (AD-DM)",'FontSize',fSize)

exportgraphics(fig,fullfile(paths.figures,'flux_DM_boxplots.png'),'Resolution',300)


function fig = generateAdBoxPlots(regressions,rxnsToPlot)


% Filter fitted regression structures
matlabRxnsToPlot = matlab.lang.makeValidName(rxnsToPlot);
regressionFields = string(fieldnames(regressions));
regressionsToPlot = rmfield(regressions, setdiff(regressionFields, matlabRxnsToPlot ) );

regressionFields = string(fieldnames(regressionsToPlot));
[~,~,ib] = intersect(matlabRxnsToPlot,regressionFields,'stable');
regressionFields = regressionFields(ib);
metNames = renameAdrcVmhToMetName(regressionFields);

% Rearrange metabolite names
matlab.lang.makeValidName(rxnsToPlot)

% Generate boxplot
fig = figure('Position',[147,432,1512,445]);
t = tiledlayout(1,4,'TileSpacing','tight','Padding','loose');
type = " fluxes in blood";
xtickLab = {'CN','AD'};

for i=1:numel(metNames)
    nexttile;
    reg = regressionsToPlot.(regressionFields(i));
    metName = metNames(i);
    createBoxplot(reg,metName,type,xtickLab, i);
end

fSize = 14;
ylabel(t,"Normalised log2 flux in mmol/day/person",'FontSize',fSize)
xlabel(t,"Healthy controls (CN) vs Alzheimer's disease (AD) patients",'FontSize',fSize)
end


function plt = createBoxplot(reg,metName,type,xtickLab, i)

% Create grouped violin for metabolite i
group = categorical(reg.Variables.AD);
plt = boxchart(group,reg.Variables.Flux,'GroupByColor',group,'BoxWidth',0.8);

% TODO: Update colours of violin plots

% Add axis labels and title
title({metName, type},'FontWeight','normal')
%ylab = 'Normalised log2 flux in mmol/day/person';
ylabel('')
xlabel('')
xticklabels(xtickLab)

% Format axis labels
set(gca,'FontSize',12)
set(gca,'TitleFontSizeMultiplier',1.2)
if matches(type," associated microbial abundances")
    set(gca,'TitleFontSizeMultiplier',1.1)
end

% Add regression p-values to plot
pVal = reg.Coefficients{'Flux','pValue'};

% Setup plot annotation with p-value
textInput = ['p = ',sprintf('%0.2E',pVal)];
xcoord = 1.5;
%ycoord = max(fluxesToVis.(metsToTest(i,1))) * 0.9;
yRange = max(ylim) - min(ylim);
ycoord = max(ylim) - yRange/25;
fSize = 14;

% Annotate plot with p-value
text(xcoord,ycoord,textInput,'HorizontalAlignment','center','FontName','Arial','FontSize',fSize,'Interpreter','none');

% Annotate plot with plot names
subPlotNames = 'ABCDEFGHIJKLMNOP';
titleProperties = get(gca,'Title');
text(0.0, titleProperties.Position(2),subPlotNames(i),'FontSize',12*1.5,'VerticalAlignment','bottom')
end

function plt = dementiaBoxPlots(preparedMetadata,preparedInputTable, type, rxnsToPlot, results, metName, i)

% Define figure inputs
group = categorical(preparedMetadata.NACCUDSD,flip({'normal cognition','MCI','Dementia'}));
values = preparedInputTable.(rxnsToPlot(i));
metRes = results(matches(results.Reaction,rxnsToPlot(i)),{'subgroups','pValue'}); % Get p-values for each regression


% Create boxplot
plt = boxchart(double(group),values,'GroupByColor',group,'BoxWidth',0.8);
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
set(gca,'XTick',[0,horlen,midPoint,endPoint])

% Add group names to xtick labels
ax.XTickLabel(1) = {'AD-DM'}; % WHY DOES THIS WORK?
ax.XTickLabel(3) = {'AD-MCI'}; % WHY DOES THIS WORK?
ax.XTickLabel(2) = {'CN'}; % WHY DOES THIS WORK?

% Add title
title({metName, type},'FontWeight','normal')
set(gca,'FontSize',12)
set(gca,'TitleFontSizeMultiplier',1.2)

% Annotate plot with plot names
subPlotNames = 'EFGHIJKLMNOP';
titleProperties = get(gca,'Title');
text(0.0, titleProperties.Position(2),subPlotNames(i),'FontSize',12*1.5,'VerticalAlignment','bottom')
end