function plt = boxplotsADRC(groups,values, varOfInterest, results_APOE, plotTitle, plotXLabel, plotAnnotation, response)

if nargin<8
    response = 'apoe';
end

% Define figure inputs
metRes = results_APOE(matches(results_APOE.Reaction,varOfInterest),{'subgroups','pValue'}); % Get p-values for each regression

% Create boxplot
plt = boxchart(groups,values,'GroupByColor',groups,'BoxWidth',0.8);
fSize = 6;
ylabel(plotXLabel,"FontSize",fSize)
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
if matches(response,'apoe')
    pcn_mci = sprintf('%0.2E',metRes.pValue(matches(metRes.subgroups,"E2-E3")) );
    pmic_dm = sprintf('%0.2E',metRes.pValue(matches(metRes.subgroups,"E3-E4")) );
    pcn_dm = sprintf('%0.2E',metRes.pValue(matches(metRes.subgroups,"E2-E4")) );
else
    pcn_mci = sprintf('%0.2E',metRes.pValue(matches(metRes.subgroups,"normal cognition-MCI")) );
    pmic_dm = sprintf('%0.2E',metRes.pValue(matches(metRes.subgroups,"MCI-Dementia")) );
    pcn_dm = sprintf('%0.2E',metRes.pValue(matches(metRes.subgroups,"normal cognition-Dementia")) );
end

% Find the p-value x coordinates
p1Xcoord = mean(hl1);
p2Xcoord = mean(hl2);
p3Xcoord = mean(hl3);

% Get the p-value y coordinates
offset = 0.025;
p1Ycoord = vl1 + (newYrange*offset);
p2Ycoord = vl2 + (newYrange*offset);
p3Ycoord = vl3 + (newYrange*offset);

% Annotate plot with p-value
fSize = 10;
text(p1Xcoord,p1Ycoord,pcn_mci,'HorizontalAlignment','center','FontName','Arial','FontSize',fSize,'Interpreter','none');
text(p2Xcoord,p2Ycoord,pmic_dm,'HorizontalAlignment','center','FontName','Arial','FontSize',fSize,'Interpreter','none');
text(p3Xcoord,p3Ycoord,pcn_dm,'HorizontalAlignment','center','FontName','Arial','FontSize',fSize,'Interpreter','none');

% Add title
title(plotTitle,'FontWeight','normal','Interpreter','none')
set(gca,'FontSize',12)
set(gca,'TitleFontSizeMultiplier',1.2)

% Annotate plot with plot names
% subPlotNames = 'abcdefghij';
titleProperties = get(gca,'Title');
text(0.0, titleProperties.Position(2),plotAnnotation,'FontSize',12*1.5,'VerticalAlignment','bottom')
end