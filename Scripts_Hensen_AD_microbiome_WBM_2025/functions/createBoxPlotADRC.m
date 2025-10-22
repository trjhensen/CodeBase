function [plt,ax] = createBoxPlotADRC(groups,values, pValTable, plotTitle, yTitle, edgeAlpha, apoeFlag, italicTitle)
% INPUTS:
% groups        Categorical vector with three ordered categories  
% values        Numerical vector with data to visualise
% pValTable     Table with p-values for each pairwise comparison.
% plotTitle     Title for plot
% yTitle        Title of y-axis
% 
% OUTPUT
% Your plot
%
% Author, Tim Hensen, July 2025


% Create boxplot
hold on
if apoeFlag 
    plt = boxchart(groups,values,'GroupByColor',groups,'BoxWidth',1,'MarkerStyle','none','BoxFaceAlpha',0.7,'BoxEdgeColor',[0 0 0]);
    newcolors = ["#66C3A6" "#F68C64" "#8DA0CA"];
    colororder(newcolors)
else
    plt = boxchart(groups,values,'GroupByColor',groups,'BoxWidth',1,'MarkerStyle','none');
end

% Find x coordinates
swarmXcoord = double(groups);
swarmXcoord(swarmXcoord==1) = swarmXcoord(swarmXcoord==1) - 1/numel(categories(groups));
swarmXcoord(swarmXcoord==3) = swarmXcoord(swarmXcoord==3) + 1/numel(categories(groups));

% Add data points
if apoeFlag 
    swarmchart(swarmXcoord,values,12,"black",'MarkerEdgeAlpha',edgeAlpha,'XJitter','density','XJitterWidth',0.5,'MarkerEdgeColor','flat')
else
    swarmchart(swarmXcoord,values,12,'filled','MarkerFaceAlpha',0,'MarkerEdgeAlpha',edgeAlpha,'XJitter','density','XJitterWidth',0.5,'MarkerEdgeColor','flat')
end

fSize = 6;
ylabel(yTitle,"FontSize",fSize,'Interpreter','tex')


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

% Find the pairwise comparisons in pValTable
groupCats = categories(groups);
cat1 = strjoin(groupCats([1 2]),' – ');
cat2 = strjoin(groupCats([2 3]),' – ');
cat3 = strjoin(groupCats([1 3]),' – ');

% Convert p-value to signifance star sign
pCat1 = convertPvalToStar( pValTable.pValue(matches(pValTable.subgroups, cat1)) );
pCat2 = convertPvalToStar( pValTable.pValue(matches(pValTable.subgroups, cat2)) );
pCat3 = convertPvalToStar( pValTable.pValue(matches(pValTable.subgroups, cat3)) );

% formatPval = @(x) append('p = ', sprintf('%0.3f',pValTable.pValue(matches(pValTable.subgroups,x)) ) );
% pCat1 = formatPval(cat1);
% pCat1 = '**';
% pCat2 = formatPval(cat2);
% pCat3 = formatPval(cat3);

% pCat1 = sprintf(strSpec,pValTable.pValue(matches(pValTable.subgroups,cat1)) );
% pCat2 = sprintf(strSpec,pValTable.pValue(matches(pValTable.subgroups,cat2)) );
% pCat3 = sprintf(strSpec,pValTable.pValue(matches(pValTable.subgroups,cat3)) );

% Find the p-value x coordinates
p1Xcoord = mean(hl1);
p2Xcoord = mean(hl2);
p3Xcoord = mean(hl3);

% Get the p-value y coordinates
%offset = 0.025;
offset = 0.03;
if apoeFlag
    offset = 0.025;
end
p1Ycoord = vl1 + (newYrange*offset);
p2Ycoord = vl2 + (newYrange*offset);
p3Ycoord = vl3 + (newYrange*offset);

% Annotate plot with p-value
fSize = 14;
text(p1Xcoord,p1Ycoord,pCat1,'HorizontalAlignment','center','FontName','Arial','FontSize',fSize,'Interpreter','none');
text(p2Xcoord,p2Ycoord,pCat2,'HorizontalAlignment','center','FontName','Arial','FontSize',fSize,'Interpreter','none');
text(p3Xcoord,p3Ycoord,pCat3,'HorizontalAlignment','center','FontName','Arial','FontSize',fSize,'Interpreter','none');

% Add title
if italicTitle == true
    plotTitle{1} = append('{\it ',char(plotTitle{1}),'}');
    title(plotTitle,'FontWeight','normal','Interpreter','tex')
else
    title(plotTitle,'FontWeight','normal','Interpreter','none')
end

set(gca,'FontSize',12)
set(gca,'TitleFontSizeMultiplier',1.2)

% Get axis properties
ax = gca;
end

function star = convertPvalToStar(pval)
thresholds = [0.001, 0.01, 0.05, inf]; % significance levels
labels     = {'***','**','*','ns'}; % corresponding stars
idx = find(pval <= thresholds, 1, 'first'); % first threshold it passes
star = labels{idx};
end