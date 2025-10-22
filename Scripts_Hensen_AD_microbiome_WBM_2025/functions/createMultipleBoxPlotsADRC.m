function [plt,ax] = createMultipleBoxPlotsADRC(preparedInputTable, preparedMetadata, rxnsOfInterest, regressionTable, param)
% Generate boxplots for APOE or AD progression results.

% Unpack parameters
response = param.response;
titleAnnotation = param.titleAnnotation;
yTitle = param.yTitle;

if ~any(matches(fields(param),'apoeFlag'))
    apoeFlag = 0;
else
    apoeFlag = param.apoeFlag;
end

if ~any(matches(fields(param),'edgeAlpha'))
    edgeAlpha = 0.4;
else
    edgeAlpha = param.edgeAlpha;
end

if ~any(matches(fields(param),'italicTitle'))
    italicTitle = false;
else
    italicTitle = param.italicTitle;
end





% Set default parameters:
tileAnnotationFSize = 18;

% Merge metadata with data table
mergedData = outerjoin( preparedInputTable(:,["ID", rxnsOfInterest']), preparedMetadata(:,["ID", response]),'MergeKeys',true );

% Get metabolite names
metNames = renameAdrcVmhToMetName(rxnsOfInterest); 
metNames = replace(metNames,'_',' ');

% Get groups to visualise
groups = mergedData.(response);

for i=1:length(rxnsOfInterest)

    % Add next tile
    tile = nexttile;
    
    % Get values from data
    values = mergedData.(rxnsOfInterest(i));
    
    % Get regression p values for metabolite of interest
    pValTable = regressionTable( matches(regressionTable.Reaction,rxnsOfInterest(i)), : ); 
    
    % Generate title for plot
    plotTitle = {metNames(i), titleAnnotation};
    
    % Create plot
    [plt,ax] = createBoxPlotADRC(groups,values, pValTable, plotTitle, yTitle, edgeAlpha, apoeFlag, italicTitle);
    
    % Annotate plot with plot names
    text(0.0, ... x-axis position
        ax.Title.Position(2),... y-axis position
        char(tilenum(tile) + 64),... % plot label letter in alphabet corresponding with current tile
        'FontSize',tileAnnotationFSize,...
        'VerticalAlignment','bottom')
end

% Hack for visualisation
if isfield(param,'addEmptyTile')

    if param.addEmptyTile == true
        tile = nexttile;
        ax = gca;
        title({'','','',''})
        % Annotate plot with plot names
        text(0.0, ... x-axis position
            ax.Title.Position(2),... y-axis position
            char(tilenum(tile) + 64),... % plot label letter in alphabet corresponding with current tile
            'FontSize',tileAnnotationFSize,...
            'VerticalAlignment','bottom')
    end
end

end