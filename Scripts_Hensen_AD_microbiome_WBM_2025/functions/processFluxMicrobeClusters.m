function [topCorrelations,fileName] = processFluxMicrobeClusters(bestCorrPerSizeList, minbestCorrPerSizeList, saveDir)
% This function ensures that each cell array in the input tables has
% identical width and table names. Then, the function concatinates all
% tables in a single tall table.

bestCorrPerSizeList = processClustCorrTabs(bestCorrPerSizeList);
minbestCorrPerSizeList = processClustCorrTabs(minbestCorrPerSizeList);

% Create a single table with all top correlations
topCorrelations = vertcat(bestCorrPerSizeList{:},minbestCorrPerSizeList{:});
topCorrelations = removevars(topCorrelations,'clusterIndex');

% Save results to table
fileName = fullfile(saveDir, 'top_microbialClusterCorrelates.csv');
writetable(topCorrelations,fileName)
end


function topCorrClust = processClustCorrTabs(clustCorrCell)
% Ensure that all tables have the same width

% Find the widest table in the cell array
[~,tabIdx] = max(cellfun(@width, clustCorrCell));

% Identify the variables that miss in the smaller tables in the cell array
fun = @(x) setdiff(clustCorrCell{tabIdx}.Properties.VariableNames, x.Properties.VariableNames);
colNamesAdd = cellfun(fun, clustCorrCell,'UniformOutput',false);

% Get the index of the non-max tables 
idxToChange = ~cellfun(@isempty,colNamesAdd); 

% Add the variables Microbe_X1,... Microbe_XN to each of the tables in the
% cell array with the non-max widths 
clustCorrCell(idxToChange) = cellfun(@iterateAddVars, clustCorrCell(idxToChange), colNamesAdd(idxToChange),'UniformOutput',false);
topCorrClust = clustCorrCell; % Define output
end

function output = iterateAddVars(input, names)
for i=1:length(names)
    varName = names(i);
    input = addvars(input,repmat("",height(input),numel(varName)),'Before','Cluster_size','NewVariableNames',varName);
end
output = input;
end
