function [prunedData, linearDependentVars, filePath] = findLinearlyDependentVariables(inputTable, tolerance, saveDir)
% findLinearlyDependentVariables
%
% USAGE:
%   [prunedData, linearDependentVars, filePath] = findLinearlyDependentVariables(inputTable, tolerance, saveDir)
%
% DESCRIPTION:
%   Identifies and removes any columns in the input table that are
%   (approximately) linear combinations of other columns, within a
%   specified numerical tolerance.
%
% INPUTS:
%   inputTable           [table]   N×m table of variables in columns; sample IDs must
%                                 reside in inputTable.Properties.RowNames (first
%                                 column must not contain the IDs).
%   tolerance            [double]  Threshold for declaring linear dependence:
%                                 ‖v1 – λ·v2‖ ≤ tolerance. Recommended: 1e-6.
%   saveDir              [char]    Path to directory in which to write the Excel file.
%
% OUTPUTS:
%   prunedData           [table]   Subset of inputTable retaining only linearly
%                                 independent columns.
%   linearDependentVars  [table]   Table summarizing which columns were found to be
%                                 dependent and their dependency coefficients.
%   filePath             [char]    Full path to the generated .xlsx file containing
%                                 prunedData and linearDependentVars.
%
% AUTHOR:
%   Tim Hensen (June 2025)

% Create data matrix
dataMatrix = table2array(inputTable);

% Perform a pearson correlation on the reaction abundances to find
% potential linear dependencies.
RHO = abs(corr(dataMatrix,'Rows','pairwise'));

% Set the diagonal to zero
RHO(1:size(RHO,1) + 1:end) = 0;

% % Find all potentially linear dependent pairs by searching for the row
% and column indices of all reaction pairs with a correlation above the
% threshold
corrThreshold = 0.99;
[row,col] = find(RHO>corrThreshold);

% This loop quantifies linear dependence by calculating the euclidean
% norm x of the difference between v1 and v2 times lambda, the ratio between 
% v2 and v1, i.e., x = ||v1 - lambda*v2||. Note that v1 and v2 are
% linearly dependent if ||v1 - lambda*v2|| = 0 

% Preallocate table for quantifying linear dependence
numVars = size(dataMatrix,2);
rxnNames = inputTable.Properties.VariableNames;
linDepTable = array2table(nan(numVars,numVars),'RowNames',rxnNames','VariableNames',rxnNames);
fprintf("> Searching for linearly dependent variables... \n")
tic 
numPairs = length(row);
for i=1:numPairs

    % Set up function input:
    vector1 = dataMatrix(:,row(i));
    vector2 = dataMatrix(:,col(i));

    % Normalise the vectors to make them robust to scale differences. 
    norm_vector1 = vector1 / norm(vector1);
    norm_vector2 = vector2 / norm(vector2);
    
    % Calculate lambda by solving a linear system of equations A*lambda = B
    lambda = norm_vector2 \ norm_vector1; 
    
    % Calculate the euclidean norm
    linNorm = norm(norm_vector1 - lambda * norm_vector2);
    
    % Add linear norm to table
    linDepTable{row(i),col(i)} = linNorm;

    if mod(i, floor(numPairs/100)) == 0
        disp(append("> ",string(round(i/numPairs*100,0)),"%"))
    end
end
toc

fprintf("> Processing results... \n")

% Apply threshold to identify linearly dependent data
adjMatrix = table2array(linDepTable);

% Find identical distributions
adjMatrix(adjMatrix<tolerance) = 1;
adjMatrix(adjMatrix~=1) = 0;

% Filter on linearly dependent metabolites
colsAndRowsToRemove = all(adjMatrix==0); % Note that the columns and rows are identical

% Find linearly dependent networks
linDepTable1 = array2table(adjMatrix,'RowNames',rxnNames,'VariableNames',rxnNames);
linDepTable1(colsAndRowsToRemove',:) = [];
linDepTable1(:,colsAndRowsToRemove) = [];

% Create graph and analyze connectivity of adjacency matrix
G = graph(table2array(linDepTable1), linDepTable1.Properties.VariableNames,'upper');
nodeDegrees = degree(G);

% Find connected components (subnetworks)
subnetworks = conncomp(G);

% Create subnetwork list with relevant information
subnetwork_list = table(subnetworks', G.Nodes.Name, nodeDegrees, ...
    'VariableNames', {'Subnetwork', 'Reaction', 'Node degree'});

% Find all linearly dependent reactions by removing the nodes with a degree
% of zero.
subnetwork_list = subnetwork_list(subnetwork_list.('Node degree') ~= 0, :);
%subnetwork_list.('Node degree') = [];
subnetwork_list1 = subnetwork_list;


% Load reaction subsystem information
rxns = cell2table(loadVMHDatabase().reactions(2:end,[1 11 12]),'VariableNames',{'Reaction','Subsystem','General subsystem'});
% Combine subsystem information with reaction information
subnetwork_list1 = outerjoin(subnetwork_list1,rxns,'Keys','Reaction','MergeKeys',true,'Type','left');

% Create output table
linearDependentVars = sortrows(subnetwork_list1,'Subnetwork');

%%% TO REMOVE?
    % First, we will add the number of non-nan values for each metabolite
    % numZeros = varfun(@(x) sum(x==0), inputTable);
    % numZeros.Properties.VariableNames = erase(numZeros.Properties.VariableNames,'Fun_');
    % numZeros = stack(numZeros,1:width(numZeros),'NewDataVariableName','numZeros','IndexVariableName','Reaction');
    % numZeros.Reaction = cellstr(numZeros.Reaction);
    % linearDependentVars = outerjoin(linearDependentVars,numZeros,'MergeKeys',true,'Type','left');
    % linearDependentVars = sortrows(linearDependentVars,'Subnetwork'); % Sort on subnetwork name
%%% 

% I can actually use the clusters to reconstruct the union of data results
% for the entire cluster. This can be done by 1) finding the mapping
% function for each reaction pair in the cluster and normalising the inputs
% according the to the vector with the lowest values. 2) Generate a new
% data vector that imputes the missing values with the scaled values of
% linearly dependent reactions. 

% Find the subnetworks in the input table
linDepFluxes = struct;
for i=1:max(linearDependentVars.Subnetwork)
    rxnsToCheck = linearDependentVars.Reaction(linearDependentVars.Subnetwork==i);
    linDepFluxes.(strcat("subnetwork_",string(i))) = inputTable(:,rxnsToCheck');
end

% Now, I want to scale all values in the linearly dependent matrices.
% First, I want to find the scaling factor for each vector and apply that.

% Find the number of subnetworks
numSubnets = max(linearDependentVars.Subnetwork); 

% Preallocate table for the combined linearly dependent results
constructedClusterFluxes = table('Size', ...
    [height(inputTable),1+numSubnets], ...
    'VariableTypes',[{'string'} repmat({'double'},1,numSubnets)], ...
    'VariableNames',[{'ID'}, cellstr(arrayfun(@(i) append("subnetwork_",string(i)),1:numSubnets))]);

% Add ID information
constructedClusterFluxes.ID = inputTable.Properties.RowNames;

% Add custom table property that includes the associated reactions
constructedClusterFluxes = addprop(constructedClusterFluxes,{'Reactions'},{'variable'}); 
constructedClusterFluxes.Properties.CustomProperties.Reactions = cell(width(constructedClusterFluxes),1);

for i=1:numel(fieldnames(linDepFluxes))

    % Take a subset of linearly dependent inputs and create a matrix
    linDepTest = linDepFluxes.(strcat("subnetwork_",string(i)));
    linDepMatrix = table2array(linDepTest);
    
    % Remove rows with NaN values in the matrix
    linDepMatrixNoNaN = linDepMatrix(all(~isnan(linDepMatrix),2),:);
    
    % Find the column with the smallest data values. Use this column as a base.
    [~, min_index] = min(median(linDepMatrixNoNaN));
    min_col = linDepMatrixNoNaN(:,min_index);
    
    % Obtain the scaling factors compared to the column with the smallest
    % values. Note that this only works because the reactions have linearly
    % dependent inputs.
    scalingFactors = arrayfun(@(x) linDepMatrixNoNaN(:,x) \ min_col, 1:width(linDepMatrixNoNaN));
    
    % Apply the scaling factors to each column. This operation will scale the
    % values of each reaction to that of the reaction with the smallest values.
    scaled_linDepMatrix = linDepMatrix .* scalingFactors;
    
    % Now, we can perform imputation of all NaN values in the matrix column
    % with the smallest values
    
    % Create new vector for integrated reaction inputs with imputed missing
    % values.
    clusterCol = linDepMatrix(:,min_index);
    
    % Find the rows with nan values
    naRows = find(isnan(clusterCol));
    
    % Check if any of the other column contain a non-nan values. If yes, then
    % impute the scaled values from the other columns in clusterCol.
    for j=1:length(naRows)
        nonNanCol = find(~isnan(scaled_linDepMatrix(naRows(j),:)));
        if ~isempty(nonNanCol)
            clusterCol(naRows(j)) = scaled_linDepMatrix(naRows(j),nonNanCol(1));
        end
    end
    
    % Add the imputed vector to the linDepFluxes variable
    linDepFluxes.(strcat("subnetwork_",string(i))).Imputed = clusterCol;

    % Also add the imputed vector to a table with the constructed cluster
    % inputs.
    constructedClusterFluxes.(strcat("subnetwork_",string(i))) = clusterCol;

    % Store the associated reactions in each column
    constructedClusterFluxes.Properties.CustomProperties.Reactions{i+1} = linDepTest.Properties.VariableNames;
end

% Now, we can remove all linearly dependent metabolites and replace them by
% the constructed clusters

% Create pruned inputs table
prunedData = removevars(inputTable, subnetwork_list1.Reaction);
prunedData = addvars(prunedData,inputTable.Properties.RowNames,'NewVariableNames','ID','Before',1);

% Add back the constructed inputs
prunedData = outerjoin(prunedData,constructedClusterFluxes,'MergeKeys',true);

% Clarify subnetwork names
linearDependentVars.Subnetwork = append("subnetwork_",string(linearDependentVars.Subnetwork));


fprintf("> Save results to disk... \n")
filePath = fullfile(saveDir,'linDepDataSubsets.xlsx');

if isfile(filePath)
    delete(filePath)
end

% Save results
writetable(linearDependentVars, filePath,'Sheet','linDepVars');
writetable(prunedData, filePath,'Sheet','uniqueData');

figure;
p = plot(G);
p.Interpreter = 'none';
title('Linearly dependent subnetworks')
end