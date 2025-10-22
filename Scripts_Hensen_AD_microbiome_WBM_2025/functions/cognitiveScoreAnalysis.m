function results_G = cognitiveScoreAnalysis(preparedInputTable,preparedMetadata, predictor, confounders, pFilter)
% Perform regressions on global cognitive scores

% Generate subgroups for the three possible comparisons
filters = {{'normal cognition'}, {'MCI'}, {'Dementia'},{'normal cognition', 'MCI'},{'normal cognition', 'MCI', 'Dementia'}};
subgroups = cellfun(@(x) preparedMetadata(matches(preparedMetadata.NACCUDSD, x),:), filters,'UniformOutput',false);

% State response variable of interest for cognitive decline regressions
% response = 'G'; % 
response = 'NACCMOCA';

% Define regression formula
regFormula = string(strcat(response,'~',predictor,'+',strjoin(confounders,'+')));

% Perform logistic regressions for each subgroup
results_G = cellfun(@(x) performRegressions(preparedInputTable,x,regFormula), subgroups,'UniformOutput',false);

% Check if all cells are non-empty and remove the non-empty cells 
cellsToKeep = cellfun(@(x) ~any(matches(fieldnames(x),'NotDefined')), results_G,'UniformOutput',true);
results_G = results_G(cellsToKeep);
filters = filters(cellsToKeep);

% Find flux results and label results with subgroup types
addCohortFun = @(x,y) addvars(x.Flux,repmat( string( strjoin(y,'-') ) ,height(x.Flux), 1), 'Before','N','NewVariableNames','subgroups');
results_G = cellfun(addCohortFun,results_G,filters,'UniformOutput',false);

% Concatinate results into a single table
results_G = vertcat(results_G{:});

% Find groups for each reaction
groups = findgroups(results_G.Reaction);
% Test for each reaction if the results are consistent
groupConsistent = splitapply(@(x) all(x > 0) || all(x < 0), results_G.estimate, groups);
% Add consistency labels to regression results
results_G.consistent = groupConsistent(groups);

% Filter on metabolites with at least on significant result
results_G = groupfilter(results_G,"Reaction",@(x) any(x<pFilter),"pValue");

end