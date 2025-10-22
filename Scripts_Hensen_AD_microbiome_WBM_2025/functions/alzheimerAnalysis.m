function [results_AD, rxnsOfInterest, results_DM, regressions_DM] = alzheimerAnalysis(preparedInputTable,preparedMetadata,predictor,confounders)
% input:
% preparedInputTable
% preparedMetadata
% predictor
% confounders

% output: 
% AD results
% AD cognitive decline results
% rxnsOfInterest


% 1) Associate fluxes with AD status

% State response variable of interest
response = 'AD';

% Define regression formula
ADForm = strcat(response,'~',predictor,'+',strjoin(confounders,'+')); 

% Perform regressions
results_AD = performRegressions(preparedInputTable,preparedMetadata,ADForm); 
results_AD = results_AD.Flux;

% Find metabolites that associate with AD.
rxnsOfInterest = results_AD.Reaction(results_AD.FDR<0.05);

% 2) Associate fluxes with cognitive decline from control to MCI to
% dementia.

% Prune input table with reactions of interest
prunedInputTable = preparedInputTable;%(:,["ID",rxnsOfInterest']);

% Generate subgroups for the three possible comparisons
filters = { {'Healthy controls','MCI'}, {'MCI','AD dementia'}, {'Healthy controls','AD dementia'} };
subgroups = cellfun(@(x) preparedMetadata(matches(preparedMetadata.NACCUDSD, x),:), filters,'UniformOutput',false);

% State response variable of interest for cognitive decline regressions
response = 'NACCUDSD';

% For each subgroup, encode the less severe group as 0 and the more severe
% group as 1
convertBinaryResp = @(x,y) grp2idx(categorical(x.(response),y,'Ordinal',true)); % Changed!
addBinaryResp = @(x,y) addvars(x,convertBinaryResp(x,y), 'After',response,'NewVariableNames','status');
subgroups = cellfun(addBinaryResp, subgroups,filters, 'UniformOutput',false);

subgroups = cellfun(@(x) removevars(x, response), subgroups,'UniformOutput',false);
subgroups = cellfun(@(x) renamevars(x, 'status',response), subgroups,'UniformOutput',false);

% Define regression formula
cDeclForm = string(strcat(response,'~','Flux+',strjoin(confounders,'+')));

% Perform logistic regressions for each subgroup
[results_DM,regressions_DM] = cellfun(@(x) performRegressions(prunedInputTable,x,cDeclForm), subgroups,'UniformOutput',false);

% Find flux results and label results with subgroup types
addCohortFun = @(x,y) addvars(x.Flux,repmat( string( strjoin(y,' – ') ) ,height(x.Flux), 1), 'Before','N','NewVariableNames','subgroups');
results_DM = cellfun(addCohortFun,results_DM,filters,'UniformOutput',false);
% Concatinate results into a single table
results_DM = vertcat(results_DM{:});

% Check if the reactions changed consistently: normal -> MCI -> Dementia

% Find groups for each reaction
groups = findgroups(results_DM.Reaction);
% Test for each reaction if the results are consistent
groupConsistent = splitapply(@(x) all(x > 0) || all(x < 0), results_DM.estimate, groups);
% Add consistency labels to regression results
results_DM.consistent = groupConsistent(groups);

% Sort on p-value
results_DM = sortrows(results_DM,'pValue','ascend');

rxnsOfInterest = unique(results_DM.Reaction(results_DM.pValue<0.05),'stable');

end