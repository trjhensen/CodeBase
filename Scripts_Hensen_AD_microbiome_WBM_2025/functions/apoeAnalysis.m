function [results_APOE,regressions] = apoeAnalysis(preparedInputTable,preparedMetadata, confounders, pFilter)

% Define regression formula
response = 'Flux';

% Only investigate cognitively normal samples
% preparedMetadata = preparedMetadata(matches(preparedMetadata.NACCUDSD,'normal cognition'),:);


% Generate subgroups for the three possible comparisons
filters = { {'ϵ2','ϵ3'}, {'ϵ3','ϵ4'}, {'ϵ2','ϵ4'} };
subgroups = cellfun(@(x) preparedMetadata(matches(preparedMetadata.APOE_ALLELE, x),:), filters,'UniformOutput',false);

% For each subgroup, encode the less severe group as 0 and the more severe
% group as 1
predictor = 'APOE_ALLELE';
convertBinaryResp = @(x,y) grp2idx(categorical(x.(predictor),y,'Ordinal',true));%grp2idx(categorical(x.(predictor),y,'Ordinal',true)); % Changed!
addBinaryResp = @(x,y) addvars(x,convertBinaryResp(x,y), 'After',predictor,'NewVariableNames','status');
subgroups = cellfun(addBinaryResp, subgroups,filters, 'UniformOutput',false);

subgroups = cellfun(@(x) removevars(x, predictor), subgroups,'UniformOutput',false);
subgroups = cellfun(@(x) renamevars(x, 'status',predictor), subgroups,'UniformOutput',false);

% Define regression formula
regFormula = string(strcat(response,'~',strjoin(confounders,'+')));

% preparedInputTableTmp = preparedInputTable;
% sampsToChange = matches(preparedInputTableTmp.ID,{'X42510869';'X42754461';'X42547915';'X4262232'});
% preparedInputTableTmp{sampsToChange,'DM_dchac[bc]'} = nan;

% Perform logistic regressions for each subgroup
[results_APOE,regressions] = cellfun(@(x) performRegressions(preparedInputTable,x,regFormula), subgroups,'UniformOutput',false);

% Find flux results and label results with subgroup types
addCohortFun = @(x,y) addvars(x.(predictor),repmat( string( strjoin(y,' – ') ) ,height(x.('age_at_collection')), 1), 'Before','N','NewVariableNames','subgroups');
results_APOE = cellfun(addCohortFun,results_APOE,filters,'UniformOutput',false);
% Concatinate results into a single table
results_APOE = vertcat(results_APOE{:});

% Ensure the APOE allele variable is categorical

% Perform logistic regressions for APOE E2 VS E4 and E3 VS E4
% [results_APOE,regressions] = performRegressions(preparedInputTable,preparedMetadata,regFormula);
% 
% % Combine the two APOE comparisons
% fields = string(fieldnames(results_APOE));
% results_APOE = vertcat( results_APOE.(fields(end-1)), results_APOE.(fields(end)) );
% results_APOE(~matches(results_APOE.Reaction,'DM_for[bc]'),:)=[];


% Reorder categories and regress again
% preparedMetadata.APOE_ALLELE = categorical(preparedMetadata.APOE_ALLELE,{'E4','E3','E2'});
% % Perform logistic regressions for APOE E2 VS E4 and E3 VS E4
% [results_APOE_reo,regressions_reo] = performRegressions(preparedInputTable,preparedMetadata,regFormula);
% 
% % Combine the two APOE comparisons
% fields = string(fieldnames(results_APOE_reo));
% results_APOE_reo = vertcat( results_APOE_reo.(fields(end-1)), results_APOE_reo.(fields(end)) );
% results_APOE_reo(~matches(results_APOE_reo.Reaction,'DM_for[bc]'),:)=[];
% 
% tst = [results_APOE;results_APOE_reo];

% Find groups for each reaction
groups = findgroups(results_APOE.Reaction);
% Test for each reaction if the results are consistent

groupConsistent = splitapply(@(x) all(x > 0) || all(x < 0), results_APOE.estimate, groups);
% Add consistency labels to regression results
results_APOE.consistent = groupConsistent(groups);

% Sort results
results_APOE = sortrows(results_APOE,'pValue','ascend');

% Filter on metabolites with at least on significant result
results_APOE = groupfilter(results_APOE,"Reaction",@(x) any(x<pFilter),"pValue");

end