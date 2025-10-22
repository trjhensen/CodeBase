function [results, regressions] = performCogDeclineAnalysis(preparedInputTable, preparedMetadata, metabolomicsData)

% inputs
% preparedInputTable
% preparedMetadata

if nargin < 3
    metabolomicsData = false;
end

% Associate the predicted fluxes with Dementia and MCI status

% State covariates to control for
confounders = {...
    'Sex',... % Male/Female    
    'age_at_collection',... % Age in years
    'mapped_species_reads',...
    'lane',... % Sequencing run lane
    'Ethanol_added',...
    'HYPERT',... % Hypertension present
    'DAILY_ALCOHOL',...% Daily alcohol consumption
    ...'APOE_E4'...
    };

if metabolomicsData == true
    confounders = {...
    'Sex',... % Male/Female    
    'age_at_collection',... % Age in years
    'HYPERT',... % Hypertension present
    'DAILY_ALCOHOL',...% Daily alcohol consumption
    ...'APOE_E4'...
    };
end

% State response variable of interest
response = 'NACCUDSD';

% Define regression formula
regFormula = string(strcat(response,'~','Flux+',strjoin(confounders,'+')));

% Remove samples without cognitive status information
preparedMetadata(matches(preparedMetadata.NACCUDSD,{'<undefined>','impaired non-mci'}),:) = [];

% Generate subgroups
filters = { {'normal cognition','MCI'}, {'MCI','Dementia'}, {'normal cognition','Dementia'} };
subgroups = cellfun(@(x) preparedMetadata(matches(preparedMetadata.NACCUDSD, x),:), filters,'UniformOutput',false);

% For each subgroup, encode the less severe group as 0 and the more severe
% group as 1
convertBinaryResp = @(x,y) grp2idx(categorical(x.(response),y))-1;
addBinaryResp = @(x,y) addvars(x,convertBinaryResp(x,y), 'After',response,'NewVariableNames','status');
subgroups = cellfun(addBinaryResp, subgroups,filters, 'UniformOutput',false);


% Perform logistic regressions for each subgroup
[results,regressions] = cellfun(@(x) performRegressions(preparedInputTable,x,regFormula), subgroups,'UniformOutput',false);

% Find flux results and label results with subgroup types
addCohortFun = @(x,y) addvars(x.Flux,repmat( string( strjoin(y,'-') ) ,height(x.Flux), 1), 'Before','N','NewVariableNames','subgroups');
results = cellfun(addCohortFun,results,filters,'UniformOutput',false);
% Concatinate results into a single table
results = vertcat(results{:});

% Find reactions with consistent changes: normal -> MCI -> Dementia

% Find groups for each reaction
groups = findgroups(results.Reaction);

% Test for each reaction if the results are consistent
groupConsistent = splitapply(@(x) all(x > 0) || all(x < 0), results.estimate, groups);
% Add consistency labels to regression results
results.consistent = groupConsistent(groups);


if 1
% Perform moderation analysis on healthy vs dementia for AD 

% First define parameters for interaction analysis
regressionParam = struct;
% Propagate standard parameters
regressionParam.response = response;
regressionParam.confounders = confounders;

% Define parameters specific to moderation analysis
% regressionParam.modVar = 'AD';
regressionParam.interactionPvalThreshold = 1; % Only investigate if p<0.1
% 
% % Perform moderation analysis on AD status (MCI vs Dementia)
% stratResAD = moderationAnalysisADRC(preparedInputTable, subgroups{2}, regressionParam);
% No differences were found

% Perform moderation analysis on APOE status
regressionParam.modVar = 'APOE_E4';

% Perform moderation analysis on each dementia group
stratResAPOE_E4 = cellfun(@(x) moderationAnalysisADRC(preparedInputTable, x, regressionParam), subgroups,'UniformOutput',false);

% Remove empty cells
cellToRm = cellfun(@isempty,stratResAPOE_E4);
stratResAPOE_E4(cellToRm) = [];
filters(cellToRm) = [];

% Find flux results and label results with subgroup types
addCohortModFun = @(x,y) addvars(x,repmat( string( strjoin(y,'-') ) ,height(x), 1), 'Before','N','NewVariableNames','subgroups');
stratResAPOE_E4 = cellfun(addCohortModFun,stratResAPOE_E4,filters,'UniformOutput',false);
% Concatinate results into a single table
stratResAPOE_E4 = vertcat(stratResAPOE_E4{:});

% results = {results, stratResAPOE_E4, stratResAD};
results = {results, stratResAPOE_E4};
end
% 
% % Define regression types for ordinal multinomial regression
% regressionType = 'multinomial';
% responseOrdination = {'normal cognition','impaired non-mci','MCI','Dementia'};
% 
% % Preallocate variables to store the regression results
% results = cell(3,2);
% regressions = cell(3,2);
% 
% % Annotate result structures
% results(:,2) = {'Flux-Dementia ordinal associations'; 'Added Flux*AD interaction term';'Added Flux*APOE_E4 interaction term'};
% regressions(:,2) = {'Flux-Dementia ordinal associations'; 'Added Flux*AD interaction term';'Added Flux*APOE_E4 interaction term'};
% disp(regFormula)
% [results{1,1},regressions{1,1}] = performRegressionsADRC(preparedInputTable,preparedMetadata,regFormula,regressionType,responseOrdination);
% 
% % Investigate the effect of AD status on the flux-dementia associations
% regFormula = strcat(response,'~','Flux+AD+',strjoin(confounders,'+'),'+Flux:AD'); 
% disp(regFormula)
% [results{2,1},regressions{2,1}] = performRegressionsADRC(preparedInputTable,preparedMetadata,regFormula,regressionType,responseOrdination);
% 
% % Interpretation of results:
% % If a person develops dementia without AD, bile acid production are unchanged or
% % slightly decreased. However, if that person has AD, bile acid productions
% % are increased. 
% 
% 
% % Investigate the effect of APOE E4 status on the flux-dementia
% % associations
% regFormula = strcat(response,'~','Flux+APOE_E4+',strjoin(confounders,'+'),'+Flux:APOE_E4'); 
% disp(regFormula)
% [results{3,1},regressions{3,1}] = performRegressionsADRC(preparedInputTable,preparedMetadata,regFormula,regressionType,responseOrdination);

% Interpretation: Dementia patients with E4 have lower ursocholate fluxes
% (good bile acid).

end
