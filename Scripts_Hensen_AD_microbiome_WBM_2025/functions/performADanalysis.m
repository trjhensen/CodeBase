function [results, regressions] = performADanalysis(preparedInputTable, preparedMetadata, metabolomicsData)

% inputs
% preparedInputTable
% preparedMetadata

if nargin < 3
    metabolomicsData = false;
end

% State response variable of interest
response = 'AD';

% State covariates to control for
confounders = {...
    'Sex',... % Male/Female    
    'age_at_collection',... % Age in years
    'mapped_species_reads',...
    'lane',... % Sequencing run lane
    'Ethanol_added',...
    'HYPERT',... % Hypertension present
    'DAILY_ALCOHOL'...% Daily alcohol consumption
    };

if metabolomicsData == true
    confounders = {...
    'Sex',... % Male/Female    
    'age_at_collection',... % Age in years
    'EDUC',... % Education in years (UNIQUE for this analysis!)
    'HYPERT',... % Hypertension present
    'DAILY_ALCOHOL',...% Daily alcohol consumption
    ...'APOE_E4'...
    };
end


% Define regression formula
regFormula = strcat(response,'~','Flux+',strjoin(confounders,'+'));

regressionType = 'logistic';
responseOrdination = {'no_AD','AD'};

% Preallocate variables to store the regression results
results = cell(3,2);
regressions = cell(3,2);

% Annotate result structures
results(:,2) = {'Flux-AD associations'; 'Added Flux*APOE_E4 interaction term';'Added Flux*NPS interaction term'};
regressions(:,2) = {'Flux-AD associations'; 'Added Flux*APOE_E4 interaction term';'Added Flux*NPS interaction term'};


% Perform regressions
[outputs,regressions{1,1}] = performRegressionsADRC(preparedInputTable,preparedMetadata,regFormula,regressionType,responseOrdination);

% Add cohort information so that all results in the future can be
% concatinated into a single table
results{1,1} = addvars(outputs.Flux,repmat("Full",height(outputs.Flux),1), 'Before','N','NewVariableNames','Cohort');

% Also investigate the effect of APOE E4 status on AD associations

% First define parameters for interaction analysis
regressionParam = struct;
% Propagate standard parameters
regressionParam.response = response;
regressionParam.confounders = confounders;
regressionParam.regressionType = regressionType;
regressionParam.responseOrdination = responseOrdination;

% Define parameters specific to moderation analysis
regressionParam.modVar = 'APOE_E4';
regressionParam.interactionPvalThreshold = 0.05; % Only investigate if p<0.1

% Perform moderation analysis on APOE E4 status
stratResAPOE_E4 = moderationAnalysisADRC(preparedInputTable, preparedMetadata, regressionParam);
results{2,1} = stratResAPOE_E4;
% Interpretation of results: Arginine and creatine are increased in AD patients.
% However, only if an individual is APOE E4 positive, increased
% deoxycholate and decreased ursocholate fluxes are found. 

% Are there also relationships with NPS?
regressionParam.modVar = 'NPS';
stratResNPS = moderationAnalysisADRC(preparedInputTable, preparedMetadata, regressionParam);
results{3,1} = stratResNPS;

% AD patients with neuropsychiatric symptoms have much lower ursocholate
% fluxes than those without neuropsychiatric symptoms.


end

