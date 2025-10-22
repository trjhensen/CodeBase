function [results, regressions] = performGlobalCognitionAnalysis(preparedInputTable, preparedMetadata, metabolomicsData)

% inputs
% preparedInputTable
% preparedMetadata

if nargin < 3
    metabolomicsData = false;
end

% Associate the predicted fluxes with Cognitive scores

% State covariates to control for
confounders = {...
    'Sex',... % Male/Female    
    'age_at_collection',... % Age in years
    'EDUC',... % Education in years (UNIQUE for this analysis!)
    'mapped_species_reads',...
    'lane',... % Sequencing run lane
    'Ethanol_added',...
    ...'HYPERT',... % Hypertension present
    'DAILY_ALCOHOL',...% Daily alcohol consumption
    'APOE_E4'...
    };

if metabolomicsData == true
    confounders = {...
    'Sex',... % Male/Female    
    'age_at_collection',... % Age in years
    'EDUC',... % Education in years (UNIQUE for this analysis!)
    ...'HYPERT',... % Hypertension present
    'DAILY_ALCOHOL',...% Daily alcohol consumption
    'APOE_E4'...
    };
end

% State response variable of interest
response = 'G';

% Define regression formula
regFormula = string(strcat(response,'~','Flux+',strjoin(confounders,'+')));

% Remove samples without cognitive score information
preparedMetadata(isnan(preparedMetadata.G),:) = [];

% Perform regressions on the fluxes
[outputs,regressions] = performRegressions(preparedInputTable,preparedMetadata,regFormula);
% isoleucine and deoxycholate correlate with cognition scores. Not arginine
% and creatine

% Lets investigate the effect of AD diagnosis and APOE E4 status on the
% results

% First define parameters for interaction analysis
regressionParam = struct;
% Propagate standard parameters
regressionParam.response = response;
regressionParam.confounders = confounders;

% Define parameters specific to moderation analysis
regressionParam.interactionPvalThreshold = 0.05; % Only investigate if p<0.1

% Perform moderation analysis on AD status
regressionParam.modVar = 'AD';
stratResAPOE_AD = moderationAnalysisADRC(preparedInputTable, preparedMetadata, regressionParam);


% Perform moderation analysis on APOE E4 status
regressionParam.modVar = 'APOE_E4';
stratResAPOE_E4 = moderationAnalysisADRC(preparedInputTable, preparedMetadata, regressionParam);

% Generate results cell array
results = {outputs.Flux, stratResAPOE_AD, stratResAPOE_E4};
end