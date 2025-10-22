function [resCD,resAD,resCS,resE4] = compareFluxCorr(preparedInputTable, preparedMetadata)
%% ADRC flux analysis
% RQ: What is the role of gut microbiota in metabolic alterations of
% cognitive decline?

% Associate the predicted fluxes with changes in fluxes from normal cognition to mild cognitively impaired to full dementia. 
[DM_flux_regressions, ~] = performCogDeclineAnalysis(preparedInputTable, preparedMetadata);

% Prune on normal-dementia associations
resCD = DM_flux_regressions{1};
resCD = resCD(matches(resCD.subgroups,"normal cognition-Dementia"),:);

% Associate the predicted fluxes with AD status
[AD_flux_regressions, ~] = performADanalysis(preparedInputTable, preparedMetadata);
resAD = AD_flux_regressions{1,1};
resAD = removevars(resAD,'Cohort');

% Associate the predicted fluxes with Global cognition
[Cog_flux_regressions, ~] = performGlobalCognitionAnalysis(preparedInputTable, preparedMetadata);
resCS = Cog_flux_regressions{1,1};

% Perform flux regressions on APOE E4 genotype

% State response variable of interest
response = 'Flux';

% State covariates to control for
confounders = {...
    'Sex',... % Male/Female    
    'age_at_collection',... % Age in years
    'mapped_species_reads',...
    'lane',... % Sequencing run lane
    'Ethanol_added',...
    'DAILY_ALCOHOL'...
    };

% Define regression formula
regFormula = strcat(response,'~','APOE_E4+',strjoin(confounders,'+'));

% Perform regressions
outputs = performRegressions(preparedInputTable,preparedMetadata,regFormula);
resE4 = outputs.APOE_E4_NO_E4;

end

% end

%
% Which metabolites consistently changed?
% results = DM_flux_regressions{1};
% results(results.consistent==0,:)=[]; % Remove metabolites with inconsistent results
% 
% % Remove metabolites with no associations under a defined threshold
% pThreshold = 0.05; % Set threshold
% [groups,names] = findgroups(results.Reaction); % Find groups for each reaction
% rxnsToKeep = arrayfun(@(x) any(results.pValue(groups==x)<pThreshold), 1:max(groups)); % Check if one or more reactions are above the threshold
% results(~matches(results.Reaction,names(rxnsToKeep)),:) = []; % Prune list of reactions of interest
% 
% % 11 metabolites are associated with cognitive decline
% unique(results.Reaction);
% 
% %%
% 
% 
% adRes = AD_flux_regressions{1,1}; % Get ad-flux associations
% adResAlsoDem = adRes(matches(adRes.Reaction,unique(results.Reaction)),:); 
% adResNotDem = adRes(~matches(adRes.Reaction,unique(results.Reaction)),:);

% All metabolites that associated with AD diagnosis also associated with 
% cognitive decline.

% Only four of the eleven metabolites that associated with cognitive
% decline also associated with AD diagnosis. 

% 