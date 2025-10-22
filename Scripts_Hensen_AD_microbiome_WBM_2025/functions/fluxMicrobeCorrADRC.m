function [RHO_table,fluxMicrobeCorrPath,processedMicrobesPath] = fluxMicrobeCorrADRC(fluxPath, mappedMicrobePath, metadataPath, fmLinksPath, saveFolder)
% Function for assessing Spearman correlations between the predicted fluxes and
% relative microbe abundances. Correlations are only performed for
% microbe-reaction combinations where the species biomass metabolite has a
% nonzero (tol=1e-6) shadow price for the optimised reaction. Microbes are also only
% correlated if they are present in at least a certain percentage of
% samples, defined by the microbeCutoff variable. 


% Load metadata
samplesToFilter = readcell(metadataPath,'Range','A2:A9999');

% Define functions for flux and microbiome data wrangling
repairIdFun = @(x) erase(x.Properties.RowNames,{'mWBM_','_female','_male'});
filterSamps = @(x) x(matches(x.Properties.RowNames,samplesToFilter),:);

% Load flux data
fluxes = readtable(fluxPath,'VariableNamingRule','preserve','ReadRowNames',true);
fluxes.Properties.RowNames = repairIdFun(fluxes);
fluxes = filterSamps(fluxes);
fluxes = removevars(fluxes,'Sex');

% Load gut microbiome data
microbiome = readtable(mappedMicrobePath,'VariableNamingRule','preserve','ReadRowNames',true);
microbiome.Properties.RowNames = repairIdFun(microbiome);
microbiome = filterSamps(microbiome);
if any(matches(microbiome.Properties.VariableNames,'Sum of taxa'))
microbiome = removevars(microbiome,'Sum of taxa');
end

% Make sure that the flux and microbiome samples are identical and identically ordered
fluxes = fluxes(microbiome.Properties.RowNames, :);
microbiome = microbiome(fluxes.Properties.RowNames, :);

% Filter the flux data on metabolites of interest
fmAssociations = readtable(fmLinksPath,'VariableNamingRule','preserve','ReadRowNames',true);
fluxes = fluxes(:,fmAssociations.Properties.VariableNames);

% Remove microbial species that were not in the models in less than 10% of
% the samples
thresholdVal = 0.9;
microbesToRemove = sum(isnan(table2array(microbiome))) > height(microbiome)*thresholdVal;
microbiome(:,microbesToRemove) = [];

% Perform spearman correlations
fluxArray = table2array(fluxes);
microbiomeArray = table2array(microbiome);
RHO = corr(microbiomeArray,fluxArray,'type','Spearman','rows','pairwise');

% Get variable names
fluxVars = fluxes.Properties.VariableNames;
microbiomeVars = microbiome.Properties.VariableNames;

% Set all incorrect microbe-flux correlations to nan
fluxMicrobeLinks = fmAssociations{microbiomeVars, fluxVars}; % Order the known flux-microbe links according to RHO
RHO(~fluxMicrobeLinks) = nan;
RHO_table = array2table(RHO,'VariableNames',fluxVars,'RowNames',microbiomeVars);

if 0
% Conver VMH IDs to Metabolites names
RHO_table.Properties.VariableNames = renameAdrcVmhToMetName(RHO_table.Properties.VariableNames);
end
% Save results
fluxMicrobeCorrPath = fullfile(saveFolder,'flux_microbe_corr.csv');
writetable(RHO_table,fluxMicrobeCorrPath,'WriteRowNames',true);

processedMicrobesPath = fullfile(saveFolder,'processedWBM_relative_abundances.csv');
writetable(microbiome,processedMicrobesPath,'WriteRowNames',true);

end