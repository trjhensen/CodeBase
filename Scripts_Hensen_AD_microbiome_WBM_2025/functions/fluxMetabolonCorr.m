function [corrTable, pvalTable] = fluxMetabolonCorr(fluxPath,metabolonPath,metadataPath, saveDir)

% Load the input flux data
[preparedFluxes, preparedMetadata] = prepareDataForStatsADRC(fluxPath, metadataPath,true);
preparedFluxes.Properties.VariableNames = erase(preparedFluxes.Properties.VariableNames,{'DM_','[bc]'});

% Load the metabolomic data
preparedPlasma = prepareDataForStatsADRC(metabolonPath, metadataPath, false);
preparedPlasma = removevars(preparedPlasma,'DM_dchac[bc]');

% Filter on all metabolites except the metabolites of interest
preparedPlasma = preparedPlasma(:,{'ID','arg_L','creat','taur'});

% Filter the flux reactions on the mapped plasma metabolites
preparedFluxes = preparedFluxes(:,preparedPlasma.Properties.VariableNames);

% Differentiate plasma names
preparedPlasma.Properties.VariableNames(2:end) = append(preparedPlasma.Properties.VariableNames(2:end),'_plasma'); 

% Append plasma to fluxes
%preparedData = outerjoin(preparedFluxes,preparedPlasma,'Keys','ID','MergeKeys',true);
% Add metadata
%preparedData = outerjoin(preparedData,preparedMetadata,'Keys','ID','MergeKeys',true,'Type','left');
preparedMetadataPlasma = outerjoin(preparedPlasma,preparedMetadata,'Keys','ID','MergeKeys',true,'Type','left');

% Generate formulae
confounders = {'Sex','age_at_collection','mapped_species_reads','lane'};
% confounders = {'Sex','mapped_species_reads','lane','Ethanol_added','APOE_E4'};
formulae = {...
    strcat('arg_L_plasma','~','Flux','+',strjoin(confounders,'+')),...
    strcat('creat_plasma','~','Flux','+',strjoin(confounders,'+')),...
    strcat('taur_plasma','~','Flux','+',strjoin(confounders,'+'))...
    };

%[results,regressions] = performRegressions(preparedFluxes,preparedMetadataPlasma,formulae{1});
[results,~] = cellfun(@(x) performRegressions(preparedFluxes,preparedMetadataPlasma,x), formulae,'UniformOutput',false);

% Find flux results and label results with subgroup types
results = cellfun(@(x) x.Flux,results,'UniformOutput',false);
results = cellfun(@(x) convertvars(x, 'Formula','string'),results,'UniformOutput',false);
% Concatinate results into a single table
results = vertcat(results{:});

% Filter flux data
% varsToKeep = [1 find(contains(preparedFluxes.Properties.VariableNames,preparedPlasma.Properties.VariableNames(2:end)))];
% preparedFluxes = preparedFluxes(:,varsToKeep);

% Correlate flux values with metabolomics
[RHO, PVAL] = corr(preparedFluxes{:,2:end}, preparedPlasma{:,2:end}, 'type', 'Pearson','Rows','pairwise');
close all; figure; imagesc(RHO); colorbar;
% Create correlation table
corrTable = array2table(RHO',"RowNames",preparedPlasma.Properties.VariableNames(2:end),'VariableNames',preparedFluxes.Properties.VariableNames(2:end));
pvalTable = array2table(PVAL',"RowNames",preparedPlasma.Properties.VariableNames(2:end),'VariableNames',preparedFluxes.Properties.VariableNames(2:end));

% save results
fileName = fullfile(saveDir,'flux_metabolon_corr.xlsx');

writetable(results,fileName,'Sheet','Regressions','WriteRowNames',true);
writetable(corrTable,fileName,'Sheet','Spearman_rho','WriteRowNames',true);
writetable(pvalTable,fileName,'Sheet','Spearman_pval','WriteRowNames',true);
end
