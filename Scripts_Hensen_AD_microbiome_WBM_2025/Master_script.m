%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Master script %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% AUTHOR: Tim Hensen,
% INSTITUTE: University of Galway, School of medicine
%
% DESCRIPTION:
% This script contains the complete computational analysis described in 
% "Personalised whole-body modelling links gut microbiota to metabolic 
% perturbations in Alzheimer's disease". Before running this code, please
% put the inputs in a folder "input" under the project directory.

clear;clc;
restoredefaultpath
addpath(genpath('/home/tim/Documents/CodeBase'))
addpath(genpath('/home/tim/Documents/ADRC/'))
addpath(genpath('/home/tim/Documents/cobratoolbox'))
addpath(genpath('/home/tim/Documents/wbm_modelingcode'))

% ADRC_stats_sketch.m

% Set paths for analysis
paths = struct;
paths.root = what('ADRC').path; % Set working directory
paths.inputs = fullfile(paths.root,'inputs');
paths.outputs = fullfile(what('ADRC').path,'outputs'); % Set directory to save analysis results

% Set paths for new flux outputs
paths.fluxes = fullfile(paths.outputs,'fluxes');
paths.FBA  = fullfile(paths.outputs,'fluxes','FBA');
paths.figures = fullfile(paths.outputs,'figures');
paths.fluxes = fullfile(paths.outputs,'fluxes');
paths.fluxAnalysis = fullfile(paths.fluxes,'analysis');
paths.fluxPath = fullfile(paths.fluxAnalysis,'processed_fluxes.csv');

paths.microbetoflux = fullfile(paths.outputs,'microbetoflux');

% Set paths for microbiome processing
paths.microbiome = fullfile(paths.inputs,'microbiome');
paths.unprocessedMicrobes = fullfile(paths.microbiome,'unprocessed');
paths.processedMicrobes = fullfile(paths.microbiome,'processed');

% Set microbiome processing outputs
paths.processedMicrobiome = fullfile(paths.outputs,'processedMicrobiome');
paths.microbiotaPath = fullfile(paths.outputs,'resultMARS','normalized_preMapped','normalized_preMapped_species.csv'); 
paths.microbiotaPathFilt = fullfile(paths.processedMicrobiome,'filtered_normalized_preMapped_species.csv');
paths.mappedMicrobePath = fullfile(paths.fluxAnalysis,'WBM_relative_abundances.csv'); % Microbial relative abundances
paths.microbiotaWbmPathFilt = fullfile(paths.processedMicrobiome,'filtered_WBM_relative_abundances.csv');
if ~isfolder(paths.processedMicrobiome); mkdir(paths.processedMicrobiome); end

% Set metadata paths
paths.metadataOutputFolder = fullfile(paths.outputs,'metadata');
paths.metadataRaw = fullfile(paths.inputs,'NACC_NCRAD_FECALSET_COMBINED_NO_NATIVE_AMER_2024-03-04_2144_COLETTE.csv');
paths.metadataProcessed = fullfile(paths.metadataOutputFolder,'prunedProcessedMetadata.csv');
paths.metadataPrunedProcessed = fullfile(paths.metadataOutputFolder,'fluxPrunedProcessedMetadata.csv');
if ~isfolder(paths.metadataOutputFolder); mkdir(paths.metadataOutputFolder); end % Create folder if not yet available

% Set summary statistics paths
paths.metadataSummaryStats = fullfile(paths.metadataOutputFolder,'metadataSummaryStats.xlsx');
paths.mappingSummaryStats = fullfile(paths.metadataOutputFolder,'microbiomeMappingSummaryStats.xlsx');

% Set outputs for the statistical results
paths.apoe = fullfile(paths.outputs,'APOE');
paths.cognition = fullfile(paths.outputs,'cognition');
paths.AD = fullfile(paths.outputs,'AD');

if 0
    % Find current folders and remove all outputs except for the modelling data
    foldersToRm = setdiff({dir(paths.outputs).name}, {'fluxes','microbiome','Knirps','ResultMARS','ResultMARS_NEW','resultMgPipe','Sneezy','logFile_initialisation.txt','.','..'});
    foldersToRm = fullfile(paths.outputs,foldersToRm); % Generate paths
    
    for i = 1:length(foldersToRm)
        if isfolder(foldersToRm{i})
            rmdir(foldersToRm{i}, 's') % Remove previous outputs
            mkdir(foldersToRm{i}) % Create new empty folders
        end
    end
    
    % Make sure that the core output folders exist
    cellfun(@mkdir, {paths.figures,paths.fluxes,paths.processedMicrobiome, paths.microbetoflux}); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Metagenomic read count mapping %%%%%%%%%$%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clean up workspace 
clearvars -except paths; clc;

% Convert the taxonomy mapped biom files to txt files
% This step is done in docker 
% biom convert -i 194939/free.biom -o 194939/free.txt --to-tsv
% biom convert -i 194890/free.biom -o 194890/free.txt --to-tsv
%
% Decontaminate read counts in R: 
% scriptPath = which('decontamination_script.R')

% Process microbiome data for MARS
paths.Mars.inputTable = fullfile(paths.processedMicrobes,'MARS_input_combined.csv');
if ~isfile(paths.Mars.inputTable)
    [combinded, paths.Mars.inputTable] = processMicrobiomesForMars(paths.unprocessedMicrobes, paths.processedMicrobes);
end
% Fix names in gut microbiota to ensure that all microbiome samples can be
% mapped onto the metadata

% Set paths for MARS
paths.Mars.outputPathMars = fullfile(paths.outputs, 'ResultMARS');

% Perform metagenomic mapping with MARS
if ~isfolder(paths.Mars.outputPathMars) % Only needs to be run once
    OK = runMarsADRC(paths.Mars.outputPathMars, paths.Mars.inputTable);

    % Investigate the effect of mapping currently unmapped taxa 
    outputPathMars = paths.Mars.outputPathMars;
    [mappingPotential, savePath] = findPotentialReadCovIncr(outputPathMars);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% METADATA INVESTIGATION AND PRUNING %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clean up workspace 
clearvars -except paths; clc;

% Set paths for ease of use
microbiomeInputFolder = paths.microbiome;
outputPathMars = paths.Mars.outputPathMars;
metadataOutputFolder = paths.metadataOutputFolder;

% Process the patient metadata
metadata = processMetadataADRC(paths.metadataRaw, paths.Mars.inputTable);

% Add cognitive score information
createFig = false;
metadataCog = appendGlobalCognitionMetadata(metadata, metadataOutputFolder, createFig);

% Add metadata from microbiome samples
mergedMetadata = appendMicrobiomeMetadata(metadataCog, microbiomeInputFolder);

% Add microbiome diversity metrics to metadata
marsFolder = paths.Mars.outputPathMars;
metadataDiversities = appendDiversityMetricsMetadata(mergedMetadata,paths.Mars.outputPathMars);

% Add MARS statistics to metadata (diversity metrics and read coverage before and after mapping)
metadataMars = appendMappingStatsMetadata(metadataDiversities, outputPathMars);

% Append technical covariates, such as vial type and ethanol usage. 
metadataTechCov = appendTechnicalCovariatesMetadata(metadataMars, microbiomeInputFolder);
% Technical covariate analysis: adrcTechnicalCovariateAnalysis

% Remove samples with multiple timepoints and samples without dementia
% information
metadataProcessed = rmAdrcOverlapAndMissing(metadataTechCov);
groupsummary(metadataProcessed,'NACCUDSD')
% Remove low quality samples (Check IBD status, alcohol intake, pain
% medication, smoking, mood disorders, Intestinal inflammation, sleep-aid
% medication. Remove samples if almost non of the
% individuals have a yes for these metadata.)
[metadataPruned, metadataIntermedPruning] = pruneMetadataADRC(metadataProcessed);

%%
% Number of samples with BMI information
tst = metadataPruned(~isnan(metadataPruned.NACCBMI),:);
groupsummary(metadataPruned(~isnan(metadataPruned.NACCBMI),:),'NACCUDSD')
%%

% Save updated metadata file with technical co-variate information
writetable(metadataPruned,paths.metadataProcessed)

% Describe the effects of metagenomic mapping the gut microbiome relative abundances

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODELLING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Next, I will generate host-microbiome WBMs using the Persephone
% pipeline. I adapted the Persephone config file for generation of
% microbiome community models and host-microbiome WBMs.

% Clean up workspace
% clearvars -except paths; clc;
% runPersephone(which('mgPipeAdrcConfig.m')); end % Run mgpipe in Persephone

% After having created the microbiome community models, I uploaded half of
% microbiome models to the Knirps server and the other half to the Sneezy
% server. There, I generated mWBMs and generated fluxes for 38 metabolites.
% I used the following script to achieve this:
% [modelStats, summaryStats, dietInfo, dietGrowthStats] = createMWBMsADRC(); % Creation of mWBMs
% [metList, vmhList] = loadADRCMets(); % Definition of metabolites of interest
% runPersephone(which('fbaAdrcConfig.m')); % Run Persephone to predict fluxes 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Processing of FBA solutions %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clean up workspace
clearvars -except paths; clc;

% Move FBA results from the two dwarfs to a shared folder 
if ~isfolder(paths.FBA); OK = moveAdrcFbaRes(paths.root,paths.FBA); end

if ~isfolder(paths.fluxAnalysis)
    % Create new folder if not present yet
    mkdir(paths.fluxAnalysis); 

    % Next, the parameters for flux processing are set. Metabolites for which a
    % solution could be found in 5% of samples or less are removed.
    paramFluxProcessing.rxnRemovalCutoff = {'fraction',  0.05};
    % paramFluxProcessing.fluxMicrobeCorrelationMetric = 'spearman_rho';

    % Process flux results
    analyseGF = true;
    analyseWBMsol(paths.FBA,paramFluxProcessing, paths.fluxAnalysis,analyseGF);

    % Calculate microbe to flux contributon potentials; 
    disp('Extract microbial abundances and shadow prices')
    extractMicrobeContributionsADRC(paths.FBA, paths.fluxAnalysis);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Flux outlier removal %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clean up workspace
clearvars -except paths; clc;

% Lets investigate potential sample outliers. 
close all; [sampleImportanceTable, explVar] = sampleOutliersInFluxesADRC(paths.fluxPath, paths.metadataProcessed, paths.fluxAnalysis);
%close all;
% Inspecting the outliers in the fluxes found that:
% 1) The top outlier (X42680612) had schizophrenia and AD
% 2) The second top (X42816255) outlier had no AD, but had a colon resection (partly removed) in 2016.

% Load metadata, remove samples and save updated metadata file
% [paths.metadata, metadata] = pruneFluxOutliersFromMetadataADRC(sampleImportanceTable,paths.metadata,2); % Remove the top 2 outliers in the fluxes
metadata = readtable(paths.metadataProcessed,'VariableNamingRule','preserve');
metadata(matches(metadata.ID,sampleImportanceTable.ID(1:2)),:) = []; % Remove samples
writetable(metadata,paths.metadataPrunedProcessed); % Save updated metadata file

% Generate summary statistics table
summaryStats = makeADRCmetadataTable(metadata);
writetable(summaryStats,paths.metadataSummaryStats,'WriteRowNames',false,'WriteMode','replacefile')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Metabolomics processing %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Map metabolon metabolomics data and append it to the metadata variables
paths.rawMetabolonPath = fullfile(paths.inputs,'metabolomics','ADRC Metabolon Preprocessed Unblinded 05102024.xlsx'); % Metabolon samples
rxnsToMap = readcell(paths.fluxPath,'Range','1C:1ZZZ');
[~,paths.metabolonPath] = appendMetabolonToMetada(rxnsToMap, paths.rawMetabolonPath, paths.metadataPrunedProcessed, paths.outputs);

% Now, we will test how well the fluxes of the selected metabolites can
% explain the metabolomic measurements.
% Correlate all flux predictions with the metabolomics abundances
[results, ~, ~] = fluxMetabolonCorr(paths.fluxPath,paths.metabolonPath,paths.metadataPrunedProcessed); % Flux-metabolome correlation analysis

% Obtain correlations with 95% confidence intervals
fluxPath = paths.fluxPath;
metabolonPath = paths.metabolonPath;
corrTable = fluxMetabolonCorr2(fluxPath,metabolonPath);

% save results
fileName = fullfile(paths.fluxes,'flux_metabolon_corr.xlsx');
writetable(results,fileName,'Sheet','Regressions');
writetable(corrTable,fileName,'Sheet','Spearman_rho');

% Also generate metadata summary file for the plasma metabolomics dataset
% Load table with plasma metabolomics
plasmaTable = readtable(paths.metabolonPath,'VariableNamingRule','preserve');
% Remove all nan samples
plasmaTable(isnan(plasmaTable{:,2}),:)=[];
% Filter metadata table
metadata = readtable(paths.metadataPrunedProcessed,'VariableNamingRule','preserve');
metadataPlasma = metadata(matches(metadata.ID,plasmaTable.ID),:);
% Create summary statistics table
summaryStatsPlasma = makeADRCmetadataTable(metadataPlasma);
% Write table to file
writetable(summaryStatsPlasma,fullfile(paths.metadataOutputFolder,'metadataPlasmaSummaryStats.xlsx'),'WriteRowNames',false,'WriteMode','replacefile')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Microbiome processing and summary statistics %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set paths for microbiome dataset pruning:
microbiotaWbmPath = paths.mappedMicrobePath;
microbiotaPath = paths.microbiotaPath;
metadataPath = paths.metadataPrunedProcessed;

% Prune microbiome tables for further statistical descriptions
[microbiome, wbmMicrobiome] = pruneMicrobiomeDataForStats(microbiotaWbmPath, microbiotaPath, metadataPath);

% Save updated tables
cellfun(@(x,y) writetable(x,y,'WriteRowNames',true),{microbiome, wbmMicrobiome},{paths.microbiotaPathFilt,paths.microbiotaWbmPathFilt})

% Investigate the effect of mapping on each class
% Set inputs
taxonomyPath = fullfile(paths.outputs,'resultMARS','preprocessedInput_afterRenaming.csv');
microbiotaPathFilt = paths.microbiotaPathFilt;
microbiotaWbmPathFilt = paths.microbiotaWbmPathFilt;
[mappingPotential, mappingPotentialFilt,savePath] = findPotentialReadCovIncr(paths.Mars.outputPathMars, microbiotaPathFilt);

% Perform an under-representation analysis
taxonEnrichmentTables = findTaxonomicOverrepInUnmappedTaxa(microbiotaPathFilt,taxonomyPath,microbiotaWbmPathFilt); 

% Save enrichment statistics to excel table
savePath = fullfile(paths.processedMicrobiome,'mappingLossEnrichment.xlsx');
cellfun(@(x,y) writetable(x, savePath, 'Sheet',y), taxonEnrichmentTables, cellfun(@(x) x.Properties.VariableNames(1), taxonEnrichmentTables))

% Next, generate summary statistics before and after mapping:
preMappedPath = fullfile(paths.outputs,'resultMARS', 'normalized_preMapped','normalized_preMapped_species.csv');
diversityStats = mappingSummaryStats(preMappedPath, microbiotaWbmPathFilt, metadataPath);
writetable(diversityStats,paths.mappingSummaryStats,'WriteRowNames',true,'WriteMode','replacefile')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Differential flux analysis %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clean up workspace
clearvars -except paths; clc; close all;

% Create new empty folders for results on APOE, cognition, and AD progression
newFolders = {paths.apoe, paths.cognition, paths.AD};

% Create empty folder for results
for i = 1:length(newFolders)
    if ~isfolder(newFolders{i})
        mkdir(newFolders{i}); 
    else
        rmdir(newFolders{i}, 's')
        mkdir(newFolders{i})
    end
end

% Dementia differential flux analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Now, we will investigate associations between the predicted and cognitive
% decline. Is there a difference between AD and controls? Are these results 
% unique to AD status or could they be explained by dementia status?
% Is there a difference between dementia and controls?


% Load and prepare input data and metadata
[preparedInputTable, preparedMetadata] = prepareDataForStatsADRC(paths.fluxPath, paths.metadataPrunedProcessed,false);
%preparedMetadata.(param.response) = renamecats( preparedMetadata.(param.response) ,dementiaCats,{'Healthy controls','MCI','AD dementia'});

% Number of samples with BMI information
tst = preparedMetadata(~isnan(preparedMetadata.NACCBMI),:);
groupsummary(preparedMetadata(~isnan(preparedMetadata.NACCBMI),:),'NACCUDSD')

tst = preparedMetadata(~isnan(preparedMetadata.NACCBMI),:);
groupsummary(preparedMetadata(~isnan(preparedMetadata.NACCBMI),:),'NACCUDSD')

% State confounders
paths.adConfounders = {'Sex','age_at_collection','mapped_species_reads','lane','Ethanol_added'};

if 1
    % dchacSampsToRm = {'X42510869','X42754461','X42547915','X4262232'};
%    dchacSampsToRm =  {'X4273975'
% 'X42409516'
% 'X42562933'
% 'X42808065'
% 'X42487528'
% 'X42373971'};
    %preparedInputTable(matches(preparedInputTable.ID,dchacSampsToRm),:)=[];
    paths.adConfounders = {'Sex','age_at_collection','NACCBMI','EDUC','APOE_E4','lane','mapped_species_reads','Ethanol_added'}; 

    % % Define regression formula
    % confounders = {'Sex','age_at_collection','mapped_species_reads','lane','Ethanol_added'};
    % form = strcat("AD ~ Flux + ", strjoin(confounders,'+') );
    % 
    % newCov = {'NACCBMI','EDUC','APOE_E4'};
    % form1 = append(form, '+', newCov, ' + Flux:',newCov);
    % form1 = append(form, '+', newCov);
    % 
    % confounderInfluence = cellfun(@(x) performRegressions(preparedInputTable, preparedMetadata,x), form1,'UniformOutput',false);

end

% Perform analysis
[results_AD, paths.adRxnsOfInterest, results_DM, regressions_DM] = alzheimerAnalysis(preparedInputTable,preparedMetadata,'Flux',paths.adConfounders);

% Save identified metabolites
paths.dementiaRxnsOfInterest = unique(results_DM.Reaction(results_DM.pValue<0.05),'stable');

% Create regression table for boxplot figure
resTableDM = buildRegressionTables(paths.dementiaRxnsOfInterest, results_DM, regressions_DM, 'AD');

% Create figure for flux results

% Set function parameters
param.response = "NACCUDSD"; param.titleAnnotation = " blood fluxes"; param.yTitle = "Flux in mmol/day/person"; param.edgeAlpha = 0.2; param.addEmptyTile=true;

% Prepare metadata for visualisation
dementiaCats = {'Healthy controls','MCI','AD dementia'};
preparedMetadata.(param.response) = categorical(preparedMetadata.(param.response),dementiaCats);
% preparedMetadata.(param.response) = renamecats( preparedMetadata.(param.response) ,dementiaCats,{'Healthy controls','MCI','AD dementia'});

% Create tiled figure
% results_DM(:,6:end-1) = fillmissing(results_DM(:,6:end-1),'constant',1);
fig = figure('Position',[39,457/3,1878,421*2]);
tiledlayout(2,length(paths.dementiaRxnsOfInterest),'TileSpacing','tight','Padding','loose');
createMultipleBoxPlotsADRC(preparedInputTable, preparedMetadata, paths.dementiaRxnsOfInterest, results_DM, param);

% Generate paths to save results
adRegTabPath = fullfile(paths.AD,'AD_progression_results.xlsx');
adRegTabSumPath = fullfile(paths.AD,'AD_RegressionTableSummary.xlsx');
adBoxplotPath = fullfile(paths.figures,'AD_Flux_boxplots.png');

% Write results to file
writetable(results_DM,adRegTabPath,'Sheet','Fluxes')
writetable(resTableDM,adRegTabSumPath,'Sheet','Fluxes','Range','A2')
exportgraphics(fig,adBoxplotPath,'Resolution',300)

% APOE differential flux analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Previously, we also found an association between APOE allele status and
% bile acid fluxes. Here, we will try to replicate this analysis 

% Load flux data
[preparedInputTable, preparedMetadata] = prepareDataForStatsADRC(paths.fluxPath, paths.metadataPrunedProcessed, false); 

% Filter on reactions that differed in the the AD patients
preparedInputTable = preparedInputTable(:,["ID", paths.dementiaRxnsOfInterest']);


% Only investigate cognitively normal samples
preparedMetadata = preparedMetadata(matches(preparedMetadata.NACCUDSD,'Healthy controls'),:);

% State response variable of interest for cognitive decline regressions
response = 'Flux';

% State covariates to control for
paths.apoeFluxConfounders = {'Sex','age_at_collection','NACCBMI','mapped_species_reads','lane','Ethanol_added','APOE_ALLELE'};

% Perform regression analysis to find associations between APOE allele
% status and the flux predictions
pFilter = 1;
[results_APOE,regressions_APOE] = apoeAnalysis(preparedInputTable,preparedMetadata, paths.apoeFluxConfounders, pFilter);

% Save identified metabolites
paths.apoeRxnsOfInterest = paths.dementiaRxnsOfInterest;%unique(results_APOE.Reaction(results_APOE.pValue<pFilter),'stable');
% paths.apoeRxnsOfInterest = unique(results_APOE.Reaction(results_APOE.pValue<0.05),'stable');

% Prepare apoe group data
preparedMetadata.APOE_ALLELE = categorical(preparedMetadata.APOE_ALLELE,{'ϵ2','ϵ3','ϵ4'});

% Create regression table for figure 1:
resTable = buildRegressionTables(paths.apoeRxnsOfInterest, results_APOE, regressions_APOE, 'APOE');

% Create figure of flux results

% Set function parameters
param.response = "APOE_ALLELE"; param.titleAnnotation = " blood fluxes"; param.yTitle = "Flux in mmol/day/person"; param.addEmptyTile=true;
param.edgeAlpha = 0.2;
param.apoeFlag = 1;


% Create tiled figure
fig = figure('Position',[39,457/3,1878,421*2]); %39,131,1878,747
tiledlayout(2,length(paths.apoeRxnsOfInterest),'TileSpacing','tight','Padding','loose');
createMultipleBoxPlotsADRC(preparedInputTable, preparedMetadata, paths.apoeRxnsOfInterest, results_APOE, param);

% Generate paths to save results
apoeRegTabPath = fullfile(paths.apoe,'APOE_regressions.xlsx');
apoeRegTabSumPath = fullfile(paths.apoe,'apoeRegressionTableSummary.xlsx');
apoeBoxplotPath = fullfile(paths.figures,'APOE_flux_boxplots.png');

% Write results to file
writetable(results_APOE,apoeRegTabPath,'Sheet','Fluxes')
writetable(resTable,apoeRegTabSumPath,'Sheet','Fluxes')
exportgraphics(fig,apoeBoxplotPath,'Resolution',300)

%%
% Cognition differential flux analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Previously, bile acids have been associated with global cognitive
% functioning and AD dementia progression. However, we did not find an 
% association with bile acids and AD dementia status. Thus, we will test if 
% we replicate the  previous findings by associating the fluxes to global cognition
% in normal, MCI, and AD-dementia patients?

% Load and prepare input flux data and metadata
[preparedInputTable, preparedMetadata] = prepareDataForStatsADRC(paths.fluxPath, paths.metadataPrunedProcessed);

% Filter on reactions that differed in the the AD patients
preparedInputTable = preparedInputTable(:,["ID", paths.dementiaRxnsOfInterest']);

if 0 
% Only investigate cognitively normal samples
preparedMetadata = preparedMetadata(matches(preparedMetadata.NACCUDSD,'Healthy controls'),:);
end

% State predictor of interest and covariates to control for
paths.cognitionConfounders = {'Sex','age_at_collection','EDUC','NACCBMI','mapped_species_reads','lane','Ethanol_added','APOE_E4'};

% paths.cognitionConfounders = {'Sex','age_at_collection','EDUC','NACCBMI','APOE_E4'};

% Prepare regression formula
predictor = 'Flux';
confounders = paths.cognitionConfounders;
response = 'NACCMOCA'; %preparedMetadata.NACCMOCA = normalize(log(preparedMetadata.NACCMOCA));
% response = 'G';
regFormula = string(strcat(response,'~',predictor,'+',strjoin(confounders,'+')));

% Perform regressions on global cognitive scores
[results_G,regressions_G] = performRegressions(preparedInputTable,preparedMetadata,regFormula);
results_G = addvars(results_G.Flux, repmat("normal cognition", height(results_G.Flux),1),'NewVariableNames','subgroups','After','Regression type');

% pFilter = 1;
% results_G = cognitiveScoreAnalysis(preparedInputTable,preparedMetadata, 'Flux', paths.cognitionConfounders, pFilter);

% Write results to file
filePath = fullfile(paths.cognition,'Cognition_regressions.xlsx');
writetable(results_G,filePath,'Sheet','Fluxes')



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Flux-microbe associations %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% After having analysed the flux results, I want to find the associated gut
% microbes. So, first I will select metabolites of interest from the
% previous flux investigations. Then, I will run the flux-microbe pipeline.
% Lastly, the identified microbial species will be analysed against apoe,
% cognition, and AD like the fluxes. 

% Clean up workspace
clearvars -except paths; clc; close all;

% Flux-microbe associations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find the reactions of interest
paths.rxnsOfInterest = paths.adRxnsOfInterest;%unique(vertcat(paths.adRxnsOfInterest,paths.apoeRxnsOfInterest,paths.cognitionRxnsOfInterest));

% Pre-define parameters for identifying the most important microbial species
% for the flux results
param = struct;
param.rxnsOfInterest = paths.rxnsOfInterest;
param.bootSamp = 1e4; % Boot samples for obtaining the mean and 95%CI of microbe contribution potentials
param.minFreq = 0.9;
param.nBootLasso = 500;
param.enBoot = 1e5; % Elastic net bootstrap samples 

paths.shadowPriceDir = fullfile(paths.fluxAnalysis,'biomass_shadow_prices'); % Folder with pan microbe biomass shadow prices 
shadowPriceDir = paths.shadowPriceDir;
fluxPath = paths.fluxPath;
mappedMicrobePath = paths.mappedMicrobePath;
saveDir = paths.microbetoflux;

% Start parallel pool if needed
poolobj = gcp('nocreate'); if isempty(poolobj); parpool(feature('numCores')); end

mContributionDir = shadowPriceDir;

% Find the associated microbial contributors to the fluxes
tic
[elasticNetResults, elasticNetStats, lassoResTab, microbeContributionStats] = findMicrobialContributors(...
    shadowPriceDir, fluxPath, mappedMicrobePath, saveDir, param);
toc
%%
% Microbe associations with AD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Step 2: Associate the identified microbial species with AD progression
% Save the identified flux-associated microbes
microbesToTest = unique(elasticNetResults.Taxa);
microbesToTest(matches(microbesToTest,'Sex'))=[];
paths.adFluxAssociatedMicrobes = microbesToTest';

% Associated relative abundances of microbes to APOE status
[preparedInputTable, preparedMetadata] = prepareDataForStatsADRC(paths.mappedMicrobePath, paths.metadataPrunedProcessed, true); % Load microbiome read data

%%

% How many BMI samples were available


%%

if 0 
preparedInputTable = preparedInputTable(:,[{'ID'},paths.adFluxAssociatedMicrobes]); % Filter on microbes of interest
%preparedInputTable(:,2:end) = fillmissing(preparedInputTable(:,2:end),'constant',0); % Do I need this?
end
paths.adConfounders = {'Sex','age_at_collection','lane','mapped_species_reads','Ethanol_added'}; 
paths.adConfounders = {'Sex','age_at_collection','EDUC','NACCBMI','mapped_species_reads','lane','Ethanol_added','APOE_E4'};
if 0
metadata.NACCBMI = fillmissing(preparedMetadata.NACCBMI,"knn",5);
preparedMetadata.EDUC = fillmissing(preparedMetadata.EDUC,"knn",5);
end
% Test if the flux associated microbes individually associate with the cognitive scores
preparedInputTable = preparedInputTable(:,{'ID','Bacteroides_thetaiotaomicron','Bacteroides_uniformis'});
[results_AD_microbe, paths.adMicrobesOfInterest, results_DM_microbe, regressions_DM_microbe] = alzheimerAnalysis(preparedInputTable,preparedMetadata,'Flux',paths.adConfounders);

% Create regression table for figure with microbes:
resTableMicrobesAD = buildRegressionTables(paths.adMicrobesOfInterest, results_DM_microbe, regressions_DM_microbe, 'AD');

% Set function parameters
param.response = "NACCUDSD"; param.titleAnnotation = " relative abundances"; param.yTitle = {'Z-scaled and log2 transformed','relative abundance'}; param.addEmptyTile=false;
param.italicTitle = true;
% Prepare metadata for visualisation
preparedMetadata.(param.response) = categorical(preparedMetadata.(param.response),{'Healthy controls','MCI','AD dementia'});

% Create tiled figure
close all; fig = figure('Position',[686,152.3333,1231,411.6667]);
tiledlayout(1,length(paths.adMicrobesOfInterest),'TileSpacing','tight','Padding','loose');
createMultipleBoxPlotsADRC(preparedInputTable, preparedMetadata, paths.adMicrobesOfInterest, results_DM_microbe, param);

% Update plot labels
fig.Children.Children(2).Children(1).String = 'B';
fig.Children.Children(1).Children(1).String = 'C';

% Generate paths to save results
adMicrobeBoxPlotPath = fullfile(paths.figures,'adMicrobeBoxPlot.png');
adRegTabSumPath = fullfile(paths.AD,'AD_RegressionTableSummary.xlsx');
filePath = fullfile(paths.AD,'AD_progression_results.xlsx');

% Write results to file
writetable(resTableMicrobesAD,adRegTabSumPath,'Sheet','Microbes')
writetable(results_DM_microbe,filePath,'Sheet','Microbes')
exportgraphics(fig,adMicrobeBoxPlotPath,'Resolution',300)
%%

% Investigate whether these two microbial species also associate with
% cognitive status in the raw gut microbiome datasets. Are the results
% robust after mapping?
[~, preparedMetadata] = prepareDataForStatsADRC(paths.mappedMicrobePath, paths.metadataPrunedProcessed, true); % Load microbiome read data

microbiome = rows2vars(readtable(paths.microbiotaPath,'ReadRowNames',true),'VariableNamingRule','preserve');
microbiome = renamevars(microbiome,'OriginalVariableNames','ID');
microbiome = microbiome(:,{'ID','Bacteroides_thetaiotaomicron','Bacteroides_uniformis'});
microbiome = convertvars(microbiome,{'Bacteroides_thetaiotaomicron','Bacteroides_uniformis'},@normalize);

% Perform regressions
[results_AD_microbe_raw, microbesOfInterest, results_DM_microbe_raw, regressions_DM_microbe_raw] = alzheimerAnalysis(microbiome,preparedMetadata,'Flux',paths.adConfounders);

% Set function parameters
param.response = "NACCUDSD"; param.titleAnnotation = " relative abundances"; param.yTitle = {'Z-scaled','relative abundance'}; param.addEmptyTile=false;
param.italicTitle = true;
% Prepare metadata for visualisation
preparedMetadata.(param.response) = categorical(preparedMetadata.(param.response),{'Healthy controls','MCI','AD dementia'});

% Create tiled figure
close all; fig = figure('Position',[686,152.3333,1231,411.6667]);
tiledlayout(1,length(microbesOfInterest),'TileSpacing','tight','Padding','loose');
createMultipleBoxPlotsADRC(microbiome, preparedMetadata, microbesOfInterest, results_DM_microbe_raw, param);

% Generate paths to save results
adMicrobeBoxPlotPath = fullfile(paths.figures,'premapped_adMicrobeBoxPlot.png');
adRegTabSumPath = fullfile(paths.AD,'AD_RegressionTableSummary.xlsx');
filePath = fullfile(paths.AD,'AD_progression_results.xlsx');

% Write results to file
writetable(resTableMicrobesAD,adRegTabSumPath,'Sheet','premapped_Microbes')
writetable(results_DM_microbe,filePath,'Sheet','premapped_Microbes')
exportgraphics(fig,adMicrobeBoxPlotPath,'Resolution',600)


%%


%%

% Microbe associations with APOE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First, select the flux-associated microbes. Then, associate their
% relative abundances with APOE genotype. 

% Save the identified flux-associated microbes
microbesToTest = elasticNetResults.Taxa(~matches(elasticNetResults.Taxa,'Sex'));
paths.apoeFluxAssociatedMicrobes = microbesToTest';

% Associated relative abundances of microbes to APOE status
[preparedInputTable, preparedMetadata] = prepareDataForStatsADRC(paths.mappedMicrobePath, paths.metadataPrunedProcessed, true); % Load microbiome read data
preparedInputTable = preparedInputTable(:,[{'ID'},paths.apoeFluxAssociatedMicrobes]); % Filter on microbes of interest
%preparedInputTable(:,2:end) = fillmissing(preparedInputTable(:,2:end),'constant',0); % Do I need this?

% Test if the flux associated microbes individually associate with APOE
% status
results_APOE_microbe = apoeAnalysis(preparedInputTable,preparedMetadata, paths.apoeFluxConfounders, 0.05);

if 0 % No visualisation needed for now
    % Set function parameters
    param.response = "APOE_ALLELE"; param.titleAnnotation = " relative abundances"; param.yTitle = "Normalised log2 relative abundance"; param.addEmptyTile=false;
    param.rxnsOfInterest = paths.apoeFluxAssociatedMicrobes;
    
    
    % Prepare metadata for visualisation
    preparedMetadata.(param.response) = categorical(preparedMetadata.(param.response),{'E2','E3','E4'});
    
    
    % Create figure;
    microbesToPlot = unique(results_APOE_microbe.Reaction);
    [plt,ax] = createMultipleBoxPlotsADRC(preparedInputTable, preparedMetadata, microbesToPlot, results_APOE_microbe, param);
end

% Write microbe APOE results to file
apoeRegTabPath = fullfile(paths.apoe,'APOE_regressions.xlsx');
writetable(results_APOE_microbe,apoeRegTabPath,'Sheet','Microbes')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Differential plasma metabolomics %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Dementia plasma metabolomics associations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Next, we will perform the same regression analysis on the plasma
% metabolomic data. Lets check if global cognition does associate with
% plasma metabolomic measurements

% Load and prepare input data and metadata
transformData = false;
[preparedInputTable, preparedMetadata] = prepareDataForStatsADRC(paths.metabolonPath, paths.metadataPrunedProcessed, transformData);

% Filter on metabolites of interest
varsToCheck = ['ID', erase(paths.adRxnsOfInterest,{'DM_','[bc]'})'];
preparedInputTable = preparedInputTable(:, matches(preparedInputTable.Properties.VariableNames,varsToCheck) );

% Get confounders for plasma metabolomics
metabolomeConfoundersAD = paths.adConfounders;
metabolomeConfoundersAD(matches(metabolomeConfoundersAD,{'mapped_species_reads','lane','Ethanol_added'}))=[];

% Perform analysis on metabolomics
[results_AD, plasmaMetsOfInterest, results_DM_plasma,regressions_DM_plasma] = alzheimerAnalysis(preparedInputTable,preparedMetadata,'Flux',metabolomeConfoundersAD);

% Rename results for clarity
results_DM_plasma.Formula = replace(results_DM_plasma.Formula,'Flux','Plasma_metabolite');
results_DM_plasma.Predictor = replace(results_DM_plasma.Predictor,'Flux','Plasma_metabolite');

if 0 % Not included in the paper
    % Define metabolites to plot
    paths.adPlasmaRxnsOfInterest = unique(results_DM_plasma.Reaction(results_DM_plasma.pValue<0.05),'stable');
    
    % Create regression table:
    % resTablePlasma = buildRegressionTables(paths.adPlasmaRxnsOfInterest, results_DM_plasma, regressions_DM_plasma, 'AD');
    % writetable(resTablePlasma,apoeRegTabSumPath,'Sheet','Plasma')
    
    % Create figure of plasma metabolomics results
    % Set function parameters
    param.response = "NACCUDSD"; param.titleAnnotation = " plasma abundances"; param.yTitle = "Normalised log2 abundance";
    
    % Prepare apoe group data
    preparedMetadata.NACCUDSD = categorical(preparedMetadata.NACCUDSD,{'normal cognition','MCI','Dementia'});
    
    % Create tiled figure
    fig = figure('Position',[-1344,286,1051,458]);
    tiledlayout(1,length(paths.adPlasmaRxnsOfInterest),'TileSpacing','tight','Padding','loose');
    createMultipleBoxPlotsADRC(preparedInputTable, preparedMetadata, paths.adPlasmaRxnsOfInterest, results_DM_plasma, param);
    
    adPlasmaBoxplotPath = fullfile(paths.figures,'AD_plasma_boxplots.png');
    exportgraphics(fig,adPlasmaBoxplotPath,'Resolution',300)
end

% Save the flux and metabolomic results
filePath = fullfile(paths.AD,'AD_progression_results.xlsx');
writetable(results_DM_plasma,filePath,'Sheet','Plasma','WriteMode','overwritesheet')


% APOE plasma metabolomics associations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Next, we will perform the same regression analysis on the plasma
% metabolomic data. Lets check if global cognition does associate with
% plasma metabolomic measurements

% First, lets load and process the plasma metabolomic data. I want to
% investigate all mapped metabolomic markers for now.
transformData = false;
[preparedInputTable, preparedMetadata] = prepareDataForStatsADRC(paths.metabolonPath, paths.metadataPrunedProcessed, transformData);
preparedMetadata = preparedMetadata(matches(preparedMetadata.NACCUDSD,'Healthy controls'),:); % Only investigate cognitively normal samples

% Remove microbiome-specific covariates
metabolomeConfounders = paths.apoeFluxConfounders;
metabolomeConfounders(matches(metabolomeConfounders,{'mapped_species_reads','lane','Ethanol_added'}))=[];

% Perform analysis on metabolomics
pFilter=0.05;
[results_APOE_plasma,regressions_APOE_plasma] = apoeAnalysis(preparedInputTable,preparedMetadata, metabolomeConfounders, pFilter);

% Rename results for clarity
results_APOE_plasma.Formula = replace(results_APOE_plasma.Formula,'Flux','plasma_met');
results_APOE_plasma.Predictor = replace(results_APOE_plasma.Predictor,'Flux','plasma_met');

% Define metabolites to plot
paths.apoePlasmaRxnsOfInterest = unique(results_APOE_plasma.Reaction(results_APOE_plasma.pValue<0.05),'stable');

% Create regression table for figure 1:
resTablePlasma = buildRegressionTables(paths.apoePlasmaRxnsOfInterest, results_APOE_plasma, regressions_APOE_plasma, 'APOE');
if 0    
    % Create figure of plasma metabolomics results
    % Set function parameters
    param.response = "APOE_ALLELE"; param.titleAnnotation = " plasma abundances"; param.yTitle = "Normalised log2 abundance";
    
    % Prepare apoe group data
    preparedMetadata.APOE_ALLELE = categorical(preparedMetadata.APOE_ALLELE,{'E2','E3','E4'});
    
    % Create tiled figure
    fig = figure('Position',[39,457,1878,421]);
    tiledlayout(1,length(paths.apoePlasmaRxnsOfInterest),'TileSpacing','tight','Padding','loose');
    createMultipleBoxPlotsADRC(preparedInputTable, preparedMetadata, paths.apoePlasmaRxnsOfInterest, results_APOE_plasma, param);
    
    % Generate path to save results
    apoePlasmaBoxplotPath = fullfile(paths.figures,'APOE_plasma_boxplots.png');

    exportgraphics(fig,apoePlasmaBoxplotPath,'Resolution',300)

end    

% Write results to file

% Save the flux and metabolomic results
apoeRegTabPath = fullfile(paths.apoe,'APOE_regressions.xlsx');
apoeRegTabSumPath = fullfile(paths.apoe,'apoeRegressionTableSummary.xlsx');
writetable(results_APOE_plasma,apoeRegTabPath,'Sheet','Plasma')
writetable(resTablePlasma,apoeRegTabSumPath,'Sheet','Plasma')
    


%% Supplementary tables

% Get the supplementary tables 
% Clean up workspace
clearvars -except paths; clc;

% Table 1: Mapped and unmapped taxa + mapping coverages

% Load the present and absent species and concatinate the results. 
% Add a one for all the mapped species and a zero for all the unmapped
% ones. Then, add taxonomic and phylogenetic information.
clc
suplFilePath = fullfile(paths.outputs,'Supplementary_tables.xlsx');
if isfile(suplFilePath); delete(suplFilePath); end % For script testing

% Preset variable with table information


%%%% MICROBIOTA AND MAPPING STATISTICS ->

% Set inputs for species-level summary statistics
paths.species = fullfile(paths.outputs,'resultMARS', 'metrics','Species'); % Set directory with MARS results
fileNames = {'preMapping_abundanceMetrics_Species.csv',... % Define the files to load
    'mapped_abundanceMetrics_Species.csv',...
    'unmapped_abundanceMetrics_Species.csv'};

% Generate supplementary table
mergedTaxa = collectAbundanceStatsADRC(paths, fileNames);

% Set inputs for phylum-level summary statistics
taxonomyPath = fullfile(paths.outputs,'resultMARS','preprocessedInput_afterRenaming.csv');
microbiotaPath = paths.microbiotaPathFilt;
microbiotaWbmPath = paths.microbiotaWbmPathFilt;

% Filter on all gut microbes in the filtered gut microbiome table
microbiomeFiltNames = readcell(paths.microbiotaPathFilt,'Range','B1:ZZZ1')';
mergedTaxa = mergedTaxa(matches(mergedTaxa.("Microbial species"),microbiomeFiltNames),:);

% Collect summary statistics on phylum-level mapping effects
summaryMerged = collectPhylumStatsADRC(microbiotaPath, taxonomyPath, microbiotaWbmPath);

% Load the unmapped characterised microbial taxa and their increase in read
% coverage if mapped

unmappedCharTaxaPath = fullfile(paths.outputs,'ResultMars','metrics','Species','potential_higher_readCov_when_mapping_more_microbes.xlsx');
unmappedCharTaxa = readtable(unmappedCharTaxaPath,'Sheet','characterised_taxa');
unmappedCharTaxa = removevars(unmappedCharTaxa, {'potential_Incr','minimum','maximum','non_zero_count'});
unmappedCharTaxa = renamevars(unmappedCharTaxa,{'Taxon','mean','SD','curr_readCov','potential_readCov','cumIncrease'},{'Unmapped microbial species','Mean relative abundance','Standard deviation relative abundance','Current mapped read coverage','Read coverage if taxa was mapped','Cumulative increase in read mapping coverage'});

% Save summary statistics to table
description = cell(2,1); 
description{1} = 'Mean relative abundances of mapped and unmapped microbial phyla and species'; % Header
description{2} = ['Left table: Microbial species mapping status and mean (SD) relative abundances across the cohort. ',...
    'Middle table: Microbial phylum-level mean (SD) relative abundances across the cohort for the total, mapped, and unmapped microbial species.',...
    'Right table: Unmapped, but fully characterised microbial species and the increase in mapping coverage if they were mapped.']; % Details
folder = paths.outputs;
suplTable = mergedTaxa;
writeSupplementADRC(mergedTaxa, description, paths.outputs, summaryMerged, unmappedCharTaxa)
%%

%%%% Microbiome-WBM content statistics ->

% Get model summary statistics
if 1
    tic
    disp('Start model content analysis')
    outputFolder = paths.outputs;
    numWorkersCreation = feature('numCores');
    [modelStats, summaryStats, paths.mWbms] = getPdMicrobiomeWBMstatsADRC(paths.outputs, feature('numCores')); % 
    toc % 168 seconds with 18 parallel workers on an intel core i9-10890xe
else
    summaryStats = array2table({'Hello world'},'VariableNames',{'Variable'});
end


% Save model content information to table
description = cell(2,1); 
description{1} = 'Summary statistics of microbiome-WBM model content'; % Header
description{2} = ''; % Details
writeSupplementADRC(summaryStats, description, paths.outputs)


%%%% DIET ->

% Find model diet and save to supplementary tables
dietTable = collectDietInfoADRC(paths.mWbms);
% Save dietary information to table
description = cell(2,1); 
description{1} = 'Diet makeup of host-microbiome WBMs'; % Header
description{2} = 'Dietary metabolites taken from the average European diet (DOI: http://dx.doi.org/10.1093/nar/gky992)'; % Details
writeSupplementADRC(dietTable, description, paths.outputs)


%%%% INVESTIGATED METABOLITES ->
% paths.mWbms = mWBMPath;
% Annotate metabolites for which fluxes were predicted
metPath = fullfile(paths.inputs,'ADRC_metabolites.txt');
[metabolitesCheck, metaboliteSummaryStats, metaboliteOntologyStats] = annotateInvestigatedMetabolitesADRC(metPath,paths.mWbms);
% Save metabolite information to table
description = cell(2,1); 
description{1} = 'Description of investigated metabolites'; % Header
description{2} = 'Metabolites for which fluxes in the WBM blood compartments were predicted'; % Details
writeSupplementADRC(metabolitesCheck, description, paths.outputs, metaboliteSummaryStats)


%%%% AD REGRESSION RESULTS FOR THE FLUXES->

% Load and prepare regression tables
fluxAdRegressionTable = prepareFluxRegressionsForSM(paths);

% Save tables to supplementary file
description = cell(2,1); 
description{1} = 'Logistic regression outcomes of predicted fluxes against AD status on the predicted fluxes. '; % Header
description{2} = [...
    'The regression log odds represent the estimated change in log odd probability between the control and the case group, e.g., MCI vs Dementia, with an increase in predicted flux of one standard deviation.', ...
    'Positive regression coefficients indicate positive associations between predicted fluxes from the less severe disease status to the more severe disease status, e.g., control -> AD MCI, while negative regression coefficients indicate negative correlations with the fluxes.']; % Details
writeSupplementADRC(fluxAdRegressionTable, description, paths.outputs)

%%%% PLASMA METABOLITES, FLUX-PLASMA CORRELATIONS, AND PLASMA REGRESSION RESULTS ->

% Load tables
[plasmaMetabolites,~, plasmaAdRegressionTable] = preparePlasmaTablesForSM(paths);

% Save tables to supplementary file
description = cell(2,1); 
description{1} = 'Left table: Measured plasma metabolites and their associated analysed VMH IDs. Right table: Logistic regression outcomes of predicted fluxes against AD status on the analysed plasma metabolites'; % Header
description{2} = ['The regression log odds represent the estimated change in log odd probability between the control and the case group, e.g., MCI vs Dementia, with an increase in normalised plasma levels of one standard deviation.', ...
    'Positive regression coefficients indicate positive associations between plasma levels from the less severe disease status to the more severe disease status, e.g., control -> AD MCI, while negative regression coefficients indicate negative correlations with the plasma levels.',...
    ]; % Details
writeSupplementADRC(plasmaMetabolites, description, paths.outputs, plasmaAdRegressionTable)



%%%% FLUX MICROBE ASSOCIATIONS ->

% Load summarised shadow prices and microbe selection frequencies 
[bootSpRes,selFreq] = prepareFluxMicrobeLinksForSM(paths);

% Save tables to supplementary file
description = cell(2,1); 
description{1} = ['Left table: Bootstrapped (N=10,000) mean and 95% confidence interval (CI) of shadow price values for the flux-associated microbial pan biomass compounds. Metabolite-microbe pairs with a bootstrapped 95% CI of their mean averages that crosses zero are indicated by a one in the last column.',...
    'Right table: Selection frequencies from LASSO regressions showing the the fractions of non-zero beta values over a subset of 500 random subsets comprising of 70% of the full cohort samples. Missing values indicate microbial species that were not associated with the flux predictions.']; % Header
description{2} = '';
writeSupplementADRC(bootSpRes, description, paths.outputs, selFreq)


%%%% ELASTIC NET REGRESSION RESULTS AND MICROBIAL AD ASSOCIATIONS ->

% Load and process the elastic net and microbial AD regression results
[enetRes, microbesAdRegressionTable] = prepareMicrobeResultsForSM(paths);

% Save tables to supplementary file
description = cell(2,1); 
description{1} = 'Left table: Standardised coefficients from elastic net regressions. Right table: Logistic regression associations for the microbial species relative abundances against AD status (Control vs MCI, MCI vs AD dementia, and Control vs AD dementia).';
description{2} = '';    
writeSupplementADRC(enetRes, description, paths.outputs, microbesAdRegressionTable)

%%%% APOE AND MOCA REGRESSION RESULTS FOR FLUXES->

% Load and prepare regression tables
[~, fluxapoeRegressionTable, fluxMocaRegressionTable] = prepareFluxRegressionsForSM(paths);

% Save tables to supplementary file
description = cell(2,1); 
description{1} = ['Left table: ANCOVA regression table for APOE status (E2 versus E3, E3 versus E4, E2 versus E4) against the flux predictions (response).',...
'Right table: Linear regression table of the flux predictions (predictor) against the Z-scaled MoCA scores.']; % Header
description{2} = 'Both the ANCOVA and linear regressions were performed on healthy individuals for flux predictions of metabolites with altered fluxes in AD patients compared to controls';
writeSupplementADRC(fluxapoeRegressionTable, description, paths.outputs, fluxMocaRegressionTable)


%%% Flux sensitivity upon dietary changes %%%

% Obtain flus sensitivity statistics

fbaDir = paths.FBA;
saveDir = paths.fluxes;
bootSamp = 10000;

poolobj = gcp('nocreate');
if isempty(poolobj); parpool(feature('numCores')); end

% Extract diet reaction reduced cost values
disp('Extract diet reaction reduced cost values')
tic
dietSensFolder = getDietSensitivity(fbaDir, saveDir);
toc
% Calculate the mean average reduced cost value and 95% CI using
% bootstrapping
disp('Calculate the bootstrapped mean average reduced cost for each objective and diet reaction')
tic
bootMeanTable = getDietRCstats(dietSensFolder, bootSamp);
toc

% Convert diet reaction names to metabolite names via VMH database lookup
database = loadVMHDatabase().metabolites;
dietMets = erase(bootMeanTable.('Diet reaction'), {'Diet_EX_', '[d]'});

[uniqueMets, ~, idx] = unique(dietMets);
metNames = repmat({'NA'}, size(uniqueMets));
[~, ia, ib] = intersect(uniqueMets, database(:,1));
metNames(ia) = database(ib, 2);

bootMeanTable.('Diet metabolite name') = renamecats(categorical(dietMets), uniqueMets, metNames);

% Save tables to supplementary file
description = cell(2,1); 
description{1} = ['Mean reduced cost values for each pair of dietary reactions and demand reaction objectives in the cohort. ',...
'Mean average robust cost values across all samples with corresponding 95% confidence intervals and P-values indicating whether the mean robust cost value is non-zero on average. The 95% confidence intervals and P-values were estimated by bootstrapping using 10,000 random samples with replacement.']; % Header
description{2} = '';
writeSupplementADRC(bootMeanTable, description, paths.outputs)


% Save results
paths.dietSensStats = fullfile(paths.fluxes,'diet_redCost_sampleMean_stats.csv');
writetable(bootMeanTable, paths.dietSensStats) % Write results to file
%%% END OF FILE %%%