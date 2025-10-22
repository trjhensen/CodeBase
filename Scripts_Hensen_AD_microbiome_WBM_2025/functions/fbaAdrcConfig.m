% ADRC mWBM creation configuration file

%% 0. Set required inputs

% Perform FBA on the set list of metabolites
paths.fba.flagFBA = true; 

% Do not perform the other parts of Persephone
paths.seqC.flagSeqC = false;
paths.Mars.flagMars = false;
paths.mWBM.usePersonalisedWBM = false;
paths.stats.flagStatistics = false; 

% Full path where you would like your results stored
resultPath ='/home/matlab/Documents/Tim Hensen/ADRC/outputs'; 

% Diet for constraining models (microbiome, personalized WBMs or mWBM-WBMs)
paths.General.diet = 'EUAverageDiet';

% Path to metadata
paths.General.metadataPath = fullfile(resultPath,'processedMetadata.csv');

% Choose your solver
paths.General.solver = 'ibm_cplex';

% Numeric value, amount of workers that are to be used for creation of
% microbiome models and mWBMs (we recommend using 80% of possible number of 
% available workers here) 
paths.General.numWorkersCreation = round(feature('numCores')*0.8);

% Numeric value, amount of workers that are to be used for optimisation of
% microbiome models, iWBMs, mWBMs and/or miWBMs
paths.General.numWorkersOptimisation = 4;


%% 2. MgPipe Inputs

% Boolean, indicates if Part 3 of the pipeline MgPipe/microbiome model 
% creation should be run.
paths.mgPipe.flagMgPipe = false; 

% Assign directory for mgPipe results
paths.mgPipe.outputPathMgPipe = fullfile(resultPath,'resultMgPipe');

% Character array variable to the folder where the output of MARS is stored.
paths.Mars.outputPathMars = fullfile(resultPath,'ResultMars');
paths.Mars.outputExtensionMars = 'csv';
paths.Mars.relAbunFilePath = fullfile(paths.Mars.outputPathMars,'renormalized_mapped_forModelling',['renormalized_mapped_forModelling_species.', paths.Mars.outputExtensionMars]);


% Set path to pan models (pan models are available to download from XXXXX)
paths.mgPipe.microbeReconstructionPath = fullfile(what('ADRC').path,'inputs','reconstructions','ApolloAgora2panSpecies');

% Boolean, indicates if netSecretion and netUptake are calculated for the 
% microbiome models via fastFVA. OPTIONAL, defaults to false
paths.mgPipe.computeProfiles = false;

%% 4. mWBM creation inputs
% Boolean, indicates if Part 4 of the pipeline: mWBM creation and descriptive 
% statistics should be run. OPTIONAL, defaults to true
paths.mWBM.flagMWBMCreation = false; 

% The path where the mWBM models will be stored.
paths.mWBM.outputPathMWBM = fullfile(resultPath, 'mWBMmodels');



%% 5. Flux analysis inputs

% Boolean, indicates if Part 6 of the pipeline: FBA should be run.
paths.fba.flagFBA = true; 

% String with the path where the flux results should be stored.
paths.fba.outputPathFluxResult = fullfile(resultPath,'resultFlux');

% String with the path where the analyses of the flux results are stored.
paths.fba.outputPathFluxAnalysis = fullfile(paths.fba.outputPathFluxResult,'fluxAnalysis');

% Boolean, indicates if the complete .v, .y., and .w vectors are stored in
% the result. OPTIONAL, defaults to true. It is recommended to set saveFullRes.
paths.fba.saveFullRes = true;
paths.fba
% *REQUIRED.* Character array containing reaction IDs for reactions that are to be 
% optimised. If the user wants to solve non-existing demand reactions, 
% simply add DM_ in front of the desired metabolite and add to the rxnList variable.

% Load list of metabolites to be investigated 
[metList, vmhList] = loadADRCMets();
rxnList = append('DM_',vmhList,'[bc]');
%rxnList
paths.fba.rxnList = rxnList;
% Set flux processing parameters

% .NumericalRounding defines how much the predicted flux values are
% rounded. A defined value of 1e-6 means that a flux value of
% 2 + 2.3e-8 is rounded to 2. A flux value of 0 + 1e-15 would be rounded to
% exactly zero. This rounding factor will also be applied to the shadow 
% price values. If microbiome relative abundance data is provided, the same
% rounding factor will be applied to the relative abundance data.
paths.fba.paramFluxProcessing.numericalRounding = 1.00E-06;

% .RxnRemovalCutoff defines the minimal number of samples for which a
% unique reaction flux could be obtained, before removing the reaction for 
% further analysis. This parameter can be expressed as 
% * fraction:  the fraction of samples with unique values, 
% * SD: the standard deviation across samples, and
% * count: the counted number of unique values. If microbiome relative 
% abundance data is provided, the same removal cutoff factor will be 
% applied to the relative abundance data.
paths.fba.paramFluxProcessing.rxnRemovalCutoff = {'fraction', 0.1};

% .RxnEquivalenceThreshold defines the minimal threshold of when 
% functionally identical flux values are predicted, and are thus part of
% the same linear pathways. The threshold for functional equivalence is
% expressed as the R2 (r-squared) value after performing a simple linear 
% regression between two reactions.  
paths.fba.paramFluxProcessing.rxnEquivalenceThreshold = 0.999;

% .fluxMicrobeCorrelationType defines the method for correlating the 
% predicted fluxes with microbial relative abundances. Note that this
% metric is not used if muWBMs are not present. The available correlation
% types are: 
% * regression_r2:  the R2 (r-squared) value from pairwise linear regression on the 
% predicted fluxes against microbial relative abundances.
% * spearman_rho: the correlation coefficient, rho obtained from pairwise
% Spearman nonparametric correlations between predicted fluxes and 
% microbial relative abundances. 
paths.fba.paramFluxProcessing.fluxMicrobeCorrelationMetric = 'spearman_rho';

%% 6. Statistical Analysis
% Logical variable indicating if a statistical analysis should be run on 
% the obtained fluxes and relative taxon abundances against a user-defined
% response variable.
paths.stats.flagStatistics = false; 

% Add other necessary flags
paths.persWBM.flagPersonalise = false;
paths.mWBM.alteredWBMPath = 'notempty';
paths.seqC.outputPathSeqC = '';
paths.Mars.outputPathMars = '';
%paths.mgPipe.outputPathMgPipe = '';
paths.persWBM.outputPathPersonalisation = '';
%paths.mWBM.outputPathMWBM = '';
paths.stats.outputPathStatistics  = '';
paths.Mars.relAbunFilePath = '';