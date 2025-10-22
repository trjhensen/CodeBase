% ADRC mWBM creation configuration file

%% 0. Set required inputs

% Only run mgpipe and createBatchWBM
paths.seqC.flagSeqC = false;
paths.Mars.flagMars = false;
paths.mWBM.usePersonalisedWBM = false;
paths.fba.flagFBA = false; 
paths.stats.flagStatistics = false; 

% Full path where you would like your results stored
resultPath = fullfile(what('ADRC').path,'outputs'); 

% Diet for constraining models (microbiome, personalized WBMs or mWBM-WBMs)
paths.General.diet = 'EUAverageDiet';

% Path to metadata
paths.General.metadataPath = fullfile(resultPath,'processedMetadata.csv');

% Choose your solver
paths.General.solver = 'gurobi';

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
paths.mgPipe.flagMgPipe = true; 

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


% Add other necessary flags
paths.persWBM.flagPersonalise = false;
paths.mWBM.alteredWBMPath = 'notempty';
paths.seqC.outputPathSeqC = '';
paths.Mars.outputPathMars = '';
%paths.mgPipe.outputPathMgPipe = '';
paths.persWBM.outputPathPersonalisation = '';
%paths.mWBM.outputPathMWBM = '';
paths.fba.outputPathFluxResult = '';
paths.fba.outputPathFluxAnalysis = '';
paths.stats.outputPathStatistics  = '';