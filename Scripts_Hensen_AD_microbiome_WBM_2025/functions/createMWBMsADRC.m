function [modelStats, summaryStats, dietInfo, dietGrowthStats] = createMWBMsADRC()
% Full path where you would like your results stored
resultPath = fullfile(what('ADRC').path,'outputs'); 

% Diet for constraining models (microbiome, personalized WBMs or mWBM-WBMs)
diet = 'EUAverageDiet';

% Choose your solver
solver = 'gurobi';

% Input for function
outputPathMgPipe = fullfile(resultPath,'resultMgPipe','Diet');

% Path to metadata information (used to ensure that metadata is available
% for each mWBM.
metadataPath = fullfile(resultPath,'processedMetadata_processed.csv');

% Numeric value, amount of workers that are to be used for creation of
% microbiome models and mWBMs (we recommend using 80% of possible number of 
% available workers here) 
numWorkersCreation = round(feature('numCores')*0.8);

% Numeric value, amount of workers that are to be used for optimisation of
% microbiome models, iWBMs, mWBMs and/or miWBMs
numWorkersOptimisation = 4;

% Path where the mWBMs are stored
outputPathMWBM = fullfile(resultPath, 'mWBMmodels');

% WBMs to be used
maleUnpersonalisedWBMpath = 'Harvey_1_04c';
femaleUnpersonalisedWBMpath = 'Harvetta_1_04c';

% Create mWBM folder if not present
if ~isfolder(outputPathMWBM); mkdir(outputPathMWBM); end

% Create mWBMs
[modelStats, summaryStats, dietInfo, dietGrowthStats] = createBatchMWBM(outputPathMgPipe, ...
    outputPathMWBM, ...
    metadataPath,...
    'solver', solver,...
    'Diet', diet,...
    'numWorkersCreation',numWorkersCreation,...
    'numWorkersOptimisation',numWorkersOptimisation,...
    'maleUnpersonalisedWBMpath',maleUnpersonalisedWBMpath,...
    'femaleUnpersonalisedWBMpath',femaleUnpersonalisedWBMpath);
end