function [modelStats, summaryStats, mWBMPath] = getPdMicrobiomeWBMstatsADRC(outputFolder, numWorkersCreation)
% This function loads data from the generated microbiome-WBMs and calculates the
% number of reactions, metabolites, constraints, and unique taxa per WBM.
% Mean averages + SD are also obtained for each statistic. 
%
% INPUTS: 
% mWBMPath                   Path to directory where the HM models are
%                               saved. 
% numWorkersCreation                    Number of cores used for parallelisation.
%                               Default = 4.
% 
% OUTPUTS
% modelStats                    Table with summary statistics on the generated WBMs:
%                               gender, number of reactions, metabolites,constrainsts, 
%                               and taxa.
% summaryStats                  Table with the mean and SD of the model
%                               statistics in the modelStats variable.
%
% Author: Tim Hensen, July 2024

if nargin < 2
    numWorkersCreation = 4;
end

% Create new central folder for mWBMs
mWBMPath = fullfile(outputFolder, 'mWBMmodels');
if ~isfolder(mWBMPath); mkdir(mWBMPath); end

% Move models from the knirps and sneezy folders to the centralised folder
knirpsFluxPath = fullfile(outputFolder,'Knirps','outputs','mWBMmodels');
sneezyFluxPath = fullfile(outputFolder,'Sneezy','outputs','mWBMmodels');

% Find current model paths
knirpsWbmPaths = fullfile(knirpsFluxPath,what(knirpsFluxPath).mat);
sneezyWbmPaths = fullfile(sneezyFluxPath,what(sneezyFluxPath).mat);

% Generate new model paths
newKnirpsWbmPaths = fullfile(mWBMPath,what(knirpsFluxPath).mat);
newSneezyWbmPaths = fullfile(mWBMPath,what(sneezyFluxPath).mat);

% Copy WBMs to centralised path
cellfun(@copyfile, knirpsWbmPaths, newKnirpsWbmPaths)
cellfun(@copyfile, sneezyWbmPaths, newSneezyWbmPaths)

% Get paths to models
paths = fullfile(mWBMPath,string(what(mWBMPath).mat));

% Set parallel pool
if numWorkersCreation > 1
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        parpool(numWorkersCreation)
    end
end

% Preallocate variables for WBM statistics
IDs = string(zeros(length(paths),1));
Sex = string(zeros(length(paths),1));
numRxns = nan(length(paths),1);
numMets = nan(length(paths),1);
numConstraints = nan(length(paths),1);
numTaxa = nan(length(paths),1);
numDietMets = nan(length(paths),1);

parfor i=1:length(paths)
    % load data from WBM
    data = load(paths(i),'rxns','mets','ctrs','ID','sex','lb');

    % Get data
    IDs(i) = string(data.ID);
    Sex(i) = string(data.sex);
    numRxns(i) = length(data.rxns);
    numMets(i) = length(data.mets);
    numConstraints(i) = length(data.ctrs);
    numTaxa(i) = sum(contains(data.rxns,'biomassPan'));
    numDietMets(i) = sum(contains(data.rxns,'Diet_') & data.lb<0);
end

% Remove parallel pool
delete(poolobj)

% Organise data to table
modelStats = table(IDs, Sex, numRxns, numMets, numDietMets, numConstraints, numTaxa);
modelStats.Properties.VariableNames = {'ID','sex','Reactions','Metabolites','Dietary metabolites','Constraints','Taxa'};

% Calculate summary statistics for male and female models separately

% Get numerical columns
numData = modelStats{:,vartype('numeric')};

% Find male and female samples
[groups,names] = findgroups(modelStats.sex);

% Calculate mean and SD stats
meanStats = array2table(splitapply(@mean,numData,groups)','VariableNames', append("Mean ",names, " samples"));
stdStats = array2table(splitapply(@std,numData,groups)','VariableNames', append("SD ",names, " samples"));

% Create table for summary statistics
summaryStats = [meanStats,stdStats];

% Add variable information to the first column
summaryStats = addvars(summaryStats, ...
    modelStats(:,vartype('numeric')).Properties.VariableNames', ...
    'NewVariableNames','Statistic',...
    'Before',1);
end