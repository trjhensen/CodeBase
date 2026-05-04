function dietSensFolder = getDietSensitivity(fbaDir, saveDir)
% Aim: The aim of this script is to obtain the reduced cost values of all
% dietary exchange reactions towards each predicted flux. 

% INPUTS:
% fbaDir            Path to directory with FBA results
% saveDir           Path to directory stating where to store the diet
%                   reaction reduced cost values

% Find paths to FBA solutions
solDir = what(fbaDir);
solPaths = string(append(solDir.path, filesep, solDir.mat));
solPaths(contains(solPaths,'_gf_')) = [];

% tmp
% solPaths = solPaths(1:100);

% Load the reduced cost values for each 
numSamp = length(solPaths);

rcTabs = cell(numSamp,1);
IDs = cell(numSamp,1);
parfor i=1:numSamp
    [rcTabs{i},IDs(i)] = getRdTable(solPaths(i)); % Load reduced cost values and ID
end

% Convert tables to 3D-matrix and slice so that for each predicted
% metabolite, the associated ID-W matrix is obtained
rcMat = cellfun(@table2array, rcTabs, 'UniformOutput', false);
rcostTensor = reshape(cat(3, rcMat{:} ), numSamp, size(rcMat{1},1), size(rcMat{1},2));

% Generate tables with reduced cost values to calculate the diet
% sensitivity metrics
createRcTabs = @(x) array2table(rcostTensor(:,:,x),'RowNames', IDs, 'VariableNames',rcTabs{1}.Properties.RowNames');
dietSensitivityTabs = arrayfun(createRcTabs,1:size(rcostTensor,3),'UniformOutput',false);

% Create folders for diet sensitivity metrics
dietSensFolder = fullfile(saveDir,'diet_rxn_reduced_cost');  if ~isfolder(dietSensFolder); mkdir(dietSensFolder); end

% Define paths to folder
dietSensPaths = fullfile(dietSensFolder, append(rcTabs{1}.Properties.VariableNames,'.csv') );
cellfun(@(x,y) writetable(x,y,'WriteRowNames',true), dietSensitivityTabs,dietSensPaths);
end

% Create tables on the reduced cost values for each predicted metabolite
% across all samples and save them

function [rdTAb, ID] = getRdTable(solPath)
% Load solution
solution = load(solPath,'ID','modelRxns','rxns','modelLB','modelUB','w');

% Find diet reactions
drIdx = contains(solution.modelRxns,'Diet_EX_') & (solution.modelLB ~=0 | solution.modelUB ~=0);

% Get rxns and remove unneeded fields
rxns = solution.rxns;
ID = solution.ID;
solution = rmfield(solution,{'ID','rxns','modelLB','modelUB'});

% Filter on diet reactions
solution = structfun(@(x) x(drIdx,:), solution, 'UniformOutput', false);

% Generate rdTable
w = round(full(solution.w) , 6);
rdTAb = array2table(w, 'RowNames',solution.modelRxns,'VariableNames',rxns);
end

