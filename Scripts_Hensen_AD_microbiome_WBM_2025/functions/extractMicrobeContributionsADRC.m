function [spResPaths, contrResPaths] = extractMicrobeContributionsADRC(fbaFolder, saveDir)
% Efficienty extracts the biomass shadow prices and relative abundances in
% order to calculate the potential microbe flux contributions.
% 
% INPUTS:
% fbaFolder = paths.FBA;
% saveDir = paths.fluxAnalysis;
% 
% OUTPUTS:
% spResPaths        Paths to biomass shadow price results
% contrResPaths     Paths to potential microbial contribution results

% Find paths to fba results
solPaths = fullfile(what(fbaFolder).path,what(fbaFolder).mat);
solPaths(contains(solPaths,'_gf_'))=[];

% Load fba solutions
solutions = cellfun( @(x) load(x,'ID','rxns','taxonNames','shadowPriceBIO','relAbundances'), solPaths);

% Find all microbial species
allSpecies = unique(vertcat(solutions(:).taxonNames));

% Preallocate biomShadowPrices array for shadow prices
numRxns = size(solutions(1).rxns,1);
numIDs = height(solutions);
numTaxa = length(allSpecies);

biomShadowPrices = nan(numIDs, numTaxa, numRxns);
relAbun = nan(numIDs, numTaxa);
microbeContributions = nan(numIDs, numTaxa, numRxns);

for i = 1:size(solutions, 1)
    
    % Find which species in the sample are present
    [index1, index2] = ismember(allSpecies, solutions(i).taxonNames);

    % Remove zero indices for non-present species in index2
    index2 = index2(index2 > 0);
    
    if ~isempty(index2)
        % Extract data
        samp_SP = solutions(i).shadowPriceBIO; % Extract shadow prices
        samp_ra = solutions(i).relAbundances; % Extract relative abundances

        % Filter microbes in model
        samp_SP = samp_SP(index2, :);
        samp_ra = samp_ra(index2);

        % Calculate microbe flux contribution potentials
        samp_c = samp_SP.*samp_ra;
        
        % Assign the shadow prices to the result array
        biomShadowPrices(i, index1, :) = samp_SP;
        relAbun(i,index1) = samp_ra;
        microbeContributions(i, index1, :) = samp_c * -1;
    end
end

% Obtain microbial biomass shadow prices for the predicted fluxes
createSpTab = @(x) array2table(biomShadowPrices(:,:,x),'RowNames', vertcat(solutions(:).ID), 'VariableNames',allSpecies);
shadowpriceTables = arrayfun(createSpTab,1:numRxns,'UniformOutput',false);

% Obtain total microbial component of the fluxes
createContrTab = @(x) array2table(microbeContributions(:,:,x),'RowNames', vertcat(solutions(:).ID), 'VariableNames',allSpecies);
contributionTables = arrayfun(createContrTab,1:numRxns,'UniformOutput',false);

% Create folders for shadow price and microbe contribution tables
spResFolder = fullfile(saveDir,'biomass_shadow_prices');  if ~isfolder(spResFolder); mkdir(spResFolder); end
contrResFolder = fullfile(saveDir,'potentialMicrobeContributions'); if ~isfolder(contrResFolder); mkdir(contrResFolder); end

% Save microbe contribution tables
spResPaths = fullfile(spResFolder, append(solutions(1).rxns','.csv') );
cellfun(@(x,y) writetable(x,y,'WriteRowNames',true), shadowpriceTables,spResPaths);

contrResPaths = fullfile(contrResFolder, append(solutions(1).rxns','.csv') );
cellfun(@(x,y) writetable(x,y,'WriteRowNames',true), contributionTables,contrResPaths);
end




