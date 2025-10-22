function [OK, allPaths, fbaFolder] = moveAdrcFbaRes(outputFolder,fbaPath)

OK = false;

% Set paths to knirps and sneezy flux results
knirpsFluxPath = fullfile(outputFolder,'Knirps','outputs','resultFlux');
sneezyFluxPath = fullfile(outputFolder,'Sneezy','outputs','resultFlux');

%tst =slimDownFBAresults(knirpsFluxPath);

% Copy the results into a new folder with only the FBA information of interest 
% slimDownFBAresults(knirpsFluxPath);
% slimDownFBAresults(sneezyFluxPath);
% 
% knirpsFluxPath = [knirpsFluxPath '_SLIM'];
% sneezyFluxPath = [sneezyFluxPath '_SLIM'];

% Create flux processing output folder
if ~isfolder(fbaPath); mkdir(fbaPath); end

% Find .mat files in the FBA folder
knirpsFluxes = what(knirpsFluxPath);
knirpsFluxPaths = string(append(knirpsFluxPath, filesep, knirpsFluxes.mat));
newknirpsPaths = string(append(fbaPath, filesep, knirpsFluxes.mat));

% Copy knirps fluxes
arrayfun(@copyfile, knirpsFluxPaths, newknirpsPaths)

% Find .mat files in the FBA folder
sneezyFluxes = what(sneezyFluxPath);
sneezyFluxes.mat(contains(sneezyFluxes.mat,'_gf_')) = []; % Remove gf models

sneezyFluxPaths = string(append(sneezyFluxPath, filesep, sneezyFluxes.mat));
newsneezyPaths = string(append(fbaPath, filesep, sneezyFluxes.mat));

% Copy knirps fluxes
arrayfun(@copyfile, sneezyFluxPaths, newsneezyPaths)

allPaths = vertcat(newknirpsPaths,newsneezyPaths);
fbaFolder = fbaPath;

OK = true;

end