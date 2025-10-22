function [bootSpRes,selFreq] = prepareFluxMicrobeLinksForSM(paths)
% Function aim: Load bootstrapped microbial biomass shadow price values for
% the predicted fluxes and load the microbial selection frequencies from
% the LASSO regression stability selection. 
%
% input:
% paths
%
% output:

%%% Shadow prices %%%

% Load shadow price data
bootSpPath = fullfile(paths.microbetoflux,'bootciMicrobeSPvalues.csv');
bootSpRes = readtable(bootSpPath,'VariableNamingRule','preserve');

% Rename variables
bootSpRes = renamevars(bootSpRes,{'Metabolite','Taxa','Mean'},{'Maximised WBM reaction','Microbial species biomass','Bootstrap mean shadow price'});

%%% Microbe selection frequencies %%%

% Load selection frequency data
selFreqPath = fullfile(paths.microbetoflux,'lassoMicrobeSelectionFreq.csv');
selFreq = readtable(selFreqPath,'VariableNamingRule','preserve');

% Rename variable
selFreq = renamevars(selFreq,'Species','Microbial species');
end
