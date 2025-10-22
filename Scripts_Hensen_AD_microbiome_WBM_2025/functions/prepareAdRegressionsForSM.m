function [fluxAdRegressionTable, plasmaAdRegressionTable] = prepareAdRegressionsForSM(paths)
% Function aim: Load the flux regression results for AD and the serum
% metabolomic results for AD and prepare the tables for supplementary
% materials.
%
% input:
% paths
%
% output:
% fluxAdRegressionTable
% plasmaAdRegressionTable


% Generate path to .xlsx file
adPath = fullfile(paths.AD,"AD_progression_results.xlsx");

%%% FLUXES %%%

% Load flux data
fluxAdRegressionTable = readtable(adPath,'VariableNamingRule','preserve','Sheet','Fluxes');

% Convert last column into logical
fluxAdRegressionTable{:,end} = logical(fluxAdRegressionTable{:,end});

% Store variable names
oldNames = fluxAdRegressionTable.Properties.VariableNames;

% Define new names
newNames = {...
    'Maximised WBM reaction',...
    'Regression formula',...
    'Predictor of interest',...
    'Distribution of dependent variable',...
    'Analysed groups',...
    'Analysed samples',...
    'Estimated log odds coefficient',...
    '2.5% CI of log odds coefficient',...
    '97.5% CI of log odds coefficient',...
    'Standard error',...
    'Regression t-statistic',...
    'P-value',...
    'FDR BH',...
    'Regression R-squared',...
    'Consistent change from CN to MCI and DEMENTIA'};

% Apply improved names
fluxAdRegressionTable = renamevars(fluxAdRegressionTable,oldNames,newNames);

%%% Metabolomics %%%

% Load table
plasmaAdRegressionTable = readtable(adPath,'VariableNamingRule','preserve','Sheet','Plasma');

% Convert last column into logical
plasmaAdRegressionTable{:,end} = logical(plasmaAdRegressionTable{:,end});

% Get variable names
oldNames = plasmaAdRegressionTable.Properties.VariableNames;

% Update new variable names for plasma metabolites
newNamesPlasma = newNames;
newNamesPlasma(1) = {'Plasma metabolite'};

% Update table variable names
plasmaAdRegressionTable = renamevars(plasmaAdRegressionTable,oldNames,newNamesPlasma);

end
