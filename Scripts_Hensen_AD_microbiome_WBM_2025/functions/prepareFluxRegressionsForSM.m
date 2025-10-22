function [fluxAdRegressionTable, fluxapoeRegressionTable, fluxMocaRegressionTable] = prepareFluxRegressionsForSM(paths)
%function fluxAdRegressionTable = prepareFluxRegressionsForSM(paths)
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


%%% Plasma metabolomics AD %%%
% 
% % Load plasma metabolomics
% plasmaAdRegressionTable = readtable(adPath,'VariableNamingRule','preserve','Sheet','Plasma');
% 
% % Convert last column into logical
% plasmaAdRegressionTable{:,end} = logical(plasmaAdRegressionTable{:,end});
% 
% % Store variable names
% oldNames = plasmaAdRegressionTable.Properties.VariableNames;
% 
% % Rename first column
% newNames{1} = 'VMH ID of plasma metabolite';
% 
% % Update table variable names
% plasmaAdRegressionTable = renamevars(plasmaAdRegressionTable,oldNames,newNames);



%%% APOE %%%

% Generate path to .xlsx file
apoePath = fullfile(paths.apoe,"APOE_regressions.xlsx");

%%% FLUXES %%%

% Load flux data
fluxapoeRegressionTable = readtable(apoePath,'VariableNamingRule','preserve','Sheet','Fluxes');

% Convert last column into logical
fluxapoeRegressionTable{:,end} = logical(fluxapoeRegressionTable{:,end});

% Store variable names
oldNames = fluxapoeRegressionTable.Properties.VariableNames;

% Update table variable names
fluxapoeRegressionTable = renamevars(fluxapoeRegressionTable,oldNames,newNames);


%%% MoCA scores %%%

% Set path
mocaPath = fullfile(paths.cognition,"Cognition_regressions.xlsx");


% Load flux data
fluxMocaRegressionTable = readtable(mocaPath,'VariableNamingRule','preserve','Sheet','Fluxes');

% Convert last column into logical
fluxMocaRegressionTable{:,end} = logical(fluxMocaRegressionTable{:,end});

% Store variable names
oldNames = fluxMocaRegressionTable.Properties.VariableNames;

% Remove last column from new names
newNames = newNames(1:end-1);

% Update table variable names
fluxMocaRegressionTable = renamevars(fluxMocaRegressionTable,oldNames,newNames);

end
