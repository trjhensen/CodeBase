function [enetRes, microbesAdRegressionTable] = prepareMicrobeResultsForSM(paths)
% This function loads and processes the elastic net regression results and
% the microbial relative abundance associations against AD status for the
% supplementary tables

% Load elastic net results
enetResPath = fullfile(paths.microbetoflux,'enetFluxMicrobeRes.csv');
enetRes = readtable(enetResPath,'VariableNamingRule','preserve');

% Load microbial associations with ad status
adPath = fullfile(paths.AD,"AD_progression_results.xlsx");

% Load flux data
microbesAdRegressionTable = readtable(adPath,'VariableNamingRule','preserve','Sheet','Microbes');

% Convert last column into logical
microbesAdRegressionTable{:,end} = logical(microbesAdRegressionTable{:,end});

% Store variable names
oldNames = microbesAdRegressionTable.Properties.VariableNames;

% Define new names
newNames = {...
    'Microbial species',...
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
microbesAdRegressionTable = renamevars(microbesAdRegressionTable,oldNames,newNames);

% Rename predictor of interest cells from 'Flux' to relative abundance
numRows = height(microbesAdRegressionTable.("Predictor of interest"));
microbesAdRegressionTable.("Predictor of interest") = repmat({'Relative abundance'},numRows,1);

% Filter regression results on microbes that associated with the flux
% predictions.
[~,ia,ib] = intersect(microbesAdRegressionTable.('Microbial species'), enetRes.Taxa, 'stable');
microbesAdRegressionTable = microbesAdRegressionTable(ia,:); % Reorder rows
enetRes = enetRes(ib,:); % Reorder rows of elastic net coefficients to be in line with the microbe regression results
end