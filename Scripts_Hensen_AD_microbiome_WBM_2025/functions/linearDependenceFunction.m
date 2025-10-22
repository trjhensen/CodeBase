% Linear dependence of reaction abundances.
clear;clc;close all;
% Set path to reaction abundance information:
path = 'C:\Users\mspg\Documents\ADRC\outputs\resultMgPipe\ReactionAbundance.csv';
saveDir = 'C:\Users\mspg\Documents';


% Load table into memory
inputTable = readtable(path,'VariableNamingRule','preserve','ReadRowNames',true);

% Transpose table so that the reactions are the variables
inputTable = rows2vars(inputTable,'VariableNamingRule','preserve');

% Set sample IDs to row names
inputTable.Properties.RowNames = inputTable.OriginalVariableNames;
inputTable = removevars(inputTable,'OriginalVariableNames');

tolerance = 1e-3;


[prunedData, linearDependentVars, filePath] = findLinearlyDependentVariables(inputTable, tolerance, saveDir);

