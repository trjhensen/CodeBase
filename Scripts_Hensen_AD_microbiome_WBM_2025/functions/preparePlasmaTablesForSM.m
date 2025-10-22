function [plasmaMetabolites,plasmaAdRegressionTable] = preparePlasmaTablesForSM(paths)
% Function goal: Load and prepare the measured plasma metabolites and the
% plasma regression results. 

% Get the analysed reactions in the flux results
rxnsToMap = readcell(paths.fluxPath,'Range','1C:1ZZZ');

% Manually map plasma metabolites to VMH ID
[~,plasmaMetabolites,vmhTable3] = mapPlasmaMetsToVmhADRC(rxnsToMap,paths.rawMetabolonPath);

% Add VMH metabolites to plasmaMetabolites table
[~,indexA,indexB] = intersect(plasmaMetabolites.CHEM_ID,vmhTable3.CHEM_ID,'stable');
plasmaMetabolites.('VMH_ID_of_analysed_metabolite')(indexA) = vmhTable3.VMHID(indexB);

% Process metabolite names for supplementary table
plasmaMetabolites = removevars(plasmaMetabolites,{'Var1','TYPE'});
plasmaMetabolites = movevars(plasmaMetabolites,'VMH_ID_of_analysed_metabolite','After','CHEMICAL_NAME');
plasmaMetabolites = sortrows(plasmaMetabolites,"VMH_ID_of_analysed_metabolite","ascend");


%%% AD regression results %%%

% Generate path to .xlsx file
adPath = fullfile(paths.AD,"AD_progression_results.xlsx");

% Load table
plasmaAdRegressionTable = readtable(adPath,'VariableNamingRule','preserve','Sheet','Plasma');

% Convert last column into logical
plasmaAdRegressionTable{:,end} = logical(plasmaAdRegressionTable{:,end});

% Get variable names
oldNames = plasmaAdRegressionTable.Properties.VariableNames;

% Define new names
newNames = {...
    'VMH ID of plasma metabolite',...
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

% Update table variable names
plasmaAdRegressionTable = renamevars(plasmaAdRegressionTable,oldNames,newNames);
end
