function corrTable = fluxMetabolonCorr2(fluxPath,metabolonPath)
% Correlate all flux results with metabolomics data
% INPUT:
% metabolomicsPath
% fluxpath

% Load data and process
readFun = @(x) readtable(x, 'VariableNamingRule', 'preserve','ReadRowNames', false);
idFun = @(x) convertvars(x, 'ID', @(y) string(erase(y,{'X','mWBM_', '_female', '_male'})));
rmSexFun = @(x) x(:,setdiff(x.Properties.VariableNames,'Sex','stable'));
firstCol2RowFun = @(x) array2table(x{:,2:end},'VariableNames',x.Properties.VariableNames(2:end),'RowNames',x{:,1});
normFun = @(x) convertvars(x, x.Properties.VariableNames, @normalize);

data = cellfun(@(x) normFun(firstCol2RowFun(rmSexFun(idFun(readFun(x))))), {fluxPath, metabolonPath},'UniformOutput',false);
[fluxes, metabolomics] = deal(data{:});

% Remove rows from metabolomics
metabolomics(all(isnan(metabolomics{:,:}),2),:) = [];

% Order sample IDs
fluxes = fluxes(metabolomics.Properties.RowNames,:);
metabolomics = metabolomics(metabolomics.Properties.RowNames,:);

combinedTables = {fluxes,metabolomics};
tabVarNames = cellfun(@(x) x.Properties.VariableNames,combinedTables,'un',0); % Get all reaction and microbe names
corrTable = array2table(zeros(cellfun(@length,tabVarNames)),"VariableNames",tabVarNames{2});
corrTable = addvars(corrTable,tabVarNames{1}','NewVariableNames','Reaction ID','Before',1);
corrTable = stack(corrTable,corrTable.Properties.VariableNames(2:end),"NewDataVariableName",'Rho','IndexVariableName','Plasma metabolite');
corrTable = convertvars(corrTable,["Reaction ID","Plasma metabolite"],'string');

% Get the number of samples in each reation-microbe pair that can
% be correlated               
getSampSize = @(x,y) sum(~isnan(x) & ~isnan(y)); 
nPairs = height(corrTable);
corrN = zeros(nPairs,1);
for i=1:nPairs
    corrN(i) = getSampSize(combinedTables{1}.(corrTable.("Reaction ID")(i)) , combinedTables{2}.(corrTable.("Plasma metabolite")(i)) );
end

% Add sample sizes and preallocate columns for the correlation
% results
corrTable = addvars(corrTable,corrN,zeros(nPairs,1),zeros(nPairs,1),zeros(nPairs,1),'NewVariableNames', {'N','2.5%CI', '97.5%CI','P-value'});
corrTable = movevars(corrTable,'N','Before','Rho');

% Fill empty rows. Needed for identifying the ranks
%preparedTables = cellfun(@(x) fillmissing(x,'constant',0),combinedTables,'UniformOutput',false);
preparedTables = combinedTables;

% Convert matrices to rank for spearman correlations
convToRank = @(x) array2table(tiedrank(table2array(x)),'RowNames',x.Properties.RowNames,'VariableNames',x.Properties.VariableNames);
rankedTables = cellfun(@(x) convToRank(x), preparedTables,'UniformOutput',false);

for i=1:nPairs
    % Perform pairwise spearman correlations
    [rho,p,lower,upper] = corrcoef( rankedTables{1}.(corrTable.("Reaction ID")(i)), rankedTables{2}.(corrTable.("Plasma metabolite")(i)) , 'Rows','pairwise'); 
    corrTable{i,{'Rho','2.5%CI', '97.5%CI','P-value'}} = cellfun(@(x) x(1,2), {rho,p,lower,upper}); % Extract rho, 95%CI, and P-value
end

% If no result could be found for the 95% CI or the P-value, make
% sure that there result for Rho is also NaN
resArray = corrTable{:,{'Rho','2.5%CI', '97.5%CI','P-value'}};
resArray(any(isnan(resArray),2),:) = nan;
corrTable{:,{'Rho','2.5%CI', '97.5%CI','P-value'}} = resArray;

% [newrowOrder, newcolOrder] = cellfun(@(x) clustRowsAndCols(x), {RHO, RHOsig},'UniformOutput',false);
% tabs = cellfun(@(x,y,z) x(y,z), {RHOTab, RHOsigTab}, newrowOrder, newcolOrder,'UniformOutput',false);
% [RHOTab, RHOsigTab] = deal(tabs{:});
% pValTab = pValTab(newrowOrder{1}, newcolOrder{1});
