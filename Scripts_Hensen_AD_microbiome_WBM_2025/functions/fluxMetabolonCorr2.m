function [RHOTab, RHOsigTab, pValTab] = fluxMetabolonCorr2(fluxPath,metabolonPath)
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

% Correlate flux values with metabolomics
[RHO, PVAL] = corr(fluxes{:,:}, metabolomics{:,:}, 'type', 'Pearson','Rows','pairwise');

% Filter RHO and p<0.05
RHOsig = RHO;
RHOsig(PVAL>0.05) = 0;

tabFun = @(x) array2table(x,"RowNames",fluxes.Properties.VariableNames, "VariableNames",metabolomics.Properties.VariableNames);
tabsRaw = cellfun(tabFun, {RHO, RHOsig, PVAL},'UniformOutput',false);
[RHOTab, RHOsigTab, pValTab] = deal(tabsRaw{:});

% [newrowOrder, newcolOrder] = cellfun(@(x) clustRowsAndCols(x), {RHO, RHOsig},'UniformOutput',false);
% tabs = cellfun(@(x,y,z) x(y,z), {RHOTab, RHOsigTab}, newrowOrder, newcolOrder,'UniformOutput',false);
% [RHOTab, RHOsigTab] = deal(tabs{:});
% pValTab = pValTab(newrowOrder{1}, newcolOrder{1});
