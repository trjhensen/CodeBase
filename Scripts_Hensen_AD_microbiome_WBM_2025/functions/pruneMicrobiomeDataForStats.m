function [microbiome, wbmMicrobiome] = pruneMicrobiomeDataForStats(microbiotaWbmPath, microbiotaPath, metadataPath)
% Goal: Process the original microbiome dataset and the mapped WBM
% microbiome dataset by filtering lowly abundant microbes and samples not
% analysed in the fluxes. 

% Load the raw microbiome dataset and transpose data
microbiome = readtable(microbiotaPath,'VariableNamingRule','preserve','ReadRowNames',true);
microbiome = array2table(microbiome{:,:}','RowNames',microbiome.Properties.VariableNames,'VariableNames',microbiome.Properties.RowNames'); % Transpose table

% Load the WBM microbiome table and process data
wbmMicrobiome = readtable(microbiotaWbmPath,'VariableNamingRule','preserve','ReadRowNames',true);
wbmMicrobiome.Properties.RowNames = erase(wbmMicrobiome.Properties.RowNames,{'mWBM_','_female','_male'});

% Load metadata table and filter samples in the microbiome tables
metadata = readtable(metadataPath,'VariableNamingRule','preserve','ReadRowNames',false);
microbiome = microbiome(metadata.ID,:);
wbmMicrobiome = wbmMicrobiome(metadata.ID,:);

% Remove taxa with relative abundancese below 1e-6
rmTaxaFun = @(x) removevars(x, find(all( x{:,:} < 1e-6, 1)));
[microbiome,wbmMicrobiome] = deal(rmTaxaFun(microbiome),rmTaxaFun(wbmMicrobiome));

end
