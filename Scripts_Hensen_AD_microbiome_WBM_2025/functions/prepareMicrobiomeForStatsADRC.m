function preparedMicrobiome = prepareMicrobiomeForStatsADRC(rxnsOfInterest, fmCorrPath, mappedMicrobePath)
% Aim: This function finds the significant flux associations with clinical
% metadata. Then, the associated microbes are found in the flux-microbe
% correlation table (fmCorrPath). Once the associated microbial species are
% found, we filter the microbial relative abundances on the microbes of
% interest. The output will be the prepared metadata. The relative
% abundances  will also be pruned on samples in the metadata.

% input:
% rxnNames
% pvals: Reaction associated p-values
% fmCorrPath: Path to flux-microbe correlations
% mappedMicrobePath: Processed relative abundances

% Find the associated microbes for each reaction
fmCorr = readtable(fmCorrPath,"ReadRowNames",true,'VariableNamingRule','preserve');
fmCorr = fmCorr(:,rxnsOfInterest');
microbesToTest = fmCorr;
microbesToTest{:,:} = ~isnan(microbesToTest{:,:});

% Find union of all microbes that associate with the found metabolites
microbesToTest = microbesToTest.Properties.RowNames(any(table2array(microbesToTest),2))';

% Load gut microbiome data
microbiome = readtable(mappedMicrobePath,'VariableNamingRule','preserve','ReadRowNames',true);

% Filter on microbes to test
microbiome = microbiome(:,microbesToTest);

% Add back ID column
microbiome = addvars(microbiome,microbiome.Properties.RowNames,'Before',1,'NewVariableNames','ID');
microbiome.Row = [];

% Log2 transform relative abundances and normalise them
inputMatrix = table2array(microbiome(:,2:end));
log2Data = log2(inputMatrix);
log2Data(isinf(log2Data)) = nan;
microbiome{:,2:end}=normalize(log2Data);

preparedMicrobiome = microbiome;
end