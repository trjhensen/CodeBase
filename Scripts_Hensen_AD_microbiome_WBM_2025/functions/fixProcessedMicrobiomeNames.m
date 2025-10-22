function [marsFixedInputTablePath, microbiomes] = fixProcessedMicrobiomeNames(processedMicrobes)
% input:
% Path to processed microbiome read table folder

% This function fixes the gut microbiota sample names so that every sample
% can be mapped to their corresponding metadata.

% Get path to MARS input data
marsInputTablePath = fullfile(processedMicrobes,'MARS_input_combined.csv');

% Load data and extract sample names
microbiomes = readtable(marsInputTablePath,'VariableNamingRule','preserve');
sampleNamesToFix = microbiomes.Properties.VariableNames(2:end); % Do not include the taxa

% Fix names
fixedNames = fixADRCsampleNames(sampleNamesToFix);

% Add fixed names to microbiome samples and save updated table
microbiomes.Properties.VariableNames(2:end) = fixedNames;

marsFixedInputTablePath = fullfile(processedMicrobes,'MARS_input_combined_fixed_names.csv');
writetable(microbiomes,marsFixedInputTablePath);

end