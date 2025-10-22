function [metadataPath, metadata] = pruneFluxOutliersFromMetadataADRC(sampleImportanceTable,metadataPath,top)
% Prune outliers from metadata file

% Get the top 3 samples to remove
samplesToRemove = sampleImportanceTable.ID(1:top);

% Load metadata
metadata = readtable(metadataPath,'VariableNamingRule','preserve');

% Remove samples
metadata(matches(metadata.ID,samplesToRemove),:) = [];

% Save updated metadata
writetable(metadata,metadataPath);
end