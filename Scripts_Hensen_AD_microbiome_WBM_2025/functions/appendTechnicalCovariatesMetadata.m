function metadataTechCov = appendTechnicalCovariatesMetadata(metadata, microbiomeInputFolder)
% Append technical covariate information to metadata file

techCovPath = fullfile(microbiomeInputFolder,'fecal_quality_filter_copy.txt');

% Load technical covariates
opts = detectImportOptions(techCovPath); 
opts = setvartype(opts, opts.VariableNames(1), 'string'); 
opts.VariableNamingRule = 'preserve';
techCov = readtable(techCovPath,opts);

% Map table to metadata IDs
techCov = renamevars(techCov,"#SampleID",'ID');
techCov.ID = append('X',erase(techCov.ID,'15448.')); % Change IDs for mapping onto metadata

% Remove the leading zeros from techCov. 
techCovIDs = erase(techCov.ID,'X');
techCovIDs = regexprep(techCovIDs,'^0+','');
techCov.ID = append('X',techCovIDs);

% Merge technical covariate information with metadata
% !! There are 52 metadata samples without information on technical covariates 
% samplesWithMissingData = setdiff(metadata.ID,techCov.ID); 
metadata = outerjoin(metadata,techCov,'Keys','ID','MergeKeys',true,'Type','left');

% Change the name of the ethanol addition
metadata = renamevars(metadata,'Ethanol added (Y/N)','Ethanol_added');


metadataTechCov = metadata;

% % Save updated metadata file with technical co-variate information
% fullMetadataPath = fullfile(outputPath,'processedMetadata.csv');
% writetable(metadataTechCov,fullMetadataPath)
end