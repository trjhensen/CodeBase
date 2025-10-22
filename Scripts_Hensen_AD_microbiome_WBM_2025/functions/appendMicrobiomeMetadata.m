function mergedMetadata = appendMicrobiomeMetadata(metadata, microbiomeInputFolder)
% Add microbiome metadata to patient metadata file

% input
% metadata
% microbiomeInputFolder
% outputPath


% Prepare functions for efficient reading of microbiome mapping tables
pathFun = @(x) fullfile(microbiomeInputFolder,'unprocessed','mapping_files',strcat(x,'_mapping_file.txt')); % Generate path
tabFun = @(x) readtable(pathFun(x),'ReadVariableNames',false,'VariableNamingRule','preserve','Delimiter','\t'); % Read table
renameFun = @(x,y) renamevars(x, x.Properties.VariableNames,readcell(pathFun(y),'Range', '1:1')); % Add variable names to table

% Load microbiome sample tables and add variable names
files = {'194890','194939'}; % File names to load
warning('off'); mappingFiles = cellfun(tabFun,files,'UniformOutput',false); warning('on') % Read both tables
mappingFiles = cellfun(renameFun, mappingFiles,files, 'UniformOutput',false);

% Combine mapping files into a single file
microbiotaMapping = outerjoin(mappingFiles{1},mappingFiles{2},'MergeKeys',true); 
microbiotaMapping = renamevars(microbiotaMapping,'sample_name','ID');

% Merge microbiome metadata with patient metadata
metadata = outerjoin(metadata,microbiotaMapping,'Keys','ID','MergeKeys',true,'Type','left');

% Load additional microbiome metadata information and map the microbiome
% metadata on the full metadata file.
microbiomeMetadataPath = fullfile(microbiomeInputFolder,'ADRC_microbiome_metadata.xlsx');
opts = detectImportOptions(microbiomeMetadataPath); 
opts = setvartype(opts, opts.VariableNames(1), 'string'); 
opts.VariableNamingRule = 'preserve';
microbiomeMD = readtable(microbiomeMetadataPath,opts);

% Map microbiome metadata table to metadata IDs and merge tables
microbiomeMD = renamevars(microbiomeMD,"sample_name",'ID');
microbiomeMD.ID = append('X',erase(microbiomeMD.ID,'15448.')); % Change IDs for mapping onto metadata
mergedMetadata = outerjoin(metadata,microbiomeMD,'Keys','ID','MergeKeys',true,'Type','left');

% Convert date variables from cell strings to date times
mergedMetadata = convertvars(mergedMetadata,{'collection_date','accession_date'},@(x) datetime(x,'InputFormat', 'M/dd/yy'));

end
