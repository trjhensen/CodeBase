function metadataMars = appendMappingStatsMetadata(metadata, outputPathMars)
% inputs:

% 1) Total read count before mapping
rawReadPath = fullfile(outputPathMars,'preprocessedInput_afterRenaming.csv'); % Raw reads
% 2) Species level read counts before and after mapping
specReadCountPath = fullfile(outputPathMars,'metrics','Species','read_counts.csv'); % Species reads before and after mapping
% 3) pielous_evenness
specPielouStats = fullfile(outputPathMars,'metrics','Species','pielous_evenness.csv'); % Pielou evenness before and after mapping
% 4) firmicutes_bacteroidetes_ratio
fbRatioPath = fullfile(outputPathMars,'metrics','Phylum','firmicutes_bacteroidetes_ratio.csv'); % 

% Add the total read count to the metadata file
totalRC = readtable(rawReadPath,'VariableNamingRule','preserve','Range','B:ZZZ'); % Exclude the taxon names
totalRC = varfun(@sum,totalRC);
totalRC = rows2vars(totalRC);
totalRC = renamevars(totalRC,totalRC.Properties.VariableNames,{'ID','total_reads'});
totalRC.ID = erase(totalRC.ID,'sum_');

% Merge with metadata
metadata = outerjoin(metadata,totalRC,'Keys','ID','MergeKeys',true,'Type','left');

% Add the number of species reads before and after mapping
specRC = readtable(specReadCountPath,'VariableNamingRule','preserve');%,'Range','B:ZZZ'); % Exclude the taxon names
specRC = rows2vars(specRC,'VariableNamesSource','Var1','VariableNamingRule','preserve');
specRC = renamevars(specRC,specRC.Properties.VariableNames,{'ID','total_species_reads','mapped_species_reads'});
specRC.('read_coverage') =  specRC.("mapped_species_reads") ./ specRC.("total_species_reads");

% Merge with metadata
metadata = outerjoin(metadata,specRC,'Keys','ID','MergeKeys',true,'Type','left');

% Add pieloe evenness to metadata
specPie = readtable(specPielouStats,'VariableNamingRule','preserve');%,'Range','B:ZZZ'); % Exclude the taxon names
specPie = rows2vars(specPie,'VariableNamesSource','Var1','VariableNamingRule','preserve');
specPie = renamevars(specPie,specPie.Properties.VariableNames,{'ID','total_pieloe_even','mapped_pieloe_even'});

% Merge with metadata
metadata = outerjoin(metadata,specPie,'Keys','ID','MergeKeys',true,'Type','left');

% Add bacteroides/firmicutes ratio
fbRatio = readtable(fbRatioPath,'VariableNamingRule','preserve');%,'Range','B:ZZZ'); % Exclude the taxon names
fbRatio = rows2vars(fbRatio,'VariableNamesSource','Var1','VariableNamingRule','preserve');
fbRatio = renamevars(fbRatio,fbRatio.Properties.VariableNames,{'ID','total_firmicutes_bacteroidetes_ratio','mapped_firmicutes_bacteroidetes_ratio'});

% Merge with metadata
metadata = outerjoin(metadata,fbRatio,'Keys','ID','MergeKeys',true,'Type','left');

metadataMars = metadata;
end