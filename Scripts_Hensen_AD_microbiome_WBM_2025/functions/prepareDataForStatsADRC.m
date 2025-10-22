function [preparedInputTable, preparedMetadata, preparedRawInputTable] = prepareDataForStatsADRC(dataPath, metadataPath,transformData)

if nargin<3
    transformData = true;
end

% Lets test if I can test the fluxes against Dementia status
inputTable = readtable(dataPath,'VariableNamingRule','preserve');

if ~isvar(inputTable,'ID') && isvar(inputTable,'Row')
    inputTable = renamevars(inputTable,'Row','ID');
end

inputTable.ID = erase(inputTable.ID,{'mWBM_','_female','_male'});

% Remove sex information if present
if any(matches(inputTable.Properties.VariableNames,'sex','IgnoreCase',true))
    inputTable.Sex = [];
end

metadata = readtable(metadataPath,'VariableNamingRule','preserve');
clc

if 0
    if 0
        metadata(matches(metadata.sex,'female'),:) = [];
    else
        metadata(matches(metadata.sex,'male'),:) = [];
    end
end

if 0
    if 0
        metadata(~matches(metadata.APOE_E4,'E4'),:)=[];
    else
        metadata(matches(metadata.APOE_E4,'E4'),:)=[];
    end
end

% Make sure that the inputTable and metadata have the same samples
[~,idxa,idxb] = intersect(string(metadata.ID),string(inputTable.ID),'stable');

if numel(idxb)<length(metadata.ID)
    error('COBRA:BadInput', 'No overlapping samples could be found between the reads/flux table and the metadata table.')
else
    metadata = metadata(idxa,:);
    inputTable = inputTable(idxb,:);
end

preparedRawInputTable = inputTable;

% Transform the inputData by performing a log2 transformation and
% z-scaling.
if transformData == true
    inputMatrix = table2array(inputTable(:,2:end));
    log2Data = log2(inputMatrix);
    log2Data(isinf(log2Data)) = nan;
    inputTable{:,2:end}=normalize(inputMatrix);
end

% Set output
preparedInputTable = inputTable;


% Transform the metadata, age, bmi, and total  sequence count for 
% statistical analysis.
metadata.age_at_collection = normalize(metadata.age_at_collection);
metadata.NACCBMI = normalize(log10(metadata.NACCBMI));
metadata.mapped_species_reads = normalize(log10(metadata.mapped_species_reads));

% Process AD data
metadata.AD = categorical(metadata.AD,{'AD','no_AD'});


% Think about the following:
dchacSampsToRm = {'X42510869','X42754461','X42547915','X4262232'};
preparedInputTable.("DM_dchac[bc]")(matches(preparedInputTable.ID,dchacSampsToRm)) = nan;

preparedMetadata = metadata;


%preparedMetadata(matches(preparedMetadata.APOE_ALLELE,""),:)=[];

end