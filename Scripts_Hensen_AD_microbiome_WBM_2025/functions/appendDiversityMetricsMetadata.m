function metadataDiversities = appendDiversityMetricsMetadata(metadata,marsFolder)
% This function adds microbiome diversity metrics to the metadata variable

microbiomePath = fullfile(marsFolder,'renormalized_mapped_forModelling','renormalized_mapped_forModelling_species.csv');

% Lets test if I can test the fluxes against Dementia status
microbiomeData = readtable(microbiomePath,'VariableNamingRule','preserve');
microbiomeData = rows2vars(microbiomeData,'VariableNamingRule','preserve');
microbiomeData.Properties.VariableNames = microbiomeData{1,:};
microbiomeData(1,:) = [];
microbiomeData = renamevars(microbiomeData,'Taxon','ID');
% microbiomeData.ID = erase(microbiomeData.ID,{'mWBM_','_female','_male'});

% Make sure that the inputTable and metadata have the same samples
[~,idxa,idxb] = intersect(string(metadata.ID),string(microbiomeData.ID),'stable');

if numel(idxb)<length(metadata.ID)
    error('COBRA:BadInput', 'No overlapping samples could be found between the reads/flux table and the metadata table.')
else
    metadata = metadata(idxa,:);
    microbiomeData = microbiomeData(idxb,:);
end

% Calculate diversity metrics for each sample and add information to
% metadata variable
microbiomeMatrix = cell2mat(table2array(microbiomeData(:,2:end)));
microbiomeMatrix(microbiomeMatrix==0) = nan;
metadata.SpeciesRichness = sum(~isnan(microbiomeMatrix),2); % Calculate species richness
metadata.SimpsonIndex = 1-sum(microbiomeMatrix.^2,2,'omitmissing'); % Calculate simpson index
metadata.shannonIndex = -sum(microbiomeMatrix .* log(microbiomeMatrix),2,'omitmissing'); % Calculate shannon index
metadata.pielouIndex = metadata.shannonIndex ./ log(metadata.SpeciesRichness); % Calculate the pieloe evennes index

% Set function output
metadataDiversities = metadata;
end




