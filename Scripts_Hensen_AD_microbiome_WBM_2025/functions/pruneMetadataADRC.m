function [metadataPruned, prunedMetadataPath, metadataIntermedPruning] = pruneMetadataADRC(metadata, metadataOutputFolder)


% Now, we remove individuals with impaired cognitive health, but
% with no cognitive decline. 
metadata(matches(string(metadata.NACCUDSD),'impaired non-mci'),:) = [];

groupsummary(metadata,'NACCUDSD')

metadataIntermedPruning{1} = metadata;

% figure; heatmap(metadata,'AD','NACCUDSD')
% Finally, we will remove all non-AD MCI and Dementia patients. This step
% is performed to remove noise from non-AD cognitive decline samples. 
cogData = string(metadata.NACCUDSD);
isMCIorDementia = matches(cogData, {'MCI', 'AD dementia'});
noAD = matches(metadata.AD, 'no_AD');
metadata(isMCIorDementia & noAD, :) = [];

% Next, we will remove two samples which were collected much earlier than
% all other samples (June 2021). Note that there is also an argument to be made to not
% remove these samples.
if 1; metadata(matches(metadata.ID,{'X5827089', 'X5946764'}),:) = []; end


groupsummary(metadata,'NACCUDSD')

% Save intermediate data variable
metadataIntermedPruning{2} = metadata; % figure; heatmap(metadata,'AD','NACCUDSD')

% There is one sample where only 4.8% of the reads passed the quality
% filter. This sample will be removed.
metadata(metadata.fraction_passing_quality_filter<0.05,:) = []; 

% Next, there are 4 samples with less than 1 million mapped reads. Three of
% which also had less than 1 million total reads. Lets remove these
% samples.
metadata(metadata.mapped_species_reads<1000000,:) = []; 

% There are also two samples with both very low read coverages (<25%) and
% an extreme firmicutes/bacteroidetes ratio in the raw reads (>1000). These
% samples are also removed.
metadata(metadata.read_coverage<0.25,:) = [];

groupsummary(metadata,'NACCUDSD')
metadataPruned = metadata;

% Save updated metadata file with technical co-variate information
prunedMetadataPath = fullfile(metadataOutputFolder,'prunedProcessedMetadata.csv');
writetable(metadataPruned,prunedMetadataPath)
end