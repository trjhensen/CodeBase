function [mappingPotential,mappingPotentialFilt, savePath] = findPotentialReadCovIncr(outputPathMars, microbiotaPathFilt)
% Description: Find for each unmapped microbial species how much the read coverage would
% improve if included.

% Load summary statistics for total read counts before and after mapping
rcStatsPath = fullfile(outputPathMars,'metrics','Species','read_counts.csv');
readCountsSum = readtable(rcStatsPath);

% Load summary statistics for unmapped taxa
unmappedSpeciesStatsPath = fullfile(outputPathMars,'metrics','Species','unmapped_abundanceMetrics_Species.csv');
unmappedMetrics = readtable(unmappedSpeciesStatsPath);

% Load pre-mapped microbiome data
microbiome = readtable(fullfile(outputPathMars,'preprocessedInput_afterRenaming.csv'));

% Filter on analysed samples and microbes
microbiomeFilt = readtable(microbiotaPathFilt,'VariableNamingRule','preserve','ReadRowNames',true);
microbiomeFilt = renamevars(rows2vars(microbiomeFilt),'OriginalVariableNames','Taxon');
microbiomeFilt = convertvars(microbiomeFilt, 'Taxon', @(x) replace(x,"_"," "));
microbiome = microbiome(contains(microbiome.Taxon, microbiomeFilt.Taxon),:);
% microbiome = microbiome(:, microbiomeFilt.Properties.VariableNames);

% Calculate current read coverage statistics from the pre-mapped relative abundances
readCountTable = rows2vars(readCountsSum,"VariableNamesSource","Var1",'VariableNamingRule','preserve');
readCountTable.("Read coverage") = readCountTable.("Post mapping")./readCountTable.("Pre mapping");

% Preprocess unmapped taxa names
unmappedMetrics.Taxon = strrep(unmappedMetrics.Taxon, '_', ' ');
unmappedMetrics = sortrows(unmappedMetrics, 'mean', 'descend');

% Also filter the readCountTable and unmappedMetrics table
% readCountTable = readCountTable( matches(readCountTable.OriginalVariableNames, microbiomeFilt.Properties.VariableNames'),:);
unmappedMetrics = unmappedMetrics( matches(unmappedMetrics.Taxon, microbiomeFilt.Taxon) ,:);

% Calculate the mapping coverage improvements after mapping of each
% unmapped microbial species. 
rcountFun = @(x,y) sum(x{ contains(x.Taxon,y), 2:end },1); % Extract the total read counts of the unmappedMetrics taxon
rcovFun = @(x,y,r) sum([r.("Post mapping")'; rcountFun(x,y)], 1) ./ r.("Pre mapping")'; % Calculate the potential improved read coverage
summarise = @(a,b) [mean(b), mean(a), mean(a)-mean(b)]; % Calculate summary statistics
rcovIncr = @(x,y,r) summarise(rcovFun(x,y,r)', r.("Read coverage")); % Calculate the increase in average read coverage if mapped

% Perform calculations
potReadCovSum = arrayfun(@(i) rcovIncr(microbiome, unmappedMetrics.Taxon(i), readCountTable),1:height(unmappedMetrics),'UniformOutput',false);
potReadCovSum = vertcat(potReadCovSum{:}); % Extract summary statistics
potReadCovSumTable = array2table(potReadCovSum,'VariableNames',{'curr_readCov','potential_readCov','potential_Incr'}); % Create table
mappingPotential = [unmappedMetrics potReadCovSumTable]; % Append new statistics to R.A of unmapped microbial species.

if 1 
    % Remove all unmapped taxa with sp
    mappingPotentialFilt = mappingPotential(~contains(mappingPotential.Taxon,' sp'),:);
    mappingPotentialFilt.cumIncrease = mappingPotentialFilt.curr_readCov + cumsum(mappingPotentialFilt.potential_Incr);

end

% Calculate the cumulative mean increase in read coverages for the unmappped microbial species
mappingPotential.cumIncrease = mappingPotential.curr_readCov + cumsum(mappingPotential.potential_Incr);

savePath = fullfile(outputPathMars,'metrics','Species','potential_higher_readCov_when_mapping_more_microbes.xlsx');
cellfun(@(x,y) writetable(x,savePath, 'Sheet', y), {mappingPotential, mappingPotentialFilt}, {'all_taxa','characterised_taxa'})
writetable(mappingPotential,savePath)

% Create overview figure on the effect of improved mapping on the data
% plotAbsentTaxaEffectOnMARScoverage(microbiome, unmappedMetrics, readCounts, outputPathMars);
end
