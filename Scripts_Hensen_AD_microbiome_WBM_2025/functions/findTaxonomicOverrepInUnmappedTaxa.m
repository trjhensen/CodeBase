function [taxonEnrichmentTables,groupcAllLevels] = findTaxonomicOverrepInUnmappedTaxa(microbiotaPath,taxonomyPath,microbiotaWbmPath)
% Goal: To identify biases in the mapping procedure. 
% Which Phyla and genera were underrepresented in the mapped taxa?
% Identify the phyla, clades, orders, families, and genera had the highest loss in speces 


% Raw microbiota relative abundances with processed names
microbiota = readtable(microbiotaPath,'VariableNamingRule','preserve','ReadRowNames',true);
microbiota = array2table( microbiota{:,:}' , 'RowNames',microbiota.Properties.VariableNames, 'VariableNames',microbiota.Properties.RowNames');
microbiota = addvars(microbiota, microbiota.Properties.RowNames, 'NewVariableNames','Species','Before',1);
microbiota.Properties.RowNames = {};
microbiota.Species = replace(microbiota.Species,'_',' ');

% Repair names in the mapped taxa to ensure overlap
microbiota.Species(matches(microbiota.Species,'Ruminococcus chamellensis')) = {'Ruminococcus champanellensis'};
microbiota = removevars(microbiota, microbiota.Properties.VariableNames(2:end));

% Species-phylum mapping after species renaming
% Load Taxonomy information
opts = detectImportOptions(taxonomyPath);
opts.SelectedVariableNames = "Taxon";
taxonomyInfo = readmatrix(taxonomyPath, opts);
taxonomyInfo = append(taxonomyInfo,';');

% Information on which microbial species were mapped on AGOR2+APOLLO

% Split taxonomic data into multiple columns
levels = {'Kingdom','Phylum','Clade','Order','Family','Genus','Species'};
levelAbbreviations = ['k','p','c','o','f','g','s'];
expresFun = @(x) ['(?<=',x,'__)(.*?)(?=\;)']; % Get the taxonomy information of interest
regexFun = @(x) string(regexp(taxonomyInfo,expresFun(x),'match')); % Extract matches in regex
taxaInfo = arrayfun(regexFun, levelAbbreviations,'UniformOutput',false); % Run regex for all levels
taxaInfo = array2table(cellstr(horzcat(taxaInfo{:})),'VariableNames',levels); % Generate table

% Select phylum infotaxaInfo
% taxaInfo = taxaInfo(:,{'Species','Phylum'});

% Map species names on taxaInfo. 
% microbiota = microbiota(:,{'Species'});

mergedTable = outerjoin(microbiota,taxaInfo,"Keys","Species","Type","left",'MergeKeys',true);

% Remove duplicate rows (WHY DO THEY EXIST?)
[~, ~, ind] = unique(mergedTable.Species,'stable'); % Find the indices of all unique species
mergedTable = mergedTable(unique(ind),:); % Filter on the unique indices of all species

% Identify which microbial species in the raw data were mapped.
% Load processed and mapped species relative abundances 
microbiotaWBM = readtable(microbiotaWbmPath,'ReadRowNames',true,'VariableNamingRule','preserve');
% microbiotaWBM.('Sum of taxa') = [];

% Process microbiota data
microbiotaWBM = renamevars(rows2vars(microbiotaWBM),'OriginalVariableNames','Species');
microbiotaWBM.Species = replace(microbiotaWBM.Species,'_',' ');

% Repair names in the mapped taxa to ensure overlap
microbiotaWBM.Species(matches(microbiotaWBM.Species,'Ruminococcus chamellensis')) = {'Ruminococcus champanellensis'};
% microbiotaWBM.Species(matches(microbiotaWBM.Species,'Companilactobacillus farciminis')) = {'Comilactobacillus farciminis'};

% Add information on which microbial species were mapped onto the WBMs
mergedTable.mapped = matches(mergedTable.Species,microbiotaWBM.Species);
sum(mergedTable.mapped)

% Calculate the total loss in species
totalSpec = height(mergedTable);
totalLoss = sum(mergedTable.mapped==0);

% Convert binary mapped to mapped/unmapped encoding
mergedTable = convertvars(mergedTable, 'mapped', @(x) renamecats(categorical(x), {'Unmapped','Mapped'}));

% Get all unique classifications
classNames = setdiff(mergedTable.Properties.VariableNames,{'Kingdom','mapped','Species'},'stable');
% taxonTables = cellfun(@(x) unique(taxonInfoNew.(x)), classNames,'UniformOutput',false);

% Get for each phylum the number of mapped and unmapped species
grCountFun = @(x) renamevars(removevars(groupcounts(mergedTable,{x,'mapped'}),'Percent'),x,'Taxon');
procGrCountFun = @(x) convertvars(unstack(x,'GroupCount','mapped'),{'Mapped','Unmapped'}, @(x) fillmissing(x,'constant',0));

groupcAllLevels = cellfun(@(x) procGrCountFun(grCountFun(x)), classNames,'UniformOutput',false);

% Create enrichment tables for each taxonomic level
enrichmentTables = cellfun(@(x) findMappingLossEnrichment(x, totalLoss, totalSpec), groupcAllLevels,'UniformOutput',false);
taxonEnrichmentTables = cellfun(@(x,y) renamevars(y,'Taxon',x), classNames,enrichmentTables,'UniformOutput',false);

end

function enrichmentTable = findMappingLossEnrichment(grCountsLevel, totalLoss, totalSpec)
% Define variable names for enrichment table
varNames = {'Number_of_species','Reduction_in_taxa','Loss_in_Taxon','Loss_outside_Taxon','Mapped_in_Taxon','Mapped_outside_Taxon','OddsRatio','2.5%CI','97.5%CI','pValue'};
enrichmentTable = array2table( nan(length(height(grCountsLevel)),length(varNames) ),'VariableNames',varNames);

for i=1:height(grCountsLevel)
    lossInGroup = grCountsLevel.Unmapped(i);
    leftInGroup = grCountsLevel.Mapped(i);
    lossOutGroup = totalLoss - lossInGroup;
    leftOutGroup = (totalSpec - leftInGroup) - totalLoss;
    numSpecies = grCountsLevel.Unmapped(i) + grCountsLevel.Mapped(i);
    lossPerc = lossInGroup / numSpecies;
    
    crossTabulation = [lossInGroup, leftInGroup; lossOutGroup, leftOutGroup];
    [~,pVal,stats] = fishertest(crossTabulation);

    % Add results to table
    enrichmentTable{i, varNames} = [numSpecies, lossPerc, reshape(crossTabulation,1,4) , stats.OddsRatio, stats.ConfidenceInterval, pVal];
end

% Add pathway names
enrichmentTable = addvars(enrichmentTable, grCountsLevel.Taxon,'newVariableNames', 'Taxon','Before',1);
% Correct for FDR
enrichmentTable = addvars(enrichmentTable, fdrBHadjustment(enrichmentTable.pValue), 'newVariableNames', 'FDR'); 
% Sort on enriched pathways
enrichmentTable = sortrows(enrichmentTable,'pValue','ascend');

end