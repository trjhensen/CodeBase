function diversityStats = mappingSummaryStats(paths)
% Aim: Calculate the species richness, read counts, and pielou evennes
% index. 
% input:
% metadata
%
% output: Table with the species richness, shannon, pieloe, and read counts
% diversityStats.

% inputs:
mappingPaths.preMappedPath = fullfile(paths.outputs,'resultMARS', 'normalized_preMapped','normalized_preMapped_species.csv');
mappingPaths.mappedPath = fullfile(paths.outputs,'fluxes', 'analysis','WBM_relative_abundances.csv');

% Load the pre-mapped, mapped present, and mapped absent abundance data
microbiomes = structfun(@(x) readtable(x),mappingPaths,'UniformOutput',false);

% Prepare the WBM relative abundance data for further analysis
wbmAbundances = microbiomes.mappedPath;
wbmAbundances(:,2:end) = fillmissing(wbmAbundances(:,2:end),'constant',0);
wbmAbundances.Properties.RowNames = erase(wbmAbundances.Row,{'_female','_male','mWBM_'});
wbmAbundances.Row = [];
wbmAbundances = rows2vars(wbmAbundances);
wbmAbundances = renamevars(wbmAbundances,'OriginalVariableNames','Taxon');
microbiomes.mappedPath = wbmAbundances;

% Load the metadata
metadata = readtable(paths.metadata,'VariableNamingRule','preserve');

% Remove samples not in the metadata table
rmVarFun = @(x) removevars(x, setdiff( x.Properties.VariableNames(2:end)', metadata.ID)' );
microbiomes = structfun(rmVarFun,microbiomes,'UniformOutput',false);

% Remove taxa with relative abundances above 1e-6 in zero samples
rmTaxaFun = @(x) x( any( x{:,2:end} > 1e-6, 2) , :);
microbiomes = structfun(rmTaxaFun,microbiomes,'UniformOutput',false);

% Calculate diversity metrics for each sample and add information to
% metadata variable
diversityData = structfun(@getDiversityMetrics,microbiomes,'UniformOutput',false);


% Now, calculate the fold changes from the pre-mapped to the mapped metrics

% Get the variable names
varNames = string(diversityData.mappedPath.Properties.VariableNames);
foldChangeFun = @(x) diversityData.mappedPath.(x) ./ diversityData.preMappedPath.(x);
foldChanges = arrayfun(foldChangeFun, varNames,'UniformOutput',false);
diversityData.foldChanges = array2table(horzcat(foldChanges{:}),'VariableNames',varNames);

% Then, calculate the mean averages for each variable in each field
formNum = @(x) string( round(x ,3) ); % Round numbers 
funCont = @(x) append( formNum(mean(x,'omitmissing')), " (" , formNum( std(x,'omitmissing') ), ")"); % Calculate mean and SD
diversityStats = structfun( @(x) funCont(table2array(x))', diversityData, 'UniformOutput', false); % Run functions on each dataset
diversityStats = struct2table(diversityStats,'RowNames',varNames); % Convert structured array to table

% Calculate p-values for the differences before and after mapping
ttestFun = @(x) ttest2( diversityData.preMappedPath.(x), diversityData.mappedPath.(x) ); % Two-sample t-test. Differences in means?
[~,pVals,~,~] = arrayfun(ttestFun, varNames, 'UniformOutput', false);
pVals = cell2mat(pVals)';

% Add p-values to statistics
diversityStats = addvars(diversityStats, pVals,'NewVariableNames','P-value');

% Add the total number of species to the table
numSpecies = [height(microbiomes.preMappedPath) height(microbiomes.mappedPath) height(microbiomes.mappedPath)/height(microbiomes.preMappedPath) nan];
diversityStats{'Total species',:} = numSpecies;


% Calculate the difference in read counts
[~,pvalReadCounts,~,~] = ttest2(metadata.total_reads, metadata.mapped_species_reads );

% Add the differences in read counts and the associated fold change
numReads = [funCont(metadata.total_reads) funCont(metadata.mapped_species_reads) funCont(metadata.mapped_species_reads./metadata.total_reads), pvalReadCounts];
diversityStats{'Read count',:} = numReads;


end

function mappingMetadata = getDiversityMetrics(microbiomeData)

microbiomeMatrix = table2array(microbiomeData(:,2:end))';
microbiomeMatrix(microbiomeMatrix==0) = nan;

mappingMetadata = table();
mappingMetadata.SpeciesRichness = sum(~isnan(microbiomeMatrix),2); % Calculate species richness
mappingMetadata.SimpsonIndex = 1-sum(microbiomeMatrix.^2,2,'omitmissing'); % Calculate simpson index
mappingMetadata.shannonIndex = -sum(microbiomeMatrix .* log(microbiomeMatrix),2,'omitmissing'); % Calculate shannon index
mappingMetadata.pielouIndex = mappingMetadata.shannonIndex ./ log(mappingMetadata.SpeciesRichness); % Calculate the pieloe evennes index
end











