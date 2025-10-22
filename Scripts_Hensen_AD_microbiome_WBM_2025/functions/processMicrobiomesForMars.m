function [combined, filePath] = processMicrobiomesForMars(unprocessedMicrobiomePath, processedMicrobiomePath)

% Force delete all current outputs in the processedMicrobiomePath
if isfolder(processedMicrobiomePath); rmdir(processedMicrobiomePath,'s'); end

% and create a new clean folder for the processed microbiome outputs.
mkdir(processedMicrobiomePath)

files = {'194890','194939'};
microbiomeData = cell(1,2);
for i=1:length(files)
    disp(strcat('Iteration=',string(i)))
    disp('Load decontaminated microbiome data')
    % Load decontaminated read counts
    microbiome = readtable([unprocessedMicrobiomePath filesep strcat('reads_',files{i},'_decontaminated.txt')],'VariableNamingRule','preserve');
    % Transpose data
    taxa = microbiome.Properties.VariableDescriptions(2:end);
    microbiome = rows2vars(microbiome,"VariableNamesSource",'ID');
    microbiome.OriginalVariableNames=taxa';
    microbiome = renamevars(microbiome,"OriginalVariableNames","Taxon");
    % Replace domain with kingdom
    microbiome.Taxon = replace(microbiome.Taxon,'d__','k__');
    % Remove all reads not at the species level
    microbiome = microbiome(contains(microbiome.Taxon,';s__'),:);
    % Create tall arrays for dataset concatination
    microbiomeStacked = stack(microbiome,microbiome.Properties.VariableNames(2:end),'NewDataVariableName','value','IndexVariableName','sample');
    % Add microbiome to cell
    microbiomeData{i} = microbiomeStacked;
end

% Combine arrays
combined_stacked = outerjoin(microbiomeData{1},microbiomeData{2},'MergeKeys',true);

% Revert to wide array and replace NaN values with zeros
combined = unstack(combined_stacked,'value','sample'); 
combined(:,2:end) = fillmissing(combined(:,2:end),'constant',0);

% Extract species names
taxa = combined.Taxon;
species = regexp(taxa,'(?<=s__).*','match'); % Extract species names
species = regexprep(string(species),'_[A-Z]',''); % Perform name transformation

% manual taxa repair using NCBI taxonomy
species(matches(species,'Ruminiclostridium siraeum')) = "Eubacterium siraeum";

% Replace with fixed species names
taxa = regexprep(taxa,'(?<=s__).*',''); % Remove old species names
taxa = append(taxa,species); % Add fixed names

% Add back updated taxa names
combined.Taxon = taxa;

% Save
filePath = [processedMicrobiomePath filesep strcat('MARS_input_combined.csv')];
writetable(combined,filePath)

end

% filePaths = fullfile(unprocessedMicrobiomePath,strcat('reads_',files,'_decontaminated.txt'));
% microbiomeData = cellfun(@loadMicrobiome, filePaths,'UniformOutput',false);
% combined = outerjoin(microbiomeData{1},microbiomeData{2},'Keys','Taxon','MergeKeys',true);
% combined = fillmissing(combined,'constant',0,'DataVariables',vartype('numeric'));
% 
% % Process microbe names
% tax = combined.Taxon; 
% tax = string(regexp(tax,'(?<=s__).*','match')); % Extract species names
% tax = regexprep(tax,'_[A-Z]',''); % Perform name transformation
% tax = replace(tax,"Ruminiclostridium siraeum","Eubacterium siraeum"); % manual taxa repair using NCBI taxonomy
% combined.Taxon = append(regexprep(combined.Taxon,'(?<=s__).*',''),tax); % Replace original species names with fixed species names
% 
% function tab = loadMicrobiome(p)
% mName = 'Taxon'; % Set name of column with microbes
% tab = readtable(p,'ReadVariableNames',false,'ReadRowNames',true); % Read table
% tab = rows2vars(tab,'VariableNamingRule','preserve'); % Set samples to columns
% tab = removevars(tab,'OriginalVariableNames'); % Remove empty column
% tab = addvars(tab,string(readcell(p,'Range', 'B1:ZZZZ1'))','NewVariableNames',mName,'Before',1); % Add taxa as first column
% tab.(mName) = replace(tab.(mName),'d__','k__'); % Replace domain with kingdom
% tab = tab(contains(tab.(mName),';s__'),:); % Remove all reads not at the species level
% end