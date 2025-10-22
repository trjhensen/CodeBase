function [metadataPlasmaMets,savePath] = appendMetabolonToMetada(rxnsToMap, metabolonPath, metadataPath, saveFolder)
% Define inputs:
% paths.rxnsToMap

% Map plasma metabolomics onto the investigated VMH metabolites:

% Define reactions of interest and load VMH database info for the
% associated metabolites.
% For now, set all reactions with fluxes as reactions of interest
%paths.rxnsToMap = readcell(paths.fluxPath,'Range','C1:ZZ1');

% Manually map plasma metabolites to VMH ID
[vmhMets,metabolon,vmhTable3] = mapPlasmaMetsToVmhADRC(rxnsToMap,metabolonPath);

% Map metabolon info on VMH ideas
metabolonMapped = outerjoin(vmhTable3,metabolon,'Keys','CHEM_ID','MergeKeys',true,'Type','left');
metabolonMapped.Properties.RowNames = metabolonMapped.Var1;

% Load processed plasma metabolomics data
metabolomics = readtable(metabolonPath,'VariableNamingRule','preserve','Sheet','assay','ReadRowNames',true);
metabolomicsMapped = metabolomics(metabolonMapped.Properties.RowNames,:);
metabolomicsMapped.Properties.RowNames = metabolonMapped.VMHID; % Replace metabolon metabolite names with VMH IDs
metabolomicsMapped = rows2vars(metabolomicsMapped,'VariableNamingRule','preserve'); % Create tall table
metabolomicsMapped = renamevars(metabolomicsMapped,'OriginalVariableNames','PARENT_SAMPLE_NAME'); % Rename variables

% Merge metabolomics with metadata information
sampleMapping = readtable(metabolonPath,'VariableNamingRule','preserve','Sheet','colData','ReadRowNames',false);
metabolomicsMerged = outerjoin(sampleMapping,metabolomicsMapped,'Keys','PARENT_SAMPLE_NAME','MergeKeys',true,'Type','left');

% Remove samples for which no faecal samples were available
metadata = readtable(metadataPath,'VariableNamingRule','preserve');
metabolomicsMerged = metabolomicsMerged( ismember(metabolomicsMerged.NACCID,metadata.NACCID), :);

% Find individuals for which multiple plasma samples are available
duplNACCIDs = groupcounts(metabolomicsMerged,"NACCID");
multSamples = duplNACCIDs.NACCID(duplNACCIDs.GroupCount>1);

% Preallocate list for samples to be removed
sampsToRemoveList = cell(numel(multSamples),1);

for i = 1:length(sampsToRemoveList)
    % Extract participant with multiple samples
    participant = metabolomicsMerged( matches( metabolomicsMerged.NACCID,multSamples(i) ) ,:);
    % Preallocate empty cell array for samples to be removed
    samplesToRemove = {};
    
    % Define functions for removal process
    over1SampLeft = @(x) numel(setdiff(participant.PARENT_SAMPLE_NAME, x))>1; % Check if any samples are left after removal
    canSampBeRem = @(x,cat1,cat2) any(matches(x,cat1)) & any(matches(x,cat2)); % Can any sample be removed?
    appendRmSampToList = @(x,cat1) vertcat(samplesToRemove, participant.PARENT_SAMPLE_NAME(matches(x,cat1)) ); % Add sample name to list with samples to be removed
    
    % Step 1: Remove samples from individuals that did not fast before sample collection
    % Step 2: Remove samples that were taken greater than 9 months away from stool collection

    % Define variables to check
    catVarsToCheck = {...
        'PlasmaFasting','NO','Yes'; ...
        'AGMPVisitNote','Within 9 months','Greater than 9 months'...
        };

    % Identify samples to be removed
    for j=1:height(catVarsToCheck)
        if over1SampLeft(samplesToRemove) == true
            input = participant.( catVarsToCheck{j,1} );
            if canSampBeRem(input,catVarsToCheck{j,2},catVarsToCheck{j,3}) == true
                samplesToRemove = appendRmSampToList(input,catVarsToCheck{j,2}); % Add sample to be removed to list
            end
        end
    end

    % Step 3: Remove samples from individuals that were farther away in time from stool collection
    % Step 4: Remove the last samples that were collected 
    numVarsToCheck = {'ClosestAGMPVisit'; 'plasmaNACCVNUM'};

    % Define functions for removal process
    canSampBeRemNum = @(x) numel(unique(x))>1;
    appendRmSampToListNum = @(x) vertcat(samplesToRemove, participant.PARENT_SAMPLE_NAME(x == max(x)) );
    
    % Identify samples to be removed
    for j = 1:length(numVarsToCheck)
        if over1SampLeft(samplesToRemove) == true
            input = str2double(participant.(numVarsToCheck{j}));
            if canSampBeRemNum(input) == true
                samplesToRemove = appendRmSampToListNum(input);
            end
        end
    end

    % Append samples to list 
    sampsToRemoveList{i} = samplesToRemove;
end

% Expand cell array to string vector
sampsToRemove = string(vertcat(sampsToRemoveList{:}));

% Remove duplicate samples 
metabolomicsMerged(matches(metabolomicsMerged.PARENT_SAMPLE_NAME,sampsToRemove),:) = [];

% Combine metabolon samples with metadata table
metadataPlasmaMets = outerjoin(metadata,metabolomicsMerged,'Keys','NACCID','MergeKeys',true,'Type','left');

% Save metabolomics results with associated faecal sample IDs
metVars = find(matches(metadataPlasmaMets.Properties.VariableNames,vmhMets(:,1)));
varsToSave = [1, metVars];
metabolomicsProcessed = metadataPlasmaMets(:,varsToSave);


% Find the number of samples in the metabolomics dat
message = strcat("There are ",string(sum(~all(isnan(metadataPlasmaMets{:,metVars}),2))), " mapped samples with metabolomics data");
disp(message)

savePath = fullfile(saveFolder,'metabolonProcessed.csv');
writetable(metabolomicsProcessed,savePath)
end

