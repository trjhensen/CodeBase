function bootMeanTable = getDietRCstats(dietSensFolder, bootSamp)
% Get the mean average diet reaction reduced cost value
% INPUT:
% dietSensFolder
% bootSamp

% Get paths to results
metNames = {dir( fullfile(dietSensFolder,'*.csv') ).name};
rcPaths = fullfile(dietSensFolder, metNames)';

bootMeanTables = cell(1,length(rcPaths));
parfor i=1:length(rcPaths)
    bootMeanTables{i} = getDietRCstatistics(rcPaths{i}, erase(metNames(i),'.csv'), bootSamp);
end

% Create tall table with 
bootMeanTable = vertcat(bootMeanTables{:});

if 0
    % Investigate the importance of each metabolite graphically
    meanRcMat = cellfun(@(x) x.Mean,bootMeanTables,'UniformOutput',false);
    meanRc = array2table(horzcat(meanRcMat{:}), 'RowNames',bootMeanTables{1}.("Diet reaction"),'VariableNames',metNames);

    [newrowOrder, newcolOrder] = clustRowsAndCols(meanRc{:,:});
    meanRc = meanRc(newrowOrder,newcolOrder);
    
    figure; imagesc(meanRc{:,:}'); colorbar;
    xticks(1:height(meanRc))
    yticks(1:width(meanRc))
    xticklabels(meanRc.Properties.RowNames)
    yticklabels(meanRc.Properties.VariableNames)
end
end

function bootMeanTable = getDietRCstatistics(rcPath, metabolite, bootSamp)
% Load values and process table
contributions = readtable(rcPath, 'PreserveVariableNames', true,'ReadRowNames',true); 
contributions = fillmissing(contributions, 'constant', 0); % Set all nan values to zeros
contributions = convertvars(contributions,contributions.Properties.VariableNames,@(x) round(x,6));
% contributions = removevars(contributions, all(contributions{:,:}==0) ); % Remove variables with zeros in all samples

% Find the average flux contributions 
bootMeanTable = getBootMeanTable(contributions, bootSamp);

% Add metabolite names
metRep = repmat(metabolite, height(bootMeanTable),1);
bootNumRep = repmat(bootSamp, height(bootMeanTable),1);
bootMeanTable = addvars(bootMeanTable,metRep, bootNumRep,'NewVariableNames',{'Objective','Bootstrap samples'},'After',1);
% bootMeanTable = sortrows(bootMeanTable,'Mean','descend');
end


function bootMeanTable = getBootMeanTable(contributions, bootSamp)
% Calculate the mean flux contributions and the associated 95% confidence
% intervals
contributionMatrix = table2array(contributions);
[ci,bootstat] = arrayfun(@(x) bootci(bootSamp,{@mean,contributionMatrix(:,x)},'Alpha',0.05), 1:size(contributions,2),'UniformOutput',false);

% Check if the 95% CI crossses zero
pValues = cellfun(@(x) 2 * min(mean(x >= 0, 1), mean(x <= 0, 1)), bootstat)';

% Generate table with summary statistics
bootMeanTable = array2table([cellfun(@mean,bootstat)',[ci{:}]',pValues],'VariableNames',{'Mean','2.5CI','97.5CI','pValue'});

% Add information on dietary reactions and the objective reaction information
bootMeanTable = addvars(bootMeanTable,contributions.Properties.VariableNames','NewVariableNames','Diet reaction','Before',1);
end