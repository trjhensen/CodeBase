function [elasticNetResults, enet_stats, lassoResTab, microbeContributionStats] = microbeToFluxAnalysis(fluxPath, fluxLimTablePath, mappedMicrobePath, param)
% function [microbeMetPresence, metTopModels, regressionModelStats, lassoResTab] = microbeToFluxAnalysis(fluxPath, fluxLimTablePath, mappedMicrobePath, param)
% This function performs iterative backward selection on linear regressions
% to identify the smallest set of microbial species with the highest
% explained variance of the fluxes.
%

% INPUTS;
% fluxPath = paths.fluxPath;
% fluxLimTablePath = paths.apoeFM;
% mappedMicrobePath = paths.mappedMicrobePath;
% 
% OUTPUTS:
% microbeMetPresence: Most important microbes for the metabolites of
% interest.

thresholdR2 = param.thresholdR2;
minFreq = param.minFreq;
nBootLasso = param.nBootLasso;

% Filter the fluxes table on the reactions in fluxLimTable
fluxLimTable = readtable(fluxLimTablePath,'VariableNamingRule','preserve','ReadRowNames',true);
rxnsOfInterest = fluxLimTable.Properties.VariableNames;
vmhIDs = renameAdrcVmhToMetName(rxnsOfInterest, true);


% Load the processed fluxes and process sample IDs
fluxes = readtable(fluxPath,'VariableNamingRule','preserve');
fluxes = fluxes(:,[{'ID','Sex'}, vmhIDs]); % Filter on metabolites of interest
fluxes.Properties.VariableNames = [{'ID','Sex'}, rxnsOfInterest];

% Load the relative abundances
microbiome = readtable(mappedMicrobePath,'VariableNamingRule','preserve','ReadRowNames',false);
microbiome = renamevars(microbiome,'Row','ID');    
if isvar(microbiome,'Sum of taxa')
    microbiome = removevars(microbiome,'Sum of taxa');
end

% Set all nan values to zero
microbiome(:,2:end) = fillmissing(microbiome(:,2:end),'constant',0);

% Ensure that the samples are in the correct order between the fluxes and
% relative abundances
[~,ia,ib] = intersect(fluxes.ID,microbiome.ID,'stable');
fluxes = fluxes(ia,:);
microbiome = microbiome(ib,:);

% Get the associated microbial species for each reaction
fluxMicrobes = cellfun(@(x) fluxLimTable.Row(fluxLimTable.(x)==1),rxnsOfInterest,'UniformOutput',false);

% Combine tables into a single table with the fluxes at the end
mergeFun = @(x,y) outerjoin(microbiome(:,[{'ID'},y']),fluxes(:,{'ID','Sex',x}),'Keys','ID','MergeKeys',true,'Type','left');
fluxMicrobeData = cellfun(mergeFun,rxnsOfInterest,fluxMicrobes,'UniformOutput',false);

% Get the number of microbial species with any contribution to the fluxes
microbeContributionStats = cell(1,2);
microbeContributionStats{1} = getMeanSDNumMicrobes(fluxMicrobeData);

% Remove the ID column
fluxMicrobeData = cellfun(@(x) removevars(x,'ID'), fluxMicrobeData,'UniformOutput',false); 

% Find subset of microbial species that predict the associated fluxes
lassoRes = pruneMicrobialFeaturesWithLasso(fluxMicrobeData, nBootLasso);

% Create table for output
lassoResTab = vertcat(lassoRes{:});

% Create wide table for the lasso selection frequencies
lassoResTab = unstack(lassoResTab,"selection_frequency","Reaction",'VariableNamingRule','preserve');

% Filter microbes with a frequency below minFreq
lassoResFiltered = cellfun(@(x) x(x.selection_frequency>minFreq,:), lassoRes,'UniformOutput',false);

% Prune microbes from relative abundance data
microbeNames = cellfun(@(x) x.Properties.VariableNames(1:end-2), fluxMicrobeData,'UniformOutput',false); % Get the microbe names
microbesToKeep = cellfun(@(x) x.Species,lassoResFiltered,'UniformOutput',false);

% Remove microbes 
pruneMicrobes = @(x,y,z) removevars(x,setdiff(y,z) );
fluxMicrobeData = cellfun(pruneMicrobes,fluxMicrobeData, microbeNames ,microbesToKeep,'UniformOutput',false);

% Get the number of microbial contributors again after lasso-based pruning
microbeContributionStats{2} = getMeanSDNumMicrobes(fluxMicrobeData);

% Recreate list of flux-associated microbes
% fluxMicrobes = cellfun(@(x) x.Properties.VariableNames(1:end-2)',fluxMicrobeData,'UniformOutput',false);

[enetTab, enet_stats] = elasticNetMicrobeWeights(fluxMicrobeData);

%%

% Convert elastic net results to wide table with coefficients for each
% flux-microbe association
enetTabWide = unstack(vertcat(enetTab{:}),'Beta','Metabolite','VariableNamingRule','preserve');
enetTabWide(:,2:end) = fillmissing(enetTabWide(:,2:end),'constant',0);
elasticNetResults = enetTabWide;

%%


% Concatinate tables to a single table
% enetTab = vertcat(enetTab{:});

% Perform linear regression on the fluxes with covariates for all microbes
% fullMdls = cellfun(@(x) fitlm(x,'Intercept',true,'RobustOpts','on'), fluxMicrobeData,'UniformOutput',false);
% 
% % Iteratively test which microbial removal has the lowest effect on the adjusted R2
% metTopModels = cell(size(rxnsOfInterest));
% 
% parfor reactionIdx = 1:length(rxnsOfInterest)
% 
%     dataFull = fluxMicrobeData{reactionIdx};
%     microbeList = fluxMicrobes{reactionIdx};
% 
%     % Initialize data
%     currentData = dataFull;
%     remainingMicrobes = microbeList;
%     numMicrobes = numel(remainingMicrobes);
%     topAdj2 = cell(numMicrobes, 4);  % Preallocate [#microbes, adjR2, model]
% 
%     fullAdjR2 = fullMdls{reactionIdx}.Rsquared.Adjusted;
%     step = 1;
% 
%     while ~isempty(remainingMicrobes)
%         % Remove each microbe one at a time
%         removedMicrobeData = cellfun(@(microbe) removevars(currentData, microbe), ...
%                                      remainingMicrobes, 'UniformOutput', false);
% 
%         % Fit linear models
%         reducedModels = cellfun(@(data) fitlm(data, 'Intercept', true, 'RobustOpts', 'on'), ...
%                                 removedMicrobeData, 'UniformOutput', false);
% 
%         % Compute adjusted R² values
%         adjR2s = cellfun(@(model) model.Rsquared.Adjusted, reducedModels);
%         deltaAdjR2 = fullAdjR2 - adjR2s;
% 
%         % Find the best microbe to remove (least reduction in adjR²)
%         [~, bestIdx] = min(deltaAdjR2);
% 
%         % Save result
%         topAdj2{step, 1} = numel(remainingMicrobes)-1;        % Number of microbes remaining
%         topAdj2{step, 2} = adjR2s(bestIdx);                   % Best adjR² after removal
%         topAdj2{step, 3} = reducedModels{bestIdx};            % Corresponding model
% 
%         % Update current data and remove microbe from the list
%         currentData = removevars(currentData, remainingMicrobes{bestIdx});
%         remainingMicrobes(bestIdx) = [];
%         topAdj2{step, 4} = remainingMicrobes;                       % Name of microbes
%         step = step + 1;
%     end
% 
%     % Add full model as the baseline at the end
%     topAdj2 = [{numel(microbeList), fullAdjR2, fullMdls{reactionIdx}, microbeList}; ...
%         topAdj2(1:height(topAdj2)-1,:)];
% 
%     % Add the adjusted R2 value to the results
%     metTopModels{reactionIdx} = topAdj2;
% end

% Next, find the minimal set of microbial species that explained the
% fluxes. For this, I will extract the regression models with the smallest
% number of microbes that explain at least X% (80%) of the fluxes.

% convertDat = @(x) cell2mat(x(:,[1 2])); % Convert data to numerical data
% findMinCluster = @(x) min( x( x(:,2) >thresholdR2 ,1) ); % Find the smallest microbial cluster that explains most of the flux predictions
% [~,bestClusterIdx] = cellfun(@(x) findMinCluster(convertDat(x)), metTopModels,'UniformOutput',false); % Find the index of the top regression models
% topMicrobialPredictors = cellfun(@(x,y) x{y,3}, metTopModels,bestClusterIdx,'UniformOutput',false); % Extract the top regression models


% Extract the associated microbes for each flux prediction and create a
% confusion table for the microbe to flux associations.

% Find all microbes for each metabolite
% topMicrobes = cellfun(@(x) x.VariableNames(1:end-2), topMicrobialPredictors,'UniformOutput',false);
% allTopMicrobes = unique(vertcat(topMicrobes{:})); % Find the totality of microbes
% microbeMetPresence = cell2mat(cellfun(@(x) matches(allTopMicrobes,x), topMicrobes,'UniformOutput',false)); % Microbe-metabolite associations
% microbeMetPresence = array2table(microbeMetPresence,'RowNames',allTopMicrobes,'VariableNames',rxnsOfInterest); % Create table
% 
% % Annotate metTopModels with the associated reactions for easy lookup
% metTopModels = [rxnsOfInterest; metTopModels];

% Preallocate table
% numRxns = size(metTopModels,2);
% regressionModelStatsCell = cell(numRxns,1);
% for i=1:numRxns
%     fluxMicrobeModels = metTopModels{2,i}(:,3); % Extract regression model
%     rxnName = metTopModels{1,i}; % Get reaction name 
%     regressionModelStatsCell{i} = extractRegressionInfo(fluxMicrobeModels, rxnName); % Create table 
% end
% 
% % Concatinate tables to a single table
% regressionModelStats = vertcat(regressionModelStatsCell{:});

end

function regressionModelStats = extractRegressionInfo(fluxMicrobeModels, rxnName)
% Aim: To create a single regression table with
% - Reaction of interest
% - Regression formula
% - N
% - Number of microbes
% - adj R2
% - AIC
% - BIC
% - SSE

% Create empty regression table
varNames = {'Reaction','Formula','N','Taxa','R2adj','R2','AIC','BIC','RMSE','SSE'};
varTypes = [repmat({'string'},1,2),repmat({'double'},1,length(varNames)-2)];
regressionModelStats = table('Size',[length(fluxMicrobeModels),length(varNames)],'VariableTypes',varTypes,'VariableNames',varNames);

% Populate table
regressionModelStats.Reaction = repmat(string(rxnName),height(regressionModelStats),1); % Add reaction name
getForm = @(x) string(append(x.Formula.ResponseName,' ~ ', x.Formula.LinearPredictor)); % Add reaction formula
regressionModelStats.Formula = cellfun(getForm, fluxMicrobeModels); 
regressionModelStats.N = cellfun(@(x) x.NumObservations, fluxMicrobeModels); % Add the number of samples
regressionModelStats.Taxa = cellfun(@(x) x.NumPredictors-1, fluxMicrobeModels); % Add the number of microbes
regressionModelStats.R2adj = cellfun(@(x) x.Rsquared.Adjusted, fluxMicrobeModels); % Add the adjusted R-squared of regression model
regressionModelStats.R2 = cellfun(@(x) x.Rsquared.Ordinary, fluxMicrobeModels); % Add the adjusted R-squared of regression model
regressionModelStats.AIC = cellfun(@(x) x.ModelCriterion.AIC, fluxMicrobeModels); % Add the AIC index
regressionModelStats.BIC = cellfun(@(x) x.ModelCriterion.BIC, fluxMicrobeModels); % Add the BIC index
regressionModelStats.RMSE = cellfun(@(x) x.RMSE, fluxMicrobeModels); % Add the RMSE
regressionModelStats.SSE = cellfun(@(x) x.SSE, fluxMicrobeModels); % Add the SSE
end
