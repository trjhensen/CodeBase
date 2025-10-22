function lassoRes = pruneMicrobialFeaturesWithLasso(fluxMicrobeData, nBootLasso)
%
% INPUTS:
% fluxMicrobeData       Cell array with the microbe relative 
%                       abundances from the WBMs, sample sex 
%                       information, and a final column with 
%                       the associated flux predictions.
%
% OUTPUTS:
% lassoResFiltered      Cell array with a list of microbial 
%                       species for which non-zero beta coefficients 
%                       were obtained when applying a lambda 
%                       regularization estimator that resulted 
%                       in the smallest mean square error (MSE) 
%                       for the flux prediction, while controlling 
%                       for sex.


% Encode sex information into binary variable (0/1) and apply function to
% each cell.
cf = @(f,x) cellfun(f, x, 'UniformOutput',false); % Preallocate cellfun function for compact code
convertSex = @(x) convertvars(x,'Sex', @(x) grp2idx(categorical(x))-1);
microbeDataCell = cf(convertSex, fluxMicrobeData);

% Get predictors and response variables and z-scaled the data column-wise
transforDat = @(x) normalize(table2array(x));
predictors = cf(@(x) transforDat(x(:,1:end-1)), microbeDataCell);
responses = cf(@(x) transforDat(x(:,end)), microbeDataCell);

% Fetch microbe and metabolite names
microbeNames = cf(@(x) x.Properties.VariableNames(1:end-1), fluxMicrobeData); % Get the microbe names
metNames = cf(@(x) x.Properties.VariableNames(end), fluxMicrobeData); % Get the metabolite names

% Preallocate replicated metabolite names for lasso table annotation
metNameForTable = cellfun(@(x,y) repmat(x,length(y),1), metNames, microbeNames,'UniformOutput',false);

% Preallocate table for annotated selection frequencies
lassoRes = cellfun(@(x,y) table(x, y','VariableNames',{'Reaction','Species'}),metNameForTable,microbeNames, 'UniformOutput',false);

% Perform lasso regression on each cell
numMets = numel(responses);
numSamp = height(fluxMicrobeData{1});
tic
fprintf("> Microbe selection in progress... \n")
for i = 1:numMets

    % Get predictor and response vars for metabolite i
    pred = predictors{i};
    resp = responses{i};

    % Create selection of random samples that take 70% of samples each time
    rsampIdx = arrayfun(@(x) randsample(numSamp, round(0.7 * numSamp),true), 1:nBootLasso, 'UniformOutput', false);
    
    % Generate samples for predictors
    X_sub = cellfun(@(x) pred(x, :), rsampIdx,'UniformOutput', false);
    y_sub = cellfun(@(x) resp(x, :), rsampIdx,'UniformOutput', false);
    
    % Preallocate result vector
    selection_count = zeros(size(pred,2), 1); 
    parfor j = 1:nBootLasso
        % Perform lasso regression on random sample
        [b,lStat] = lasso(X_sub{j},y_sub{j},'CV',10,'Intercept',false,'Standardize',false); % Data is already normalized
    
        % Use the lambda with minimum CV error
        coef = b(:, lStat.IndexMinMSE);
        
        % Record which variables were selected
        selection_count = selection_count + (coef ~= 0);
    end
    
    % Add the selection frequencies to table
    lassoRes{i}.selection_frequency = selection_count / nBootLasso;

    % Display progress
    disp(append("> ",string(round(i/numMets*100,0)),"%"))
end
toc

end

%%
% Get indices for minimum mse 
% indexL = cf(@(x) x.IndexMinMSE,lassoStats);
% 
% % Extract coefficients for minimal mse and corresponding lambda
% microbeLassoTable = cell(1,numMets);
% for i=1:numMets
%     lsMet = lassoStats{i};
%     b = B{i}(:,indexL{i}); % Get coefficients of best fit
%     nzBeta = b~=0; % Get non-zero coefficients
%     % Create table with identified microbes
%     bPruned = b(nzBeta);
%     filteredMicrobes = lsMet.PredictorNames(nzBeta)';
%     metNameForTable = repmat(metNames{i},length(bPruned),1);
%     % Construct table
%     lassoTab = table(metNameForTable, filteredMicrobes, bPruned,'VariableNames',{'Reaction','Species','Coefficient'});
%     lassoTab(matches(lassoTab.Species,'Sex'),:)=[]; % Remove sex from table
%     microbeLassoTable{i} = lassoTab;
% end
% 
% % Then extract the associated lasso statistics for a second table
% sumStatsLasso = table(...
%     string(metNames)',...
%     cellfun(@(x,y) x.MSE(y), lassoStats, indexL)',...
%     cellfun(@(x,y) x.Lambda(y), lassoStats, indexL)',...
%     cellfun(@(x,y) x.DF(y), lassoStats, indexL)',...
%     'VariableNames',{'Reaction','MSE','Lambda','DF'});
% end
% 

