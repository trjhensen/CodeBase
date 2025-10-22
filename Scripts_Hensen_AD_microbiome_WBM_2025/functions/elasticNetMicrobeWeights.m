function [enetTab, enet_stats] = elasticNetMicrobeWeights(fluxMicrobeData)
% This function finds the weight coefficients of microbial species for flux
% predictions using elastic net regression

alpha = 0.5; % Elastic net ridge parameters

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

% Preallocate table for elastic net results
enetTab = cellfun(@(x,y) table(x, y', zeros(length(y'),1),'VariableNames',{'Metabolite','Taxa','Beta'}),metNameForTable, microbeNames, 'UniformOutput',false);

% Preallocate table with summary statistics
numMets = length(metNames);
varNames = {'MSE','R2','DoF','Alpha','minLambda'};
enet_stats = array2table(zeros(numMets,length(varNames)),'VariableNames',varNames,'RowNames',string(metNames)');

% Populate alpha variable
enet_stats.Alpha = repmat(alpha,height(enet_stats),1);

% Warm up algorithm
lasso(predictors{1}, responses{1}, 'CV',10,'Intercept',false,'Standardize',false,'Alpha',alpha);
lasso(predictors{2}, responses{2}, 'CV',10,'Intercept',false,'Standardize',false,'Alpha',alpha);

for i=1:numMets
    % Get predictor and response vars for metabolite i
    pred = predictors{i};
    resp = responses{i};
    
    % Remove nan values
    nanSamples = isnan(resp);
    y = resp(~nanSamples);
    X = pred(~nanSamples,:);
    
    % Perform lasso regression on random sample
    [b,lStat] = lasso(X, y, 'CV',10,'Intercept',false,'Standardize',false,'Alpha',alpha); % Data is already normalized
    
    % Find the coefficients for each microbe and add them to the table
    coef = b(:, lStat.IndexMinMSE);
    enetTab{i}.Beta = coef;
    
    % Populate statistics table
    enet_stats.MSE(i) = lStat.MSE(lStat.IndexMinMSE); % MSE
    enet_stats.DoF(i) = lStat.DF(lStat.IndexMinMSE); % Degrees of freedom
    enet_stats.minLambda(i) = lStat.LambdaMinMSE;
    
    % Add the R2 value
    y_pred = X * coef + lStat.Intercept(lStat.IndexMinMSE); % Predict using the selected coefficients and intercept
    SST = sum((y - mean(y)).^2);         % total sum of squares (SST) in y
    SSR = sum((y - y_pred).^2);          % residual sum of squares (SSR) variance
    enet_stats.R2(i) = 1 - SSR / SST;   % R2
end

end
