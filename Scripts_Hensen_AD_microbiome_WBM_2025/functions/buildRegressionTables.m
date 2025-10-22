function resTable = buildRegressionTables(rxnsOfInterest, regressionTable, regressionModels, varOfInterest)
% ADRC flux result regression table:
% inputs:
% resultsDM:
% regressionModels
% rxnsOfInterest

% Generate regression table

% 1) Extract the reactions of interest from resultsDM and the regression
% data
varsToKeep = {'Reaction','subgroups','estimate','low','high','pValue'};
resTable = regressionTable(matches(regressionTable.Reaction,rxnsOfInterest),varsToKeep); % Pruned results table

% Combine estimate, 95% CI, and p-value columns in a single string array

% Define function to format value as float above 1e-3 and as scientific notation below 1e-3
fmt = @(x,formats) arrayfun(@(y) sprintf(formats{1 + (abs(y) < 1e-3)}, y), x, 'UniformOutput', false);


if matches(varOfInterest,'AD')

    % Round numbers in table
    resTable.estimate = fmt(resTable.estimate,{'%#0.2g','%.2E'});
    resTable.low = fmt(resTable.low,{'%#0.2g','%.2E'});
    resTable.high = fmt(resTable.high,{'%#0.2g','%.2E'});
    resTable.pValue = fmt(resTable.pValue,{'%#0.2g','%.2E'});
    % resTable.estimate = arrayfun(@(x) string(sprintf('%#.2g',x)),resTable.estimate);
    % resTable.low = arrayfun(@(x) string(sprintf('%0.1G',x)),resTable.low);
    % resTable.high = arrayfun(@(x) string(sprintf('%0.1G',x)),resTable.high);
    % resTable.pValue = arrayfun(@(x) string(sprintf('%0.3f',x)), resTable.pValue);

    % Rename group variables
    resTable.subgroups = replace(resTable.subgroups,'normal cognition','Healthy controls');
else
    % Round numbers in table
    % resTable.estimate = fmt(resTable.estimate,{'%#0.2g','%.2E'});
    % resTable.low = fmt(resTable.low,{'%#0.2g','%.2E'});
    % resTable.high = fmt(resTable.high,{'%#0.2g','%.2E'});
    % resTable.pValue = fmt(resTable.pValue,{'%#0.2g','%.2E'});
    resTable.estimate = arrayfun(@(x) string(sprintf('%0.2f',x)),resTable.estimate);
    resTable.low = arrayfun(@(x) string(sprintf('%0.2f',x)),resTable.low);
    resTable.high = arrayfun(@(x) string(sprintf('%0.2f',x)),resTable.high);
    resTable.pValue = arrayfun(@(x) string(sprintf('%0.2g',x)), resTable.pValue);
end

% Combine strings in a single table
resTable = table( ...
    resTable.Reaction,...
    resTable.subgroups,...
    append(resTable.estimate, " (", resTable.low,", ",resTable.high,"); ",resTable.pValue),...
    'VariableNames',{'Reaction','subgroups','Statistics'});

% Convert original regression table to wide table
resTable = unstack(resTable,'Statistics','subgroups','VariableNamingRule','preserve');

% Reorder variables
switch varOfInterest
    case 'AD'
        resTable = resTable(:,{'Reaction','Healthy controls – MCI','MCI – AD dementia','Healthy controls – AD dementia'}); % –
    case 'APOE'
        resTable = resTable(:,{'Reaction','ϵ2 – ϵ3','ϵ3 – ϵ4','ϵ2 – ϵ4'});
end

% Next, lets identify the sample numbers from the regression data

% Get sample numbers of normal cognition
rxns = matlab.lang.makeValidName(rxnsOfInterest);
numReg = numel(regressionModels); numRxns = numel(rxns);
sampNums = cell(numRxns,numReg);

% Define formula to extract sample numbers 
switch varOfInterest
    case 'AD'
        getDat = @(x,y,z) x.(string(y)).Variables.NACCUDSD(z);
    case 'APOE'
        getDat = @(x,y,z) x.(string(y)).Variables.APOE_ALLELE(z);
end

for i=1:numReg
    for j=1:numRxns
        % Find the non-missing sample number indices
        getUsedSamps = @(x,y) x.(string(y)).ObservationInfo.Missing==0;
        samps = getUsedSamps( regressionModels{i}, rxns(j) ); 
        % Get the number of samples from the regressions
        groupInfo = getDat(regressionModels{i},rxns(j),samps);
        sampNums{j,i}(:,1) = sum( groupInfo == min(groupInfo) ); % Get the less severe disease sample numbers
        sampNums{j,i}(:,2) = sum( groupInfo == max(groupInfo) ); % Get the more severe disease sample numbers
    end
end

switch varOfInterest
    case 'AD'
        % Find the normal cognition sample numbers and add them to the regression
        % table
        resTable = addvars(resTable, cellfun(@(x) x(:,2), sampNums(:,3)), 'NewVariableNames','AD dementia N','After','Reaction');
        resTable = addvars(resTable, cellfun(@(x) x(:,1), sampNums(:,2)), 'NewVariableNames','MCI N','After','Reaction');
        resTable = addvars(resTable, cellfun(@(x) x(:,1), sampNums(:,1)), 'NewVariableNames','Healthy controls N','After','Reaction');
    case 'APOE'
        resTable = addvars(resTable, cellfun(@(x) x(:,2), sampNums(:,3)), 'NewVariableNames','ϵ4 group,N','After','Reaction');
        resTable = addvars(resTable, cellfun(@(x) x(:,1), sampNums(:,2)), 'NewVariableNames','ϵ3 group,N','After','Reaction');
        resTable = addvars(resTable, cellfun(@(x) x(:,1), sampNums(:,1)), 'NewVariableNames','ϵ2 group,N','After','Reaction');
end

if 0
    % Create nested table for space efficiency
    % Merge summary statistics
    resTable = mergevars(resTable, ...
        resTable.Properties.VariableNames(5:7),...
        "NewVariableName","Regression β (95% CI); P-value","MergeAsTable",true);
    
    % Merge count data statistics
    resTable = mergevars(resTable, ...
        resTable.Properties.VariableNames(2:4),...
        "NewVariableName","Sample count","MergeAsTable",true);
    end

% Change names of data variables
resTable.Properties.VariableNames(end-2:end) = append(resTable.Properties.VariableNames(end-2:end), ',Regression β (95% CI); P-value');

% Find metabolite names 
resTable.Reaction = renameAdrcVmhToMetName(resTable.Reaction);
% Rename Reaction to metabolite
resTable = renamevars(resTable,'Reaction','Predicted metabolite');
end