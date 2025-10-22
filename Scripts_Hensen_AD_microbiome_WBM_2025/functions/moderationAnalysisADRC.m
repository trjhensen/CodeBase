function stratResTable = moderationAnalysisADRC(preparedInputTable, preparedMetadata,regressionParam)


% Define parameters for moderation analysis
response = regressionParam.response;
confounders = regressionParam.confounders;
modVar = regressionParam.modVar;
interactionPvalThreshold = regressionParam.interactionPvalThreshold;

% If the modVar is already in the confounder list, remove this variable
% from the confounder list.
confounders(matches(confounders,modVar))=[];


% Define regression formula
regFormula = strcat(response,'~','Flux+',modVar,'+',strjoin(confounders,'+'),'+Flux:',modVar);

% Make sure that the moderation variable is not numeric
if 0%isnumeric(preparedMetadata.(modVar))
    modvarDat = cellstr(num2str(preparedMetadata.(modVar)));
    modvarDat(matches(modvarDat,{'0'})) = {'No'};
    modvarDat(matches(modvarDat,{'1'})) = {'Yes'};
    preparedMetadata.(modVar) = modvarDat;
end

% Perform regression with interaction effect

results = performRegressions(preparedInputTable,preparedMetadata,regFormula);

% Extract interaction effects 
resFields = string(fieldnames(results));
intResults = results.(resFields(end));

% Test if any further interaction analyses will be performed
proceed = false;
if ~matches(resFields,'NotDefined')
    proceed = true;
end

if proceed == true
    % Filter on reactions of interest
    intResults( isnan(intResults.estimate) | intResults.pValue > interactionPvalThreshold,:) = [];
    rxnsToTest = intResults.Reaction';

    % Test if any further interaction analyses will be performed
    proceed = false;
    if ~isempty(rxnsToTest)
        proceed = true;
    end
end

if proceed == true
    % Prune flux table
    interactionInputTable = preparedInputTable(:,["ID" rxnsToTest]);
    
    % Find binary groups in interaction variable
    [~,cohort] = findgroups(preparedMetadata.(modVar));

    
    % Stratify metadata based on group status
    metadataStrat = cellfun(@(x) preparedMetadata(matches(preparedMetadata.(modVar), x), :), cohort, 'UniformOutput', false);
    
    % Redefine regression formula
    regFormula = strcat(response,'~','Flux+',strjoin(confounders,'+')); 

    % Perform regressions on both strata
    stratRes = cellfun(@(md) performRegressions(interactionInputTable, md, regFormula),metadataStrat, 'UniformOutput', false);

    % Remove cell entry if no regressions could be performed on a stratum
    rmRes = cellfun(@(x) any(matches(fieldnames(x),{'NotDefined'})), stratRes,'UniformOutput',true);
    stratRes = stratRes(~rmRes);
    cohort = cohort(~rmRes);

    if isempty(stratRes)
        proceed = false;
    end
end

if proceed == true

    % Unnest stratified results and add cohort label
    addCohortFun = @(x,y) addvars(x.Flux,repmat(string(y),height(x.Flux),1), 'Before','N','NewVariableNames','Cohort');
    stratRes = cellfun(addCohortFun,stratRes,cohort,'UniformOutput',false);

    stratRes = cellfun(@(x) convertvars(x,'Formula','string'),stratRes,'UniformOutput',false);
    
    % Also add cohort label to the interaction results
    stratRes{3} = addvars(intResults,repmat("Full",height(intResults),1), 'Before','N','NewVariableNames','Cohort');
    
    % Create a single table for the moderation analysis results
    stratResTable = vertcat(stratRes{:});
    stratResTable = sortrows(stratResTable,'Reaction');
else
    stratResTable = table;
    disp('No interaction effect found for any of the reactions and moderation variables')
end

end

