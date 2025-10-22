% performADRCanalysis
clc
% inputs
% preparedInputTable, preparedMetadata, preparedRawInputTable


% RQ1: How do the predicted fluxes differ between dementia patients and
% controls, MCI patients, and controls, and dementia patients and MCI
% patients. 

% Prepare the three input tables
inputDemCn = preparedMetadata(matches(preparedMetadata.NACCUDSD,{'Dementia','normal cognition'}),:); % Filter on healthy and dementia
inputMciCn = preparedMetadata(matches(preparedMetadata.NACCUDSD,{'MCI','normal cognition'}),:); % Filter on healthy and MCI
inputDemMci = preparedMetadata(matches(preparedMetadata.NACCUDSD,{'Dementia','MCI'}),:); % Filter on Dementia and MCI

% State response variable of interest
response = 'AD';

% State covariates of interest for Dementia analysis
confounders = {...
    ... % Demographics
    'Sex',... % Male/Female    
    'age_at_collection',... % Age in years
    ... % Technical covariates
    ...%'Visit_year',...
    'mapped_species_reads',...
    'lane',... % Sequencing run lane
    ...'DaysCollectiontoProcessing',...
    'HYPERT',... % Hypertension present
    'Ethanol_added',...% Daily alcohol consumption
    'DAILY_ALCOHOL',...
    ...%'TOBAC100',...
    ...%'NACCUDSD',...
    'APOE_E4',... % E4 allele in APOE gene present
    ...%'HEART_DISEASE',... % Cardiovascular disease present.
    'NPS',... % Neuropsychiatric symptoms     
    'NACCBMI',...% body mass index
    ... % Comorbidities
    ...%'TOBAC30',...% Smoking in the past 30 days
    };

% Note: BMI reduced the explained variance. NPS increased the explained
% variance. Smoking drastically decreases the explained variance.
% Hypertension and alcohol consumption increase the explained variance.
% Heart disease increases the explained variance a bit, but I am still
% considering to remove it. The mapped reads does not seem to do anything.
% I might replace it with the total species reads (total_reads). The lane
% covariate actually improves the signal nicely.
% Adding BMI reduces the explained variance of the model and muddies the
% signal. I will not use BMI in the regressions. 

% nps correlates with AD status, which might result in it masking or
% explaining the found associations. The same might be said for E4 status.
% Hypertension however clarifies the associations with AD status.

%
% Lets investigate the influence of each confounder on the results by
% iteratively adding a confounder

% What if I remove MCI individuals?
preparedMetadata1 = preparedMetadata(matches(preparedMetadata.NACCUDSD,{'Dementia','normal cognition'}),:);

preparedMetadata1 = preparedMetadata;
%preparedMetadata1 = renamevars(preparedMetadata1,'Ethanol added (Y/N)','ETHANOL');

% Note: amet is more associated with AD if you remove all MCI cases.

regCell = cell(length(confounders)+1,3);
confounderSetup = {};
counter = 0;

for i=1:length(confounders)+1
    i
    confounderList = confounderSetup;

    % Get regression results
    [regCell{i,1},regCell{i,2}] = testConfounders(preparedInputTable, preparedMetadata1, response, confounderList);
    regCell{i,3} = regCell{i,1}.Formula(1,:);

    % update counter
    counter = counter + 1;
    if counter == length(confounders)+1
        break;
    end
    
    confounderSetup = confounders(1:counter);
end

%regCellFullSet = regCell;


%% Next, I want to perform a stepwise glm for AD status on the arginine flux.
clc

response = 'APOE_E4';
varsOfInterest = {'ID',...
    ... % Demographics
    'Sex',... % Male/Female    
    'age_at_collection',... % Age in years
    'NACCBMI',... % body mass index
    ...%'Visit_year',...
    ...%'APOE_E4',... % E4 allele in APOE gene present
    ... % Comorbidities
    'NPS',... % Neuropsychiatric symptoms 
    ...%'TOBAC30',...% Smoking in the past 30 days
    'DAILY_ALCOHOL',...% Daily alcohol consumption
    'HYPERT',... % Hypertension present
    'HEART_DISEASE',... % Cardiovascular disease present.
    ... % Technical covariates
    'mapped_species_reads',... % Total mapped reads
    'lane',... % Sequencing run lane
    ...%'NACCUDSD',...
    response...
    };



rxn = 'DM_tdca3s[bc]';

metadataVarsToCheck = preparedMetadata(:,varsOfInterest);
dataForSWglm = outerjoin(preparedInputTable(:,{'ID',rxn}), metadataVarsToCheck,'Keys','ID','MergeKeys',true);
dataForSWglm.ID = [];
dataForSWglm = movevars(dataForSWglm,rxn,'After',response);

dataForSWglm.(string(response)) = grp2idx(categorical(dataForSWglm.(string(response))))-1;

% Perform stepwise glm on arginine fluxes
mdl = stepwiseglm(dataForSWglm,'linear','Distribution','normal');
mdl
mdl.Rsquared.Ordinary



function [fluxRes,regressions] = testConfounders(input, metadata, response, confounderList)

    % Prepare regression formula
    if isempty(confounderList)
        regFormula = strcat(response,'~','Flux');
    else
        regFormula = strcat(response,'~','Flux+',strjoin(confounderList,'+'));
    end
    
    % Perform regressions
    [results,regressions] = performRegressions(input,metadata,regFormula);
    % Create output
    fluxRes = results.Flux;
    % Prune tables
    fluxRes = fluxRes(:,{'Reaction','Formula','N','estimate','pValue','FDR','SE','R2'});
end