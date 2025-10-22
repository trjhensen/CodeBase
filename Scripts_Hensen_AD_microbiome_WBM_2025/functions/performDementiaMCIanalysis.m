% function performDementiaMCIanalysis
%clc
% inputs
% preparedInputTable
% preparedMetadata
clc
% State response variable of interest
response = 'NACCUDSD';

% State covariates to control for
confounders = {...
    'Sex',... % Male/Female    
    'age_at_collection',... % Age in years
     'mapped_species_reads',...
     'lane',... % Sequencing run lane
     'Ethanol_added',...
    ...'HYPERT',... % Hypertension present
    'DAILY_ALCOHOL',...% Daily alcohol consumption
    ...'APOE_E4'...
    };

% Prepare metadata by transforming the response variable into a numerical
% variable.
metadata = preparedMetadata;

% Remove missing samples and impaired non-mci samples
metadata(matches(metadata.NACCUDSD,{'<undefined>'}),:) = [];

% Make relevant metadata categorical
metadata = convertvars(metadata,{'Sex'},'categorical'); % I do not need to specify extra if a variable is categorical in the regression

% Create single table for multinomial regression (healthy, vs MCI, vs
% Dementia).
mets = 'amet';
metadata = metadata(:,[{'ID'},confounders, cellstr(response)]);
inputData = preparedInputTable(:,{'ID',strcat('DM_',mets,'[bc]')});
inputData.Properties.VariableNames(2) = {'Flux'};
mergedData = outerjoin(metadata,inputData,'Keys','ID','MergeKeys',true,'Type','left');

% Define regression formula
regFormula = strcat(response,'~','Flux+',strjoin(confounders,'+'));

% Create categorical response variable
% cogDecline = zeros(height(mergedData),1);
% cogDecline(matches(mergedData.(response),'normal cognition')) = 0;
% cogDecline(matches(mergedData.(response),'MCI')) = 1;
% cogDecline(matches(mergedData.(response),'Dementia')) = 2;
% metadata.(response) = categorical(cogDecline);

%


% Compare multinomial regression to two separate logistic regressions. How
% do the results differ?

% mciInput = mergedData(matches(mergedData.(response),{'normal cognition','MCI'}),:);
% mciInput.(string(response)) = grp2idx(categorical( mciInput.(string(response)),{'normal cognition','MCI'}  ) )-1;
% 
% mdlMCI = fitglm(mciInput,regFormula,'Distribution','binomial','LikelihoodPenalty','jeffreys-prior'); % Firth's regression
% mdlMCI
% 
% 
% dmInput = mergedData(matches(mergedData.(response),{'normal cognition','Dementia'}),:);
% dmInput1 = dmInput;
% dmInput.(string(response)) = grp2idx(categorical( dmInput.(string(response)),{'normal cognition','Dementia'}  ) )-1;
% 
% mdldm = fitglm(dmInput,regFormula,'Distribution','binomial','LikelihoodPenalty','jeffreys-prior'); % Firth's regression
% mdldm



% Perform multinomial regression
% mergedData1 = mergedData;
mergedData.(string(response)) = categorical( mergedData.(string(response)),{'Dementia','MCI','impaired non-mci','normal cognition'}  );
% 
% multReg = fitmnr(mergedData,regFormula);
% % multReg.Rsquared.Ordinary
% 

%

MnrMdl = fitmnr(mergedData,regFormula,'ModelType','ordinal','Link','probit')
%figure; plotSlice(MnrMdl)
MnrMdl.Rsquared.Ordinary





