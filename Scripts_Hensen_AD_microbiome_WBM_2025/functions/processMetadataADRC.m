function updatedMetadata = processMetadataADRC(metadataPath, microbiomePath)
% processMetadata checks if the required variables, ID and Sex, are
% included in the metadata. If the ID column is called, sample, name,
% sample_id, or samples_name, the variable is renamed to ID. If the
% variable Sex encodes Sex information in the form, m/f, the entries are
% translated to "male" and "female". An error is thrown if the metadata
% contains names or formats outside the allowed formats.
%
% INPUT
% metadataPath          Path to the metadata file
% 
% OUTPUT
% updatedMetadataPath   Path to the processed metadata file
%
% AUTHOR: Tim Hensen, July 2024

% Read the metadata file
metadata = readtable(metadataPath,'VariableNamingRule','preserve');

% Identify sample ID column
metadata = renamevars(metadata,'Specimen.Bar.Code','ID');
metadata = movevars(metadata,'ID','Before',1);

% To ensure that the metadata sample IDs overlap with the gut microbiome
% samples, we will now remove any trailing zeros in the metadata IDs
% Get metadata sample names. Also, I will append an X to mitigate matlab
% variable name handling later in the pipeline.
metadata.ID = append('X',regexprep(string(metadata.ID),'0+$',''));

% Filter metadata on gut microbiome samples
microbiomeSamples = readcell(microbiomePath,'Range', 'B1:ZZZ1'); % Read all sample names and exclude the row name
metadata = metadata(ismember(metadata.ID,microbiomeSamples),:);

% Convert clinical marker variables to numeric values 

% First, convert the Ptau variable to numerical values
ptauNum = nan(size(metadata.("Ptau.181")));
ptauNum(matches(metadata.("Ptau.181"),'pos')) = 1;
ptauNum(matches(metadata.("Ptau.181"),'neg')) = -1;
metadata.("Ptau.181") = ptauNum;

% Also convert visit information
metadata = renamevars(metadata,'Visit','Visit_year');
metadata.Visit_year = str2double(erase(metadata.Visit_year,'Y'));

% Then convert the remaining clinical markers to numerical values
metadata = convertvars(metadata, {'HD.X.pTau181','HD.X.NFL','HD.X.GFAP','HD.X.ABeta40','HD.X.ABeta42'}, 'str2double');

% Set all clinical marker data where no information is available to nan
clinMarkers = {'AMYLPET','AMYLCSF','FDGAD','HIPPATR','TAUPETAD','CSFTAU','FDGFTLD','TPETFTLD'};
metadata{:,clinMarkers}(metadata{:,clinMarkers}==8) = nan;

% Set all variables where no data was available to nan
numericMD = metadata{:,vartype('numeric')};
numericMD(ismember(numericMD, [-4, 99, 999, 9999,888.8000])) = NaN; % [-4, 99, 999, 9999] indicate values without data
metadata{:,vartype('numeric')} = numericMD;

% Remove columns with no information 
metadata(:,all(ismissing(metadata))) = [];

% Define variables of interest to decode
varsOfInterest = {'SEX','NACCUDSD', 'NACCALZD', 'NACCAPOE'};

% Convert variables of interest to categorical variables
metadata = convertvars(metadata, varsOfInterest, 'categorical');

% Decode variables of interest
metadata.Sex = renamecats(metadata.SEX, {'1','2'}, {'male','female'}); % Rename sex variable for persephone
% metadata.NACCUDSD = renamecats(metadata.NACCUDSD, {'1','2','3','4'}, {'normal cognition','impaired non-mci','MCI','Dementia'}); % Dementia status
metadata.NACCUDSD = renamecats(metadata.NACCUDSD, {'1','2','3','4'}, {'Healthy controls','impaired non-mci','MCI','AD dementia'}); % Dementia status
metadata.NACCALZD = renamecats(metadata.NACCALZD, {'0','1','8'}, {'No AD diagnosis','AD diagnosis','No cognitive impairment'}); % AD status
metadata.NACCAPOE = renamecats(metadata.NACCAPOE, {'1','2','3','4','5','6','9'}, ...
                              {'APOE e3/e3','APOE e3/e4','APOE e2/e3','APOE e4/e4','APOE e2/e4','APOE e2/e2','NaN'}); % APOE information

% Prepare response variables of interest
metadata.AD = string(nan(height(metadata),1));
metadata.AD(matches(string(metadata.NACCALZD),'AD diagnosis')) = 'AD';
metadata.AD(~matches(string(metadata.NACCALZD),'AD diagnosis')) = 'no_AD';

% Group APOE samples by E4 status and by Allele presence. E2/4 will be set to NaN. 

metadata.NACCAPOE = string(metadata.NACCAPOE);
E4_presence = matches(metadata.NACCAPOE,{'APOE e2/e4','APOE e3/e4','APOE e4/e4'}); % E4 presence
metadata.APOE_E4 = string(nan(height(metadata),1));
metadata.APOE_E4(E4_presence) = 'E4';
metadata.APOE_E4(~E4_presence) = 'NO_E4';
metadata.APOE_E4 = categorical(metadata.APOE_E4,{'E4','NO_E4'});

% Add APOE allele status (APOE E2/4 will be excluded)
metadata.APOE_ALLELE = string(nan(height(metadata),1));
metadata.APOE_ALLELE( matches(metadata.NACCAPOE, {'APOE e2/e2','APOE e2/e3'}) ) = 'E2';
metadata.APOE_ALLELE( matches(metadata.NACCAPOE, 'APOE e3/e3') ) = 'E3';
%metadata.APOE_ALLELE( matches(metadata.NACCAPOE, {'APOE e2/e4','APOE e3/e4','APOE e4/e4'}) ) = 'E4';
metadata.APOE_ALLELE( matches(metadata.NACCAPOE, {'APOE e3/e4','APOE e4/e4'}) ) = 'E4';

% Prepare apoe group data
metadata.APOE_ALLELE = categorical(metadata.APOE_ALLELE,{'E2','E3','E4'});
metadata.APOE_ALLELE = renamecats( metadata.APOE_ALLELE,{'E2','E3','E4'},{'ϵ2','ϵ3','ϵ4'});

% Aim: Combine metadata information of relevant covariates 
% into the following groups:
% Neuropsychiatric symptoms (y/n) -> depression or anxiety or PTSD or other
% Smoking (y/n) -> ?
% Alcoholism (y/n) -> 
% Cardiovascular problems (y/n) -> 

%%% Neuropsychiatric problems %%%
% Extract and combine neuropsychiatric data

% INSOMN: Insomnia/hyposomnia
% DEP: Presumptive etiologic diagnosis - Depression
% ANXIET: Presumptive etiologic diagnosis - Anxiety
% BIPOLDX: Presumptive etiologic diagnosis - Bipolar disorder
% SCHIZOP: Presumptive etiologic diagnosis - Schizophrenia or other psychosis
% DELIR: Presumptive etiologic diagnosis - Delirium
% PTSDDX: Presumptive etiologic diagnosis - Post-traumatic stress disorder (PTSD)
% OTHPSY: Presumptive etiologic diagnosis - Other psychiatric disease
% NACCADEP: Reported use of antidepressants
% NACCAPSY: Reported use of anti-anxiety medications
% NACCAANX: Reported current use of an anxiolytic, sedative, or hypnotic agent
nps = metadata(:,{'INSOMN','DEP','DEP2YRS','OCD','ANXIET','ANXIETY','BIPOLAR', ...
    'BIPOLDX','SCHIZ','PTSD','SCHIZOP','DELIR','PTSDDX','OTHPSY','PSYCDIS', ...
    'NACCADEP','NACCAANX','NACCAPSY'});

metadata.NPS = any(nps{:,:}==1,2); % Add binary variable indicating if a participant had any neuropsychiatric problems. (1 = yes)

%%% Smoking %%%

% Extract and combine smoking information
% TOBAC30: Smoked cigarettes in last 30 days
% TOBAC100: Smoked more than 100 cigarettes in life
% SMOKYRS: Total years smoked cigarettes
% PACKSPER: Average number of packs smoked per day
% QUITSMOK: If the subject quit smoking, age at which he/she last smoked (i.e., quit)
smoking = metadata(:,{'TOBAC30','TOBAC100','SMOKYRS','PACKSPER','QUITSMOK'});
smoking = fillmissing(smoking,'constant',0);

% Lifetime smoking is associated with cognitive decline
% (https://doi.org/10.1093/aje/kwm116). 'TOBAC100' and 'SMOKYRS' should be
% included when assessing cognitive scores. 
% However, after quitting, the gut microbiome restores to the non-smoking
% state (https://doi.org/10.1038/s41598-022-10132-z). So, I want to only
% include TOBAC30 in the AD and MCI analyses.


%%% Alcohol abuse %%%

% Extract and combine alcohol information:
% ALCOCCAS: In the past three months, has the subject consumed any alcohol?
% ALCFREQ: During the past three months, how often did the subject have at least one drink of any alcoholic beverage such as wine, beer, malt liquor, or spirits?
% ALCOHOL: Alcohol abuse - clinically significant impairment occurring over a 12-month period manifested in one of the following areas: work, driving, legal, or social

alcohol = metadata(:,{'ALCOCCAS','ALCFREQ','ALCOHOL'}); %alcohol = fillmissing(alcohol,'constant',0);

% Alcohol abuse thresholds: https://www.niaaa.nih.gov/alcohols-effects-health/alcohol-drinking-patterns
% ALCFREQ = 4 -> Daily or almost daily drinking. I will use this as a
% cutoff for heavy alcohol use.
metadata.DAILY_ALCOHOL = metadata.ALCFREQ==4;


% Extract cardiovascular problems and combine data into a binary variable.
% I will split this part in hypertension related problems and broader
% cardiovascular problems. 

% Hypertension:
% NACCAHTN: Reported current use of any type of an antihypertensive or blood pressure medication
% NACCHTNC: Reported current use of an antihypertensive combination therapy
% HYPERTEN: Hypertension
% HYPERT: Hypertension present
% NACCVASD: Reported current use of a vasodilator

hypertension = metadata(:,{'NACCAHTN','NACCHTNC','HYPERTEN','HYPERT','NACCVASD'});

% Some individuals without active hypertension still use antihypertensive
% drugs. Many individuals have hypertension. I will choose to use just the
% HYPERT variable as a covariate in the results.
% HYPERT

% Cardiovascular problems:
% CVHATT: Heart attack/cardiac arrest
% HATTMULT: More than one heart attack/cardiac arrest?
% CVBYPASS: Cardiac bypass procedure
% CVOTHR: Other cardiovascular disease
% CVOTHRX: Specification for other cardiovascular disease
% CONGHRT: Congestive heart failure present
% HVALVE: Procedure: heart valve replacement or repair within the past 12 months
cardiovasc = metadata(:,{'CVHATT','HATTMULT','CVBYPASS','CVOTHR','CVOTHRX','CONGHRT','HVALVE'});
cardiovasc.CVOTHRX = [];
cardiovasc = fillmissing(cardiovasc,'constant',0);
metadata.HEART_DISEASE = any(cardiovasc{:,:}==1,2); % Add information. There are only 19 individuals with a heart disease.

% What is the overlap between heart failure and hypertension?
% cardProblems = table(hypertension.HYPERT,cardiovasc.any);
% cardProblems.Properties.VariableNames = {'hypertension','heart failure'};
% figure; heatmap(cardProblems,'hypertension','heart failure')
% Of the 19 people with heart failure, 11 of them also had hypertension.

% Conclusion: The following covariates will be controlled for, 
% AD and dementia status: NPS, TOBAC30, DAILY_ALCOHOL, HYPERT, (HEART_DISEASE)
% Cognition: NPS, TOBAC100, DAILY_ALCOHOL, HYPERT, (HEART_DISEASE)

updatedMetadata = metadata;
end

%%
% metaNum  = metadata(:,vartype('numeric'));
% %metaNum(:,{'Kit.Number','Specimen.Bar.Code','EDUC','NACCAGE','FORMVER'})
% 
% minFourMetaNum = table2array(metaNum);
% noDataCode = [-4, 99, 999, 9999, nan];
% for i = 1:length(noDataCode)
%     minFourMetaNum( minFourMetaNum == noDataCode(i)) = nan;
% end
% minFourMetaNum(isnan(minFourMetaNum)) = 0;
% minFourMetaNum( minFourMetaNum~=0) = 1;
% 
% % Find the number of samples with data for each variable 
% dataSum = table(metaNum.Properties.VariableNames',sum(minFourMetaNum,1)','VariableNames',{'Variable','Number of measurements'});
% dataSum = sortrows(dataSum,'Number of measurements','descend');
% 
% figure;
% imagesc(minFourMetaNum)
% colorbar