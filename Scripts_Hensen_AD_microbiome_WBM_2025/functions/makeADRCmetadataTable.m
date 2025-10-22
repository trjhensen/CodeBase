function tab = makeADRCmetadataTable(metadata)
% Personal function for the ADRC project by Tim Hensen.
% This function creates a summary statistics table for the manuscript.  
% Table format: 
% Columns: Variable, N, CN (N=), MCI (N=), Dementia (N=)
% Rows: AD diagnosis, Age, Sex female, no. (%), Education in years, BMI,
% DAILY_ALCOHOL, Hypertension, NPS, APOE genotype (E4), Global cognition, p-value

% Overall: [mean(metadata.age_at_collection),std(metadata.age_at_collection)]

% Get variables of interest
vars = {'NACCUDSD','Sex','age_at_collection','EDUC','NACCBMI','DAILY_ALCOHOL','HYPERT','NPS','NACCAPOE','NACCMOCA','G','APOE_E4'};
metadata = metadata(:,vars);

% Generate a binary variable for each APOE genotype
% Assume T is your table and T.APOE is a string array or cellstr
categories = ["APOE e2/e2", "APOE e2/e3", "APOE e2/e4", ...
              "APOE e3/e3", "APOE e3/e4", "APOE e4/e4"];

% Loop through each category and create a binary variable
for i = 1:numel(categories)
    genotype = double(matches(metadata.NACCAPOE,categories(i)));
    genotype(matches(metadata.NACCAPOE,'NaN')) = nan;
    metadata.(categories(i)) = matches(metadata.NACCAPOE,categories(i));
end

% Remove samples for which no NACCUDSD is available
%metadata = metadata(~matches(metadata.NACCUDSD,'<undefined>'),:);

% Generate binary variables
metadata.Sex = matches(metadata.Sex,"female");
metadata.APOE_E4 = matches(metadata.APOE_E4,"E4");
% metadata.AD = matches(metadata.AD,"AD");
binaryVars = {'Sex','DAILY_ALCOHOL','HYPERT','NPS'};
metadata = convertvars(metadata,binaryVars,'logical');

% Identify all binary variables
binaryVars = metadata(:,vartype('logical')).Properties.VariableNames;

% Generate empty table as follows
% Columns: Variable, N, CN (N=), MCI (N=), Dementia (N=)
% Rows: AD diagnosis, Age, Sex female, no. (%), Education in years, BMI,
% DAILY_ALCOHOL, Hypertension, NPS, APOE genotype (E4), Global cognition, p-value

% Generate the table variable names
[groups,names] = findgroups(string(metadata.NACCUDSD)); % Find dementia groups
gcounts = string(groupcounts(groups)); % Get the number of samples in each group
demNames = append(names, " (N = ", gcounts," )")';
varNames = ["Variable", "N", demNames, "P-value"]; % Create variable names array

% Generate the table row names
rowNames = string(metadata.Properties.VariableNames(2:end));


% Preallocate table with variables of interest
tab = table('Size',[numel(rowNames), numel(varNames)], ...
    'VariableTypes',[{'string'},{'double'},repmat({'string'},1,numel(names)),{'double'}],...
    'VariableNames',varNames, 'RowNames',rowNames');

% Populate first column
tab.Variable = rowNames';

metadataTmp = metadata;
% Loop through each category and create a binary variable
for i = 1:numel(categories)
    genotype = double(matches(metadataTmp.NACCAPOE,categories(i)));
    genotype(matches(metadataTmp.NACCAPOE,'NaN')) = nan;
    metadataTmp.(categories(i)) = genotype;
end

tab('NACCAPOE',:) = [];
metadata = removevars(metadata,'NACCAPOE');
metadataTmp = removevars(metadataTmp,'NACCAPOE');

% Populate the second column with the number of samples per variable
tab.N = table2array(varfun(@(x) sum(~isnan(x)), metadataTmp(:,2:end)))';

%tst = table2array(varfun(@(x) sum(~isnan(x)), metadataTmp(:,2:end)))';




% Define function for data type processing
formNum = @(x) string( round(x ,2) );

% Populate the table rows of binary variables with the number of samples
% and percentage of positive hits
fun1 = @(x) append( string(sum(x))," (" , formNum( sum(x) / numel(x) *100), "%)"); % Get the number of samples and percentages
cellRes = cellfun( @(y) splitapply(fun1,metadata.(y),groups)',binaryVars,'UniformOutput',false); % Apply this function for each dementia group and for each variable
tab{binaryVars,demNames} = vertcat(cellRes{:}); % Add the generated statistics to the summary table

% Add p-values from chi2 test of independence
[~, ~, p] = cellfun(@(x) crosstab(metadata.NACCUDSD,metadata.(x)), binaryVars,'UniformOutput',false);
tab{binaryVars,"P-value"} = cell2mat(p)';

% Populate the table rows of all continues variables with their mean (SD)
% statistics.
contVars = setdiff(metadata.Properties.VariableNames(:,2:end),binaryVars); % Find the variables with continuous data
funCont = @(x) append( formNum(mean(x,'omitmissing')), " (" , formNum( std(x,'omitmissing') ), ")"); % Define function for calculating mean (SD)
cellRes = cellfun( @(y) splitapply(funCont,metadata.(y),groups)',contVars,'UniformOutput',false); % Apply function for each dementia group and for each variable
tab{contVars,demNames} = vertcat(cellRes{:}); % Add the generated statistics to the summary table

% Generate p-values for the continues variables using a one-way anova test
% (why not kruskall-wallis?)
tab{contVars,'P-value'} = cellfun(@(x) anova1(metadata.(x),groups,"off") , contVars)';


% Process table for inclusion in manuscript:

pval = tab.("P-value"); % Process p-value formatting
vsmall = pval<1e-16; % Find very small p-values 
pval = arrayfun(@(x) string(sprintf('%.2e', x)), pval); % Format p-values to scientific notation
pval(vsmall) = "<<1e-16"; % Set very small p-values to <<1e-16
tab.("P-value") = pval;



% Update the row names
% tab.Variable(matches(tab.Variable,'AD')) = 'AD diagnosis';
tab.Variable(matches(tab.Variable,'Sex')) = 'Female samples';
tab.Variable(matches(tab.Variable,'age_at_collection')) = 'Age in years';
tab.Variable(matches(tab.Variable,'NACCBMI')) = 'BMI in kg/m2';
tab.Variable(matches(tab.Variable,'EDUC')) = 'Years of education';
tab.Variable(matches(tab.Variable,'DAILY_ALCOHOL')) = 'Alcohol abuse';
tab.Variable(matches(tab.Variable,'HYPERT')) = 'Hypertension';
tab.Variable(matches(tab.Variable,'NPS')) = 'Neuropsychiatric problems';
tab.Variable(matches(tab.Variable,'NACCMOCA')) = 'MOCA test score';
tab.Variable(matches(tab.Variable,'G')) = 'Global cognition score';


countVarIdx = ~matches(tab.Variable,{'Age in years','BMI in kg/m2','Years of education','MOCA test score','Global cognition score'});
tab.Variable(countVarIdx) = append(tab.Variable(countVarIdx), ', no. (%)');

% Reorder variables
varNamesReorder = [varNames(1:2), flip(varNames(3:end-1)), varNames(end)];
tab = tab(:,varNamesReorder);
end