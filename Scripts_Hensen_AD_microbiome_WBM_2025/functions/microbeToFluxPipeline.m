function [elasticNetResults, elasticNetStats, lassoResTab, microbeContributionStats] = microbeToFluxPipeline(mContributionDir, fluxPath, mappedMicrobePath, saveDir, param)
%function [bootMeanTableFiltered, microbeMetPresence, metTopModels, regressionModelStats, lassoResTab] = microbeToFluxPipeline(mContributionDir, fluxPath, mappedMicrobePath, saveDir, param)

% Wrapper function for microbial contribution analysis ADRC project.
%%
% Extract the function parameters
rxnsOfInterest = param.rxnsOfInterest;
bootSamp = param.bootSamp;
parTrue = param.parTrue;
cumulativeFraction = param.cumulativeFraction;

% Explictely state parameters for microbe pruning
paramPruning.thresholdR2 = param.thresholdR2;
paramPruning.minFreq = param.minFreq;
paramPruning.nBootLasso = param.nBootLasso;

% Find the microbial species with the highest average flux contribution
% potentials
bootMeanTable = fluxMicrobeSensitivityAnalysisADRC_2(rxnsOfInterest, mContributionDir, bootSamp, parTrue);

% Remove microbial species with an average shadow price of zero within the
% 95% CI.
bootMeanTableFiltered = bootMeanTable( bootMeanTable.('p>0.05')==0 , :);

% Generate intersection table 

% Get microbe-metabolite pairs
metMicrobePairs = bootMeanTableFiltered(:,{'Metabolite','Taxa'}); 

% Create intersection table with taxa names as row names
metMicrobePairs.value = ones(height(metMicrobePairs),1);
metMicrobePairsWide = unstack(metMicrobePairs,"value","Metabolite",'VariableNamingRule','preserve');
metMicrobePairsWide.Properties.RowNames = metMicrobePairsWide.Taxa;

% Replace nan with zero
metMicrobePairsWide = fillmissing(metMicrobePairsWide(:,2:end),'constant',0);

% Save intersection table
fluxLimTablePath = fullfile(saveDir,'fluxMicrobeAssociations.csv');
writetable(metMicrobePairsWide,fluxLimTablePath,'WriteMode','overwrite','WriteRowNames',true)

% Filter on the most important microbial species based on the distribution
% of their bootstrapped mean average.
% [~, bootMeanTableFiltered,fluxLimTablePath,~] = filterMostImportantMicrobes_2(bootMeanTable, cumulativeFraction, saveDir);

% Start new parallel pool for microbe to flux ML analysis
% delete(gcp('nocreate'))
% parpool(feature('numCores'))

% Next, we will find the smallest subset of microbial species from within
% the filtered flux-associated microbes that best explain the variance in
% the flux predictions. 
% [microbeMetPresence, metTopModels, regressionModelStats, lassoResTab] = microbeToFluxAnalysis(fluxPath, fluxLimTablePath, mappedMicrobePath, paramPruning);
[elasticNetResults, elasticNetStats, lassoResTab, microbeContributionStats] = microbeToFluxAnalysis(fluxPath, fluxLimTablePath, mappedMicrobePath, param);

% Get the flux-associated microbes after each step
% numMicrobesLT = getMicrobeNumbers(bootMeanTable, bootMeanTableFiltered, lassoResTab, paramPruning, microbeMetPresence);

% Save the number of microbes after each step
% filePath = fullfile(saveDir,'numMicrobesForMets.csv');
% writetable(numMicrobesLT,filePath,'WriteRowNames',true)
% 
% % Save flux-microbe association results
% filePath = fullfile(saveDir,'topMicrobesForMets.csv');
% writetable(microbeMetPresence,filePath,'WriteRowNames',true)

% Save the shadow price averages for each flux-microbe pair
filePath = fullfile(saveDir,'bootciMicrobeSPvalues.csv');
writetable(bootMeanTable,filePath,'WriteRowNames',true)

% Save lasso estimation results
filePath = fullfile(saveDir,'lassoMicrobeSelectionFreq.csv');
writetable(lassoResTab,filePath,'WriteRowNames',true)

% Save elastic net estimation results
filePath = fullfile(saveDir,'enetFluxMicrobeRes.csv');
writetable(elasticNetResults,filePath,'WriteRowNames',true)

% Save summary statistics
filePath = fullfile(saveDir,'fluxMicrobeSummaryStats.xlsx');
writetable(microbeContributionStats{1},filePath,'WriteRowNames',true,'Sheet','Contributors');
writetable(microbeContributionStats{2},filePath,'WriteRowNames',true,'Sheet','Analysed');
writetable(elasticNetStats,filePath,'WriteRowNames',true,'Sheet','EnetStats');

% Save regression results
% filePath = fullfile(saveDir,'fluxMicrobeRegressionStats.csv');
% writetable(regressionModelStats,filePath)
end


function numMicrobesLT = getMicrobeNumbers(bootMeanTable, bootMeanTableFiltered, lassoResTab, paramPruning, microbeMetPresence)


% Get the number of microbes for each metabolite
numMicrobesL1 = groupcounts(bootMeanTable,"Metabolite");
numMicrobesL1 = removevars(numMicrobesL1,'Percent');
numMicrobesL1 = renamevars(numMicrobesL1,"GroupCount","L1");

% Get the number of microbes after filtering
numMicrobesL2 = groupcounts(bootMeanTableFiltered,"Metabolite");
numMicrobesL2 = removevars(numMicrobesL2,'Percent');
numMicrobesL2 = renamevars(numMicrobesL2,"GroupCount","L2");

% Combine tables
numMicrobesLT = outerjoin(numMicrobesL1,numMicrobesL2,'Keys','Metabolite','MergeKeys',true,'Type','left');

% Rename reactions to metabolite names
numMicrobesLT.Metabolite = renameAdrcVmhToMetName(numMicrobesLT.Metabolite);

% Get the number of microbes after lasso regression
filterFun = @(x) x(x.selection_frequency>paramPruning.minFreq,:);
numMicrobesL3 = groupcounts(filterFun(lassoResTab),"Reaction");
numMicrobesL3 = removevars(numMicrobesL3,'Percent');
numMicrobesL3 = renamevars(numMicrobesL3,{'Reaction','GroupCount'},{'Metabolite','L3'});

% Combine tables
numMicrobesLT = outerjoin(numMicrobesLT,numMicrobesL3,'Keys','Metabolite','MergeKeys',true,'Type','left');

% Get the final number of microbes
numMicrobesL4 = table(microbeMetPresence.Properties.VariableNames',sum(table2array(microbeMetPresence))','VariableNames',{'Metabolite','L4'});

% Combine tables
numMicrobesLT = outerjoin(numMicrobesLT,numMicrobesL4,'Keys','Metabolite','MergeKeys',true,'Type','left');
end