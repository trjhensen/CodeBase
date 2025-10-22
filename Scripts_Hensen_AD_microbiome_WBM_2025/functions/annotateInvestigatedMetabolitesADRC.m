function [metabolitesCheck, metaboliteSummaryStats, metaboliteOntologyStats, capturedSubsystems] = annotateInvestigatedMetabolitesADRC(metPath,mWBMPath)
% function annotateInvestigatedMetabolitesADRC

% input:
% metPath
% mWBMPath

% Goal: Describe the selected metabolites in the manuscript.
% Requirement: Identify the metabolite superclasses, classes, and
% subclasses.
% Requirement 2: Find the inchi key for the metabolites and obtain classyfire identifications at https://cfb.fiehnlab.ucdavis.edu/

% Load a random mWBM model
mWbmInfo = what(mWBMPath);
model = load(fullfile(mWBMPath,mWbmInfo.mat{1}),'rxns','mets','lb','ub');

% Load metabolites to investigate
bloodReactions = readcell(metPath,'Delimiter','\n');
metsToInvestigate = erase(bloodReactions,{'DM_','[bc]'});

% Obtain common name and subsystem:
database = loadVMHDatabase;
[~,~,ib] = intersect(metsToInvestigate, database.metabolites(:,1),'stable');
metAnnotation = string(database.metabolites(ib,[1 2 9 13]));
metTab = array2table(metAnnotation,'VariableNames',{'VMHID','Common name','InChIKey','Subsystem'});

%%%% Add metabolite ontology information 

% Not all metabolites have associated Inchi Keys. Manually add this
% information:

metTab.InChIKey(matches(metTab.VMHID,'3dhchol')) = "OEKUSRBIIZNLHZ-DJDNIQJZSA-N"; % https://www.hmdb.ca/metabolites/HMDB0000502
metTab.InChIKey(matches(metTab.VMHID,'7ocholate')) = "RHCPKKNRWFXMAT-RRWYKFPJSA-N";% https://www.hmdb.ca/metabolites/HMDB0000391
metTab.InChIKey(matches(metTab.VMHID,'cdca3g')) = "GDNGOAUIUTXUES-BWGRGVIUSA-N"; % https://pubchem.ncbi.nlm.nih.gov/compound/Chenodeoxycholic-acid-3-glucuronide
metTab.InChIKey(matches(metTab.VMHID,'isochol')) = "RUDATBOHQWOJDD-JGFDLHJZSA-N"; % https://pubchem.ncbi.nlm.nih.gov/compound/164673#section=IUPAC-Name
metTab.InChIKey(matches(metTab.VMHID,'uchol')) = "BHQCQFFYRZLCQQ-UTLSPDKDSA-N"; % https://pubchem.ncbi.nlm.nih.gov/compound/Ursocholic-acid#section=Names-and-Identifiers

% Note: The remaining missing inchiKeys are all from lithocholate and
% deoxycholate derivatives. The subclasses of deoxy and lithocholate can
% just be appended to these metabolites

% Manual step: Copy inchiKeys to https://cfb.fiehnlab.ucdavis.edu/
% (Accessed at 29 May 2025) and download the classyfire annotations. The
% annotations are then downloaded and manually placed in ~ADRC/inputs/

% Load metabolite ontology table
ontologyPath = which('classyfire_20250529160838.csv');
metOnt = cell2table(readcell(ontologyPath));
metOnt.Properties.VariableNames = metOnt{1,:};
metOnt(1,:) = [];

% Merge metabolite ontology information with the VMH annotation table
metTab = outerjoin(metTab,metOnt,'Keys','InChIKey','MergeKeys',true,'Type','left');

% Remove unneeded variables
metTab = removevars(metTab,{'ClassyFy Status','Parent Level 1','Parent Level 2','Parent Level 3'});

% Fill missing metabolite ontology information 
%metTab(matches(metTab.VMHID,'dchac'),{'Kingdom','Superclass','Class','Subclass'})
metTab(:,'Kingdom') = fillmissing(metTab(:,'Kingdom'),'constant',{'Organic compounds'});
metTab(:,'Superclass') = fillmissing(metTab(:,'Superclass'),'constant',{'Lipids and lipid-like molecules'});
metTab(:,'Class') = fillmissing(metTab(:,'Class'),'constant',{'Steroids and steroid derivatives'});
metTab(:,'Subclass') = fillmissing(metTab(:,'Subclass'),'constant',{'Bile acids, alcohols and derivatives'});

% Prune table
metTab.InChIKey = [];
metabolitesCheck = metTab;

% Create table with metabolite WBM location information

% Find which metabolites are present in the model diet
dietRxnsToCheck = append('Diet_EX_',metAnnotation(:,1),'[d]');
metsInDiet = model.rxns(model.lb<0 & contains(model.rxns,dietRxnsToCheck));
metsInDiet = extractBetween(metsInDiet,'Diet_EX_','[d]');
%dietMetIdx = matches(metAnnotation(:,1),metsInDiet);
dietMetIdx = intersect(metAnnotation(:,1),metsInDiet,'stable');

% Find which of the investigated metabolites are present in the microbiota
% compartment (Think about discriminating between exchanged metabolites and
% only produced metabolites?)
microbiotaMetsToCheck = append(metAnnotation(:,1),'[luM]');
metsInLumen = model.mets(matches(model.mets,microbiotaMetsToCheck));
% metsInLumen(contains(metsInLumen,'pan'))=[];
metsInLumen = erase(metsInLumen,'[luM]');
microMetIdx = intersect(metAnnotation(:,1),metsInLumen,'stable');

% Find metabolites present in just the microbiota lumen
onlyMicrobeMetIdx = setdiff(microMetIdx,dietMetIdx);

% Create table for metabolite presence in the WBMs
% metabolitesCheck = array2table(metAnnotation,'VariableNames',{'VMH ID','Common name','Subsystem'});

% Populate table with metabolite location info
metabolitesCheck = addvars(metabolitesCheck,true(height(metabolitesCheck),1),...
    matches(metabolitesCheck.VMHID,dietMetIdx),...
    matches(metabolitesCheck.VMHID,microMetIdx),...
    matches(metabolitesCheck.VMHID,onlyMicrobeMetIdx),...
    'NewVariableNames',{'Present in blood','Present in diet','Present in microbiota lumen','Non-dietary microbial metabolite'});

% Get summary statistics for each metabolite category
metaboliteSummaryStats = array2table(sum(metabolitesCheck{:,end-3:end})', 'VariableNames',{'Number of metabolites:'});
metaboliteSummaryStats = addvars(metaboliteSummaryStats,metabolitesCheck.Properties.VariableNames(end-3:end)','NewVariableNames',{'Metabolite presence'},'Before',1);

% Find the number of metabolites in each category
metaboliteInfo = metabolitesCheck(:,{'Kingdom','Superclass','Class','Subclass','Subsystem'});
fun = @(x) sortrows(groupcounts(metaboliteInfo(:,x),metaboliteInfo.Properties.VariableNames(x)),'GroupCount','descend');
metaboliteOntologyStats = arrayfun(fun, 1:width(metaboliteInfo),'UniformOutput',false);


% the 116 metabolites.
metabolites = database.metabolites;
metabolites(cellfun(@isempty,metabolites(:,13)),:)=[];
capturedSubsystems = numel(unique(metabolitesCheck.Subsystem)) /  numel(unique(metabolites(:,13))) * 100;


end