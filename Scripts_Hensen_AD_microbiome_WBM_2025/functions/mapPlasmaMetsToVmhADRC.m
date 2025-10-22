function [vmhMets,metabolon,vmhTable3] = mapPlasmaMetsToVmhADRC(rxnsToMap,metabolonPath)
% This function hard codes manual mappings from metabolon plasma metabolite
% IDs to VMH ids for the ADRC study.

% Process reactions of interest to VMH metabolites
vmhMets = erase(rxnsToMap,{'DM_','[bc]'})';

% Add VMH information and create table
database = loadVMHDatabase().metabolites;
vmhMets = database(matches(database(:,1),vmhMets),:);
vmhMets(:,2) = lower(vmhMets(:,2));
vmhTable = array2table(vmhMets(:,[1, 2, 6, 7, 8, 9, 11, 13]));
vmhTable.Properties.VariableNames = {'VMHID','CHEMICAL_NAME','KEGG','PUBCHEM','CHEBI','INCHIKEY','HMDB','SUBSYSTEM'};
vmhTable = convertvars(vmhTable,vmhTable.Properties.VariableNames,'string');

% Load metabolite annotations:
metabolon = readtable(metabolonPath,'VariableNamingRule','preserve','Sheet','rowData','ReadRowNames',false);
metabolon = movevars(metabolon,'CHEMICAL_NAME','Before',1);
metabolon = convertvars(metabolon,metabolon.Properties.VariableNames,'string');

% Map metabolon names on VMH metabolites
resource = {'PUBCHEM','KEGG','INCHIKEY','HMDB','CHEMICAL_NAME'};
vmhTable2 = vmhTable;
for i=1:length(resource)
    [~,ia,ib] = intersect(vmhTable.(resource{i}),metabolon.(resource{i}),'stable');
    vmhTable2.CHEM_ID(ia) = metabolon.CHEM_ID(ib);
    vmhTable2.CHEMICAL_NAME(ia) = metabolon.CHEMICAL_NAME(ib);
    vmhTable2.SUB_PATHWAY(ia) = metabolon.SUB_PATHWAY(ib);
end

% Manually repair wrong mapping for 3-dehydro-deoxycholate
vmhTable2.CHEM_ID(matches(vmhTable2.CHEMICAL_NAME,'aconitate [cis or trans]')) = missing; %  <- Wrong mapping on '3-dehydro-deoxycholate'

% Manually map the vmh names on their associated metabolon names

% The following metabolites could be mapped after manual inspection
vmhTable2 = addvars(vmhTable2,repmat("",height(vmhTable2),1) ,'After','CHEMICAL_NAME','NewVariableNames',"alternative_name");
vmhTable2.alternative_name(matches(vmhTable2.CHEMICAL_NAME,"2-hydroxybutyrate")) = "2-hydroxybutyrate/2-hydroxyisobutyrate";
vmhTable2.alternative_name(matches(vmhTable2.CHEMICAL_NAME,"3-methyl-2-oxopentanoate")) = "3-methyl-2-oxovalerate";
vmhTable2.alternative_name(matches(vmhTable2.CHEMICAL_NAME,"adrenaline")) = "NA";
vmhTable2.alternative_name(matches(vmhTable2.CHEMICAL_NAME,"deoxycholic acid-24glucuronide, cda-24g")) = "deoxycholic acid glucuronide";
vmhTable2.alternative_name(matches(vmhTable2.CHEMICAL_NAME,"glycodeoxycholic acid 3-sulfate")) = "glycodeoxycholate 3-sulfate"; 
vmhTable2.alternative_name(matches(vmhTable2.CHEMICAL_NAME,"lithocholic acid 3-sulfate")) = "lithocholate sulfate (1)"; 

% The following metabolites were no perfect matches, but fluxes for the
% alternative metabolites are necessarily identical to the predicted
% metabolites (I think)
vmhTable2.alternative_name(matches(vmhTable2.CHEMICAL_NAME,"3-dehydrocholate")) = "3-dehydrochenodeoxycholate"; % NOTE: Although not the same metabolite, their fluxes are by definition identical.
vmhTable2.alternative_name(matches(vmhTable2.CHEMICAL_NAME,"taurolithocholate")) = "taurolithocholate 3-sulfate"; % Lets think about this
vmhTable2.alternative_name(matches(vmhTable2.CHEMICAL_NAME,"deoxycholic acid 3-sulfate")) = "deoxycholic acid 12-sulfate*"; % Lets think about this. Would there always be identical fluxes?

% The following metabolites could not be mapped after manual inspection
vmhTable2.alternative_name(matches(vmhTable2.CHEMICAL_NAME,"3,4-dihydroxyphenylethyleneglycol")) = "NA";
vmhTable2.alternative_name(matches(vmhTable2.CHEMICAL_NAME,"3alpha,12alpha-dihydroxy-7-oxo-5beta-cholanate")) = "NA";
vmhTable2.alternative_name(matches(vmhTable2.CHEMICAL_NAME,"chenodeoxycholic acid-3glucuronide, cdca-3g")) = "NA"; 
vmhTable2.alternative_name(matches(vmhTable2.CHEMICAL_NAME,"lithocholic acid-3glucuronide, cdca-3g")) = "NA"; 
vmhTable2.alternative_name(matches(vmhTable2.CHEMICAL_NAME,"hyodeoxycholic acid-6glucuronide, hdca-6g")) = "NA"; 
vmhTable2.alternative_name(matches(vmhTable2.CHEMICAL_NAME,"isocholate")) = "NA"; % Closest match isoursodeoxycholate
vmhTable2.alternative_name(matches(vmhTable2.CHEMICAL_NAME,"ursocholate")) = "NA"; % Closest match glycoursodeoxycholic acid sulfate (1)
vmhTable2.alternative_name(matches(vmhTable2.CHEMICAL_NAME,"propionate")) = "NA"; % Closest match: 3-(2-methoxyethoxy)propanoic acid
vmhTable2.alternative_name(matches(vmhTable2.CHEMICAL_NAME,"dopamine")) = "NA"; % Closest match: dopamine 3-O-sulfate
vmhTable2.alternative_name(matches(vmhTable2.CHEMICAL_NAME,"formate")) = "NA"; % Closest match: 3-formylindole
vmhTable2.alternative_name(matches(vmhTable2.CHEMICAL_NAME,"formaldehyde")) = "NA"; % Closest match: ?
vmhTable2.alternative_name(matches(vmhTable2.CHEMICAL_NAME,"homovanillate")) = "NA"; % Closest match: valine 
vmhTable2.alternative_name(matches(vmhTable2.CHEMICAL_NAME,"norepinephrine")) = "NA"; % Closest match:  
vmhTable2.alternative_name(matches(vmhTable2.CHEMICAL_NAME,"metanephrine")) = "NA"; % Closest match:  
vmhTable2.alternative_name(matches(vmhTable2.CHEMICAL_NAME,"s-adenosyl-l-methionine")) = "NA"; % Closest match: methionine sulfone 

% Perform mapping using the manually added alternative metabolite names
[~,ia,ib] = intersect(vmhTable2.('alternative_name'),metabolon.('name'),'stable');
vmhTable2.CHEM_ID(ia) = metabolon.CHEM_ID(ib);
vmhTable2.CHEMICAL_NAME(ia) = metabolon.CHEMICAL_NAME(ib);
vmhTable2.SUB_PATHWAY(ia) = metabolon.SUB_PATHWAY(ib);

% For the metabolites that are not yet mapped, perform a more manual
% investigation:
% vmhMissing = vmhTable2(ismissing(vmhTable2.CHEM_ID),:);
% vmhMissing(matches(vmhMissing.alternative_name,"NA"),:) = [];
% vmhMissing(~matches(vmhMissing.alternative_name,""),:) = [];

% Filter on mapped VMH IDs 
vmhTable3 = vmhTable2(:,{'CHEM_ID','VMHID'});
vmhTable3(ismissing(vmhTable3.CHEM_ID),:)=[];
end