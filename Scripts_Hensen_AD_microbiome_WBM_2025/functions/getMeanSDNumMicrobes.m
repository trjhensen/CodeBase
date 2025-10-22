function meanSdNumMicrobesForMets = getMeanSDNumMicrobes(fluxMicrobeData)
% This function obtains the mean average and SD of microbial contributors.
% 
% INPUTS:
% fluxMicrobeData           Cell array with tables with gut microbiome data
%
% OUTPUTS:
% meanSdNumMicrobesForMets  Table with the mean average number of microbes
%                           present for each metabolite.

if matches('ID',fluxMicrobeData{1}.Properties.VariableNames)
    fluxMicrobeData = cellfun(@(x) removevars(x,'ID'), fluxMicrobeData,'UniformOutput',false); 
end
microbeNames = cellfun(@(x) x.Properties.VariableNames(1:end-2), fluxMicrobeData,'UniformOutput',false);  % Get microbe names
metNames = cellfun(@(x) x.Properties.VariableNames(end), fluxMicrobeData,'UniformOutput',true); % Get the metabolite names

% Filter on microbes in the microbiome
getMicrobesFun = @(x,y) table2array(x(:,y'));
subsets = cellfun(getMicrobesFun,fluxMicrobeData, microbeNames,'UniformOutput',false);

% Get the number of microbes per sample
msdFun = @(x) [sum(any(x,1)), mean( sum(x~=0,2) ), std( sum(x~=0,2) )];
meanSdForMets = cellfun(msdFun,subsets,'UniformOutput',false);
meanSdForMets = vertcat(meanSdForMets{:});
meanSdNumMicrobesForMets = array2table(meanSdForMets,'VariableNames',{'Total','Mean','SD'},'RowNames',metNames');
end
