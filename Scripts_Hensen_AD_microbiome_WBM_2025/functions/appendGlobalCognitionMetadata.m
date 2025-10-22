function metadataCog = appendGlobalCognitionMetadata(metadata, metadataOutputFolder,createFig)

% First, set all invalid scores to nans
invalidScores = [88, 95, 96, 97, 98];
numericMD = metadata{:,vartype('numeric')};
numericMD(ismember(numericMD, invalidScores)) = NaN; 
metadata{:,vartype('numeric')} = numericMD;

% Extract cognitive variables of interest and also define the response validity scores (RESPVAL)
varsToUse = {'ID','MOCATOTS','CRAFTVRS','CRAFTURS','CRAFTDVR','CRAFTDRE','DIGFORCT','DIGBACCT','DIGFORSL','DIGBACLS','MINTTOTW','RESPVAL'};
cogData = metadata(:,varsToUse);

% Set all cognitive values where the response is likely not valid, i.e.,
% RESPVAL>1 to nan. 
cogData{cogData.RESPVAL>1,2:end-1} = missing;
cogData.RESPVAL=[]; % Remove variable

% Remove all samples with any nan values
nanRows = any(isnan(cogData{:,vartype('numeric')}),2);
cogData(nanRows,:)=[];

% Generate a univariate cognitive score by z-transforming the individuals
% scores and then performing dimensionality reduction using robust pca.

cogScores = normalize(cogData{:,vartype('numeric')});

% Perform robust PCA
[~,PCscores,~,~,explained,~] = pca(cogScores, 'Centered', false, 'Algorithm', 'als');

% After performing PCA, extract PC1 and merge the data with the metadata
gScore = table(cogData.ID,PCscores(:,1),'VariableNames',{'ID','G'});

% Add global cognition to metadata
metadataCog = outerjoin(metadata,gScore,'Keys','ID','MergeKeys',true,'Type','left');


% Describe cognition data
if createFig == true
    % Correlate cognitive scores
    varNames = cogData.Properties.VariableNames(2:end);
    [RHO, ~] = corr(cogScores);
    RHO = round(RHO,2);
    
    f1 = figure; 
    h = heatmap(RHO);
    ax = gca;
    ax.YDisplayLabels = varNames';
    ax.XDisplayLabels = varNames';
    ax.Title = 'Pearson correlations between cognitive test scores';
    h.ColorLimits = [0 1];
    
    % Save figure
    path = fullfile(metadataOutputFolder,'cognitive_test_correlations.png');
    exportgraphics(f1,path,'Resolution',600)
    
    % Visualise cognitive scores from PCA
    f2 = figure;
    scatter(PCscores(:,1), PCscores(:,2));
    % Axis titles
    xlabel(strcat("PC1 explained variance: ",string(explained(1)),"% "))
    ylabel(strcat("PC2 explained variance: ",string(explained(2)),"% "))
    title('PCA of cognitive scores')
    
    % Save figure
    path = fullfile(metadataOutputFolder,'cognitive_test_PCA.png');
    exportgraphics(f2,path,'Resolution',600)
    
    % Close figures
    close all
end

end
