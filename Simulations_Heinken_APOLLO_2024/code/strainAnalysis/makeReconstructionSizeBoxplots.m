
% make boxplots with reaction and metabolite numbers for Figure 2b

clear all
rootDir = pwd;

resources={
    'All_statistics_Pasolli_refined.mat','Pasolli'
    'All_statistics_Almeida_refined.mat','Almeida'
    '','APOLLO'
    'All_statistics_AGORA2.mat','AGORA2'
    };

table={};

cnt=1;
for i=1:length(resources)
    if i~=3
    load([rootDir filesep 'data' filesep 'plot_ModelStatistics' filesep resources{i,1}])
    else
        load([rootDir filesep 'data' filesep 'plot_ModelStatistics' filesep resources{1,1}])
        statsTmp=stats;
        load([rootDir filesep 'data' filesep 'plot_ModelStatistics' filesep resources{2,1}])
        stats=vertcat(statsTmp,stats(2:end,:));
    end
    for j=2:length(stats)
        table{cnt,1}='Reactions';
        table{cnt,2}=resources{i,2};
        table{cnt,3}=stats{j,12};
        cnt=cnt+1;
        table{cnt,1}='Metabolites';
        table{cnt,2}=resources{i,2};
        table{cnt,3}=stats{j,13};
        cnt=cnt+1;
        table{cnt,1}='Genes';
        table{cnt,2}=resources{i,2};
        table{cnt,3}=stats{j,14};
        cnt=cnt+1;
    end
end

table=cell2table(table,'VariableNames',{'Feature','Resource','Value'});
table.Resource=categorical(table.Resource,resources(:,2));

figure
boxchart(table.Resource,table.Value,'GroupByColor',table.Feature)
legend('Location','best')
title('Reconstruction properties compared by resource')
set(gca,'FontSize',18)
print([rootDir filesep 'results' filesep 'strains' filesep 'Stats_overview_boxchart'],'-dpng','-r300')
