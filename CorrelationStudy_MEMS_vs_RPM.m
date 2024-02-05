close all
clear all

datapath = 'D:\Seafile\Minmotion\Algorithm Development and Data Analysis\Mojtaba\MinMotion_Micromachine paper project\Respiration curves/'

filepattern=fullfile(datapath,'*.mat')
matfiles=dir(filepattern)

for k=1:length(matfiles)
    %     try
    baseFileName = matfiles(k).name;
    load([fullfile(datapath, baseFileName)]);   % load the respiration data
    
    
    %% Calculate Pearson Correlation and P-values
    count(k,:)=k;
    [Sim_Score_ADRx(k,:),Pval_ADRx(k,:)]=corr(RespirationCurves.ADRx,RespirationCurves.RPM,'Type','Pearson')
    
    [Sim_Score_ADRy(k,:),Pval_ADRy(k,:)]=corr(RespirationCurves.ADRy,RespirationCurves.RPM,'Type','Pearson')
    
    [Sim_Score_GDRx(k,:),Pval_GDRx(k,:)]=corr(RespirationCurves.GDRx,RespirationCurves.RPM,'Type','Pearson')
    
    [Sim_Score_GDRy(k,:),Pval_GDRy(k,:)]=corr(RespirationCurves.GDRy,RespirationCurves.RPM,'Type','Pearson')
    
    [Sim_Score_pca(k,:),Pval_pca(k,:)]=corr(RespirationCurves.PCA,RespirationCurves.RPM,'Type','Pearson')
    
    %%plot the signals
    figure
    plot(RespirationCurves.RPM,'LineWidth',2,'Linestyle','-')
    hold on
    plot(RespirationCurves.ADRx,'LineWidth',2,'Linestyle','-')
    hold on
    plot(RespirationCurves.ADRy,'LineWidth',2,'Linestyle','-')
    hold on
    plot(RespirationCurves.GDRx,'LineWidth',2,'Linestyle','-')
    hold on
    plot(RespirationCurves.GDRy,'LineWidth',1.5,'Linestyle','-')
    hold on
    plot(RespirationCurves.PCA,'LineWidth',1.5,'Linestyle','-')
    hold on
    
    ylabel('angular Displacement(rad)')
    
    legend('RPM','ADR (acc_X)','ADR (acc_Y)','GDR (gyro_X)','GDR (gyro_Y)','FUSION I (PCA)');
    
end


T_respiration = table(count,Sim_Score_ADRx,Sim_Score_ADRy,Sim_Score_GDRx,Sim_Score_GDRy,Sim_Score_pca,...
    'VariableNames',{'Subject','ADRx','ADRy ','GDRx','GDRy','FusionIPCA'})

Pvals_respiration = table(count,Pval_ADRx,Pval_ADRy,Pval_GDRx,Pval_GDRy,Pval_pca,...
    'VariableNames',{'Subject','ADRx','ADRy ','GDRx','GDRy','FusionIPCA'})


correlations=table2array(T_respiration(:,2:end));
col=@(x)reshape(x,numel(x),1);
boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),...
    cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)),varargin{:},...
    'Labels',{'ADR_X','ADR_Y','GDR_x','GDR_y','PCA'},'Whisker',1);
figure
boxplot2({abs(correlations(:,1)),abs(correlations(:,2)),abs(correlations(:,3)),abs(correlations(:,4)),abs(correlations(:,5))});
title('Correlation Study [MEMS vs. RPM] ')
xlabel('Respiration Derived Modality (Including Data Fusion)');
ylabel('Pearsons Correlation Coefficient (r)');

