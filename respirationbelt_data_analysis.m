% close all
clear all
cd
datapath = 'D:\BIOSIGNAL_PROCESSING\Physiological Recordings\SCG and Respirationbelt/'

filepattern=fullfile(datapath,'*.txt')
txtfiles=dir(filepattern)



for k=1:length(txtfiles)
  
    %%
    data=read_scg_respirationbelt([datapath   txtfiles(k).name])
    fs=800;
    
    %%%%%%%%%%% for respiratory gating
    st=5*fs;
    sp=numel(data.accZ)-fs*5;
    
    AccX=double(data.accX(st:sp))*(2/32760);
    AccY=double(data.accY(st:sp))*(2/32760);
    AccZ=double(data.accZ(st:sp))*(2/32760);
    
    
    %%%%% For cardiac gating
    
    SCGX=double(data.accX)*(2/32760);
    SCGY=double(data.accY)*(2/32760);
    SCGZ=double(data.accZ)*(2/32760);
    
      
%     hp_A=4;
%     lp_A=40;
%     
%     SCGx=fft_filter(SCGX,fs,hp_A,lp_A);
%     SCGy=fft_filter(SCGY,fs,hp_A,lp_A);
%     SCGz=fft_filter(SCGZ,fs,hp_A,lp_A);
    % %
   
   
    %%
    %%%% respiration extraction from ACCELEROMETER X and Y using ARCTAN
    thetaA=atan(AccY./(sqrt(AccX.^2+AccZ.^2)));
    thetaB=atan(AccX./(sqrt(AccY.^2+AccZ.^2)));
    
    hpf=0.1;
    lpf=2;
    ThetaA=fft_filter(thetaA,fs,hpf,lpf);
    ThetaB=fft_filter(thetaB,fs,hpf,lpf);
    
    %%%%%%%%%%%%%%%%%%%%%%%  respiration extraction from ACCELEROMETER X and Y
    %%%%%%%%%%%%%%%%%%%%%%%  WITHOUT ARCTAN (APPROXIMATION)
    thetaA_aprx=(AccY./(sqrt(AccX.^2+AccZ.^2)));
    thetaB_aprx=(AccX./(sqrt(AccY.^2+AccZ.^2)));
    ThetaA_aprx=fft_filter(thetaA_aprx,fs,hpf,lpf);
    ThetaB_aprx=fft_filter(thetaB_aprx,fs,hpf,lpf);
    %%%%%%%%%%%%%%%%%%%%%%%%% respiration extraction from ACCELEROMETER X and
    %%%%%%%%%%%%%%%%%%%%%%%%% Y using ARCSIN
    alphaA=asin(AccY);
    alphaB=asin(AccX);
    AlphaA=fft_filter(alphaA,fs,hpf,lpf);
    AlphaB=fft_filter(alphaB,fs,hpf,lpf);
    
    AlphaB_aprx=fft_filter(AccX,fs,hpf,lpf);
    AlphaA_aprx=fft_filter(AccY,fs,hpf,lpf);
    
    %%%%% ADR signals from Acc X and Y
%     figure
%     title('Accelerometeric Derived Respiration Rotation')
%     ax(1)=subplot(2,1,1);plot(ThetaA)
%     xlabel('X axis')
%     ylabel('angular Displacement')
%     ax(2)=subplot(2,1,2);plot(ThetaB)
%     ylabel('angular Displacement')
%     xlabel('Y axis')
%     linkaxes([ax(2) ax(1)],'x');

    
    %% Principal component analysis (PCA)
    temp=[ThetaA';ThetaB'];
    [coef,resp_pca,vars]=pca(temp');
    
    resp_pca=fft_filter(resp_pca(:,1),fs,hpf,lpf);
    
    
%     figure
%     title('Accelerometeric Derived Respiration Rotation')
%     ax(1)=subplot(2,1,1);plot(AlphaA)
%     xlabel('Alpha A Y axis')
%     ylabel('angular Displacement')
%     ax(2)=subplot(2,1,2);plot(AlphaB)
%     ylabel('angular Displacement')
%     xlabel('Alpha B X axis')
%     linkaxes([ax(2) ax(1)],'x');
    
    
    %% Median Data Fusion
    
    %% downsampling
    ADR_accX=downsample(ThetaB,32);
    ADR_accY=downsample(ThetaA,32);
    Resp_PCAFusion=downsample(resp_pca(:,1),32);
    
    %%%%%%%%%%%% Synching MEMS derived Respiration with respirationbelt
    r1=downsample(data.respirationBelt(st:sp),2);
    [respiration_filt]=moving_average(r1,200);
    Respbelt=fft_filter(respiration_filt,fs,0.2,5);

%     respiration_belt=medfilt1(respiration_filt,30);
%     figure(33)
%     plot(respiration_belt)

%     Respbelt=fft_filter(data.respirationBelt(st:sp),fs,0.2,3);
    Respbelt=downsample(Respbelt,16);
    samples=1:numel(Respbelt);
    trpm=samples/25;
    

ADR_accX_syn=zscore(ADR_accX);
ADR_accY_syn=zscore(ADR_accY);
PCAFusion_syn=zscore(Resp_PCAFusion);      
Respbelt_syn=zscore(Respbelt);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT Sychnced signals   
    
    figure
    plot(Respbelt_syn,'LineWidth',2,'Linestyle','-')
    hold on
    plot(ADR_accX_syn,'LineWidth',2,'Linestyle','-')
    hold on
    plot(ADR_accY_syn,'LineWidth',2,'Linestyle','-')
    hold on
    plot(PCAFusion_syn,'LineWidth',1.5,'Linestyle','-')
    hold on
    %plot(normalize(MedianFusion_syn,1),'LineWidth',1.5,'Linestyle','-')

  ylabel('angular Displacement(rad)')
    
    legend('Respiration Belt','ADR (acc_X)','ADR (acc_Y)','FUSION (PCA)');
    
    
%     jitter_1=floor(min(length(ADR_accX_syn),length(RPM_syn)));
%     jitter_2=floor(min(length(ADR_accY_syn),length(RPM_syn)));
%     jitter_3=floor(min(length(GDR_gyrX_syn),length(RPM_syn)));
%     jitter_4=floor(min(length(GDR_gyrY_syn),length(RPM_syn)));
%     jitter_5=floor(min(length(PCAFusion_syn),length(RPM_syn)));
% 
%     shimmer_1=find(diff(ADR_accX_syn)>0);
%     shimmer_2=find(diff(ADR_accY_syn)>0);
%     shimmer_3=find(diff(GDR_gyrX_syn)>0);
%     shimmer_4=find(diff(GDR_gyrY_syn)>0);
%     shimmer_5=find(diff(MedianFusion_syn)>0);
%     shimmer_6=find(diff(PCAFusion_syn)>0);
%     pause
%    continue 
  count(k,:)=k;  
%   [Sim_Score_ADRx(k,:),Pval_ADRx(k,:)]=corr(normalize(ADR_accX_syn(shimmer_1(1):jitter_1),1),normalize(RPM_syn(shimmer_1(1):jitter_1),1),'Type','Pearson')
%   
%   [Sim_Score_ADRy(k,:),Pval_ADRy(k,:)]=corr(normalize(ADR_accY_syn(shimmer_2(1):jitter_2),1),normalize(RPM_syn(shimmer_2(1):jitter_2),1),'Type','Pearson')
%   
%   [Sim_Score_GDRx(k,:),Pval_GDRx(k,:)]=corr(normalize(GDR_gyrX_syn(shimmer_3(1):jitter_3),1),normalize(RPM_syn(shimmer_3(1):jitter_3),1),'Type','Pearson')
%   
%   [Sim_Score_GDRy(k,:),Pval_GDRy(k,:)]=corr(normalize(GDR_gyrY_syn(shimmer_4(1):jitter_4),1),normalize(RPM_syn(shimmer_4(1):jitter_4),1),'Type','Pearson')
%   
%   [Sim_Score_median(k,:),Pval_median(k,:)]=corr(normalize(MedianFusion_syn(shimmer_5(1):jitter_5),1),normalize(RPM_syn(shimmer_5(1):jitter_5),1),'Type','Pearson')
%   
%   [Sim_Score_pca(k,:),Pval_pca(k,:)]=corr(normalize(PCAFusion_syn(shimmer_6(1):jitter_6),1),normalize(RPM_syn(shimmer_6(1):jitter_6),1),'Type','Pearson')
%   
%   
%   RespirationCurves=struct('RPM',RPM_syn(shimmer_1(1):jitter_1) ,'ADRx',ADR_accX_syn(shimmer_1(1):jitter_1),'ADRy',ADR_accY_syn(shimmer_2(1):jitter_2),'GDRx',GDR_gyrX_syn(shimmer_3(1):jitter_3),'GDRy',GDR_gyrY_syn(shimmer_4(1):jitter_4),'PCA',PCAFusion_syn(shimmer_6(1):jitter_6));
%                               
                                
shimmer=1
jitter=3000

     [Sim_Score_ADRx_N(k,:),Pval_ADRx(k,:)]=corr(ADR_accX_syn(shimmer:jitter),Respbelt_syn(shimmer:jitter),'Type','Pearson')
       
     [Sim_Score_ADRy_N(k,:),Pval_ADRy(k,:)]=corr(ADR_accY_syn(shimmer:jitter),Respbelt_syn(shimmer:jitter),'Type','Pearson')

     [Sim_Score_pca_N(k,:),Pval_pca(k,:)]=corr(PCAFusion_syn(shimmer:jitter),Respbelt_syn(shimmer:jitter),'Type','Pearson')

     T_respiration_NORM = table(count,Sim_Score_ADRx_N,Sim_Score_ADRy_N,Sim_Score_pca_N,...
    'VariableNames',{'Subject','ADRx','ADRy ','FusionIPCA'})



shimmer=3001;
jitter=4500;

     [Sim_Score_ADRx_S(k,:),Pval_ADRx(k,:)]=corr(ADR_accX_syn(shimmer:jitter),Respbelt_syn(shimmer:jitter),'Type','Pearson')
       
     [Sim_Score_ADRy_S(k,:),Pval_ADRy(k,:)]=corr(ADR_accY_syn(shimmer:jitter),Respbelt_syn(shimmer:jitter),'Type','Pearson')

     [Sim_Score_pca_S(k,:),Pval_pca(k,:)]=corr(PCAFusion_syn(shimmer:jitter),Respbelt_syn(shimmer:jitter),'Type','Pearson')

     T_respiration_SLOW = table(count,Sim_Score_ADRx_S,Sim_Score_ADRy_S,Sim_Score_pca_S,...
    'VariableNames',{'Subject','ADRx','ADRy ','FusionIPCA'})

shimmer=4501;
jitter=numel(ADR_accX_syn);

     [Sim_Score_ADRx_F(k,:),Pval_ADRx(k,:)]=corr(ADR_accX_syn(shimmer:jitter),Respbelt_syn(shimmer:jitter),'Type','Pearson')
       
     [Sim_Score_ADRy_F(k,:),Pval_ADRy(k,:)]=corr(ADR_accY_syn(shimmer:jitter),Respbelt_syn(shimmer:jitter),'Type','Pearson')

     [Sim_Score_pca_F(k,:),Pval_pca(k,:)]=corr(PCAFusion_syn(shimmer:jitter),Respbelt_syn(shimmer:jitter),'Type','Pearson')

     T_respiration_FAST = table(count,Sim_Score_ADRx_F,Sim_Score_ADRy_F,Sim_Score_pca_F,...
    'VariableNames',{'Subject','ADRx','ADRy ','FusionIPCA'})

%     RespirationCurves=struct('RPM',Respbelt_syn(shimmer_pca(1):jitter_pca),'ADRx',ADR_accX_syn(shimmer_pca(1):jitter_pca)*sign(Sim_Score_ADRx(k)),'ADRy',ADR_accY_syn(shimmer_pca(1):jitter_pca)*sign(Sim_Score_ADRy(k,:)),...
%         'GDRx',GDR_gyrX_syn(shimmer_pca(1):jitter_pca)*sign(Sim_Score_GDRx(k)),'GDRy',GDR_gyrY_syn(shimmer_pca(1):jitter_pca)*sign(Sim_Score_GDRy(k)),'PCA',PCAFusion_syn(shimmer_pca(1):jitter_pca)*sign(Sim_Score_pca(k)));                     
%                   save(['RespirationCurves_Subject_ID' num2str(k) '.mat'],'RespirationCurves')
         
                  
                  
                  
   
%  Pvals_respiration = table(count,Pval_ADRx,Pval_ADRy,Pval_GDRx,Pval_GDRy,Pval_pca,...
%     'VariableNames',{'Subject','ADRx','ADRy ','GDRx','GDRy','FusionIPCA'})


            
            
%                [ LM_Cardiac_5bin, LM_Cardiac_3bin]=peakdetector_cardiacgating(data,fs,1);
%    pause
 
  close all
   continue
 
correlations=table2array(T_respiration(:,2:end));
 col=@(x)reshape(x,numel(x),1);
            boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),...
                cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)),varargin{:},...
                'Labels',{'ADR_X','ADR_Y','GDR_x','GDR_y','PCA'},'Whisker',1);
            figure
            boxplot2({abs(correlations(:,1)),abs(correlations(:,2)),abs(correlations(:,3)),abs(correlations(:,4)),abs(correlations(:,5))});
            % title('Rotational Complexity Ratio ')
            %  title(' Similarity Score Between Healthy Subjects  (Serial Measurement)')
            %   title(' Similarity Score Between STEMI Subjects  (Pre PCI vs Post PCI)')
            xlabel('Respiration Derived Modality (Including Data Fusion)');
            ylabel('Pearsons Correlation Coefficient (r)');
 continue
 
 
 
%     %% Gating Data Generation
%    [ minlocs,maxlocs ] = MinMaxLocFinder( ADR_accX_syn,ADR_accY_syn,GDR_gyrX_syn,GDR_gyrY_syn,Resp_KalmanFusion,MedianFusion_syn, Resp_PCAFusion );
%     [ listmode_resp_Phase_5bin ] = PhaseBasedGating(minlocs,MedianFusion_syn);
%     [ listmode_resp_AMP_5bin ] = MEMS_Amplitude_GATER( MedianFusion_syn,minlocs,maxlocs);
%     read=0;
%     if read
%         savdir1 ='D:\BIOSIGNAL_PROCESSING\MEDICAL IMAGING PROJECT\MIC2018 project\Gating Files' ;
%         RespGating_listmode_AMP = int64(listmode_resp_AMP_5bin);
%         RespLM_gyroBindataT=table(RespGating_listmode_AMP);
%         RespGating_listmode_PHASE= int64(listmode_resp_Phase_5bin);
%         RespLM_PHASE=table(RespGating_listmode_PHASE);
% %         writetable(RespLM_gyroBindataT,'LM_RespGatingLM_5Bins_withoutAMP_01062016.txt','Delimiter','\t');
% %         writetable(RespLM_PHASE,'LM_RespGatingLM_5Bins_PHASE_01062016.txt','Delimiter','\t');
%         writetable(RespLM_gyroBindataT,'LM_RespGatingLM_5Bins_withoutAMP_23032016.txt','Delimiter','\t');
%         writetable(RespLM_PHASE,'LM_RespGatingLM_5Bins_PHASE_23032016.txt','Delimiter','\t');
%     end
%     
%     
%     %% cardiac gating Generation
%     
%      [ LM_Cardiac_5bin, LM_Cardiac_3bin]=peakdetector_cardiacgating(data,fs,1);
%     read=0
%     if read
%         
% %         savdir1 ='D:\BIOSIGNAL_PROCESSING\MEDICAL IMAGING PROJECT\MIC2018 project\Gating Files' ;
%         %  gyro_gatingbins=(listmode_data);
% %           save(fullfile(savdir1,sprintf('GyroGateBins_sub_%02d.txt',2)),'listmode_data','-ascii','-tabs');
%         
%         writetable(LM_Cardiac_5bin,'LMcardiac_MEMS_5Bins_01062016.txt','Delimiter','\t')
%         writetable(LM_Cardiac_3bin,'LMcardiac_MEMS_3Bins_01062016.txt','Delimiter','\t')
%     end
%     
%     %  [ LM_CardiacData ] = cardiac_gate_generator( signal, fd_points,gr )
%     
%     % [peaks_scg,fd_gyro]=peakdetector_dualgating(data,fs,detrend(tmptrigy),0,1);
%     
%     pause
%     close all
end