% close all
clear all

cd

%%%%% HEALTHY GROUP
% datapath = 'D:\BIOSIGNAL_PROCESSING\MEDICAL IMAGING PROJECT\PilotStudy_2015\RPM and MEMS data (HealthySubjects)\MEMS/'
% datapath_rpm = 'D:\BIOSIGNAL_PROCESSING\MEDICAL IMAGING PROJECT\PilotStudy_2015\RPM and MEMS data (HealthySubjects)\RPM/'


% %%%%%%%%%%%%%%%% Patient group
datapath = 'D:\BIOSIGNAL_PROCESSING\MEDICAL IMAGING PROJECT\MIC2018 project\Micromachines project\MINMOTION_NEWPROJECT\MEMS DATA/'
datapath_rpm = 'D:\BIOSIGNAL_PROCESSING\MEDICAL IMAGING PROJECT\MIC2018 project\Micromachines project\MINMOTION_NEWPROJECT\RPMDATA_SYNC/'

filepattern=fullfile(datapath,'*.txt')
txtfiles=dir(filepattern)

filepattern_rpm=fullfile(datapath_rpm,'*.mat')
rpmfiles=dir(filepattern_rpm)
locs_min_mems=[];locs_min_rpm=[];locs_max_mems=[];locs_max_rpm=[];
MCG_HR=[];MCG_STI=[];MCG_CyCPerc=[];HR_ECG=[];
offset_peaks=[];offset_peaks_std=[];
for k=1:length(txtfiles)
    %% LOAD RPM 
    
    load([datapath_rpm   rpmfiles(k).name]);
    amp=RPMDataAfterSync.amp;
    phase=RPMDataAfterSync.phase;
    timestamp=RPMDataAfterSync.time;
    validflag=RPMDataAfterSync.valid;
    ttlin=RPMDataAfterSync.tin;
    mark=RPMDataAfterSync.mark;
    ttlout=RPMDataAfterSync.tout;
    RPM=(amp-mean(amp))/std(amp);
    RPM=fft_filter(-amp,25,0.1,12);

    %% LOAD MEMS
    data=load2_txt([datapath   txtfiles(k).name])
    fs=800;
    [ RespirationCurves_imu ] = respiration_extractor_from_imu( data, fs );
   
    
    
    
    %% downsampling ADR AND GDR
    ADR_accX=downsample(RespirationCurves_imu.ADRx,32);
    ADR_accY=downsample(RespirationCurves_imu.ADRy,32);
    GDR_gyrX=downsample(RespirationCurves_imu.GDRx,32);
    GDR_gyrY=downsample(RespirationCurves_imu.GDRy,32);
    GDR_gyrz=downsample(RespirationCurves_imu.GDRz,32);
    Resp_PCAFusion=downsample(RespirationCurves_imu.PCA,32);
    samples=1:numel(RPM);
    trpm=samples/25;
    
    
    %% Synching MEMS derived Respiration with RPM
%     RPM=(RPM-mean(RPM))/std(RPM);
   
    s7=Resp_PCAFusion;
    s6=GDR_gyrz;
    s5=GDR_gyrY;
    s4=GDR_gyrX;
    s3=ADR_accY;
    s2=ADR_accX;
    s1=RPM;

    %%%%check the polarity of rpm and pca sigals before synchronizing
    [RPM_syn,PCAFusion_syn] = alignsignals(s1,s7);
    jitter_pca=floor(min(length(PCAFusion_syn),length(RPM_syn)));
if corr(RPM_syn(1:jitter_pca),PCAFusion_syn(1:jitter_pca),'Type','Pearson')>0
    shimmer_pca=find(diff(PCAFusion_syn)>0);
     if numel(s2)>(numel(PCAFusion_syn)-shimmer_pca(1))
        diff_pad=abs(numel(s2)-(numel(PCAFusion_syn)-shimmer_pca(1)));
        zero_padding=zeros(1,shimmer_pca(1)-diff_pad);
     else
            zero_padding=zeros(1,shimmer_pca(1));
     end
else
    [RPM_syn,PCAFusion_syn] = alignsignals(s1,-s7);
    jitter_pca=floor(min(length(PCAFusion_syn),length(RPM_syn)));
    shimmer_pca=find(diff(PCAFusion_syn)>0);
    if numel(s2)>(numel(PCAFusion_syn)-shimmer_pca(1))
        diff_pad=abs(numel(s2)-(numel(PCAFusion_syn)-shimmer_pca(1)));
        zero_padding=zeros(1,shimmer_pca(1)-diff_pad);
    else
        zero_padding=zeros(1,shimmer_pca(1));
    end
end

    ADR_accX_syn=[zero_padding';  (s2)];
    ADR_accY_syn=[zero_padding';  (s3)];
    GDR_gyrX_syn=[zero_padding';  (s4)];
    GDR_gyrY_syn=[zero_padding';  (s5)];
    GDR_gyrZ_syn=[zero_padding';  (s6)];

    PCAFusion_syn=zscore(PCAFusion_syn);
    RPM_syn=zscore(RPM_syn);
    
    RPM_flag_cor=validflag(shimmer_pca(1):jitter_pca);
    
    RPM_sync_cor=RPM_syn(shimmer_pca(1):jitter_pca);
    ADRx_sync_cor=ADR_accX_syn(shimmer_pca(1):jitter_pca);
    ADRy_sync_cor=ADR_accY_syn(shimmer_pca(1):jitter_pca);
    GDRx_sync_cor=GDR_gyrX_syn(shimmer_pca(1):jitter_pca);
    GDRy_sync_cor=GDR_gyrY_syn(shimmer_pca(1):jitter_pca);
    GDRz_sync_cor=GDR_gyrZ_syn(shimmer_pca(1):jitter_pca);

    PCADR_sync_cor=PCAFusion_syn(shimmer_pca(1):jitter_pca);
   
   %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT Sychnced signals
    
    figure
    plot(RPM_syn,'LineWidth',2,'Linestyle','--')
    hold on
    plot(ADR_accX_syn,'LineWidth',1.5,'Linestyle','-')
    hold on
    plot(ADR_accY_syn,'LineWidth',1.5,'Linestyle','-')
    hold on
    plot(GDR_gyrX_syn,'LineWidth',1.5,'Linestyle','-')
    hold on
    plot(GDR_gyrY_syn,'LineWidth',1.5,'Linestyle','-')
    hold on
%     plot(GDR_gyrZ_syn,'LineWidth',1.5,'Linestyle','-')
    hold on
    plot(PCAFusion_syn,'LineWidth',1.5,'Linestyle',':')
    hold on
    
    ylabel('angular Displacement(rad)')
    
    legend('RPM','ADR (acc_X)','ADR (acc_Y)','GDR (gyro_X)','GDR (gyro_Y)',' PCA (acc+gyro)');
    

    %%  Correlation Analysis => rejection mode OFF
    
    [Sim_Score_ADRx(k,:),Pval_ADRx(k,:)]=corr(ADRx_sync_cor,RPM_sync_cor,'Type','Pearson') 
    [Sim_Score_ADRy(k,:),Pval_ADRy(k,:)]=corr(ADRy_sync_cor,RPM_sync_cor,'Type','Pearson')
    
    [Sim_Score_GDRx(k,:),Pval_GDRx(k,:)]=corr(GDRx_sync_cor,RPM_sync_cor,'Type','Pearson')
    
    [Sim_Score_GDRy(k,:),Pval_GDRy(k,:)]=corr(GDRy_sync_cor,RPM_sync_cor,'Type','Pearson')
    
    [Sim_Score_GDRz(k,:),Pval_pca(k,:)]=corr(GDRz_sync_cor,RPM_sync_cor,'Type','Pearson')

    [Sim_Score_pca(k,:),Pval_pca(k,:)]=corr(PCADR_sync_cor,RPM_sync_cor,'Type','Pearson')
    
        %% REJECTION MODE ON
%     [Sim_Score_ADRx(k,:),Pval_ADRx(k,:)]=corr(ADRx_sync_cor(BeamOn),RPM_sync_cor(BeamOn),'Type','Pearson')
%     
%     [Sim_Score_ADRy(k,:),Pval_ADRy(k,:)]=corr(ADRy_sync_cor(BeamOn),RPM_sync_cor(BeamOn),'Type','Pearson')
%     
%     [Sim_Score_GDRx(k,:),Pval_GDRx(k,:)]=corr(GDRx_sync_cor(BeamOn),RPM_sync_cor(BeamOn),'Type','Pearson')
%     
%     [Sim_Score_GDRy(k,:),Pval_GDRy(k,:)]=corr(GDRy_sync_cor(BeamOn),RPM_sync_cor(BeamOn),'Type','Pearson')
%     
%     [Sim_Score_pca(k,:),Pval_pca(k,:)]=corr(PCADR_sync_cor(BeamOn),RPM_sync_cor(BeamOn),'Type','Pearson')
% 
%     if corr(PCADR_sync_cor,RPM_sync_cor,'Type','Pearson')
%     [Sim_Score_MO(k,:),Pval_pca(k,:)]=corr(MEC_OPT,RPM_sync_cor,'Type','Pearson')
%     else
%     [Sim_Score_MO(k,:),Pval_pca(k,:)]=corr(-MEC_OPT,RPM_sync_cor,'Type','Pearson')
%     end
%% selecting the best candidate signals for respiratory gating based on the correlation values
    corr_vec=[Sim_Score_ADRx(k,:),Sim_Score_ADRy(k,:),Sim_Score_GDRx(k,:),Sim_Score_GDRy(k,:),Sim_Score_GDRz(k,:),Sim_Score_pca(k,:)];
    resp_cor_mat=[ADRx_sync_cor,ADRy_sync_cor,GDRx_sync_cor,GDRy_sync_cor,GDRz_sync_cor,PCADR_sync_cor];
    resp_org_mat=[ADR_accX_syn,ADR_accY_syn,GDR_gyrX_syn,GDR_gyrY_syn,GDR_gyrZ_syn,PCAFusion_syn];

    sign_val=sign(corr_vec);
    sign_mat=repmat(sign_val,[size(resp_cor_mat,1),1]);
    resp_cor_mat=sign_mat.*resp_cor_mat;
    [best_resp,axis]=max(abs(corr_vec));
     
     MEMS_resp_candidate=resp_cor_mat(:,axis);
     MEMS_resp_candidate_org=sign_val(axis)*resp_org_mat(:,axis);

  %% Saving respiration signals (Synchronized with RPM)
     RespirationCurves=struct('RPM',-amp,'ADRx',sign_val(1)*resp_org_mat(:,1),'ADRy',sign_val(2)*resp_org_mat(:,3),...
    'GDRx',sign_val(3)*resp_org_mat(:,3),'GDRy',sign_val(4)*resp_org_mat(:,4),'GDRz',sign_val(5)*resp_org_mat(:,5),'PCA',sign_val(6)*resp_org_mat(:,6));
    save(['RespirationCurves_' txtfiles(k).name '.mat'],'RespirationCurves')
      %% Analyzing the rejected cycles
    
        [BeamOff,bb]=find(RPM_flag_cor<0);
     
        [BeamOn,bb]=find(RPM_flag_cor>=0);

        
        figure
        ax(1)=subplot(211);
        plot(RPM_sync_cor,'LineWidth',2)
        hold on
        plot(RPM_flag_cor)
        plot(BeamOff,RPM_sync_cor(BeamOff),'ro','MarkerSize',3,...
    'MarkerEdgeColor','red',...
    'MarkerFaceColor',[1 .6 .6])
        ylabel('Respiratory cycles')
        
        legend('RPM derived respiration','Navigator Pulse','Irregular Respiration (Rejected) ');
        ax(2)=subplot(212);
        plot(MEMS_resp_candidate,'LineWidth',2)
        hold on
        plot(RPM_flag_cor)
        hold on
        plot(BeamOff,MEMS_resp_candidate(BeamOff),'go','MarkerSize',3,...
    'MarkerEdgeColor','red',...
    'MarkerFaceColor',[1 .6 .6])
        legend('PCA derived respiration','Navigator Pulse','Irregular Respiration (Rejected) ');
        ylabel('Respiratory cycles')
        linkaxes([ax(2) ax(1)],'x');
        
        
%% Localizing Maxima and Minima points
     
    [ minloc_mems, maxloc_mems, minloc_rpm, maxloc_rpm ] = find_peaks_rpm_guided( MEMS_resp_candidate,RPM_sync_cor, 1 );
    
    offset_peaks(k,:)= mean(abs(maxloc_rpm-maxloc_mems))/25;
    offset_peaks_std(k,:)=std(abs(maxloc_rpm-maxloc_mems))/25;
    
    count(k,:)=k;
    
  

%% table of correlations    
    T_respiration = table(count,Sim_Score_ADRx,Sim_Score_ADRy,Sim_Score_GDRx,Sim_Score_GDRy,Sim_Score_GDRz,Sim_Score_pca,offset_peaks,offset_peaks_std,...
        'VariableNames',{'Subject','ADRx','ADRy ','GDRx','GDRy','GDRz','FusionPCA','Offset','offsetSTD'})
    
    
    %% Boxplot generation
    correlations=table2array(T_respiration(:,2:end));
    col=@(x)reshape(x,numel(x),1);
    boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),...
        cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)),varargin{:},...
        'Labels',{'ADR_X','ADR_Y','GDR_x','GDR_y','GDR_z','PCA'},'Whisker',1);
    figure
    boxplot2({abs(correlations(:,1)),abs(correlations(:,2)),abs(correlations(:,3)),abs(correlations(:,4)),abs(correlations(:,5)),abs(correlations(:,6))});
 
    xlabel('Respiration Derived Modality (Including Data Fusion)');
    ylabel('Pearsons Correlation Coefficient (r)');
    
    
%     pause
%     continue
     close all
    
    
        %% Gating Data Generation
        [val_rpm,pos_rpm]=setdiff(RPM_syn,RPM_sync_cor);
        [val_mems,pos_mems]=setdiff(MEMS_resp_candidate_org,MEMS_resp_candidate);
        RPM_syn(pos_rpm)=0;
        MEMS_resp_candidate_org(pos_mems)=0;
        fs_cardiac=100;
        fs_resp=25;
        shifter=(shimmer_pca(1)/fs_resp)*fs_cardiac % zeropadding cardiac signals according to respiration signals (fs=100)
       
        read=1;
        if read
          [ LM_Cardiac_5bin, LM_Cardiac_3bin]=peakdetector_cardiacgating_sync(data,fs,shifter,1);
          [ listmode_resp_Phase_5bin ] = PhaseBasedGating(maxloc_mems+shimmer_pca(1),MEMS_resp_candidate_org);
          [ listmode_resp_AMP_5bin ] = MEMS_Amplitude_GATER( MEMS_resp_candidate_org,minloc_mems+shimmer_pca(1),maxloc_mems+shimmer_pca(1));
          [ rpm_listmode_resp_AMP_5bin ] = RPM_Amplitude_GATER(RPM_syn,minloc_rpm+shimmer_pca(1),maxloc_rpm+shimmer_pca(1));
          [ rpm_resp_Phase_5bin ] = PhaseBasedGating(maxloc_rpm+shimmer_pca(1),RPM_syn);

            savdir1 ='D:\BIOSIGNAL_PROCESSING\MEDICAL IMAGING PROJECT\MIC2018 project\Gating Files' ;
           
            
            timeshift=abs(int64(timestamp(1,1)));%% time shift from RPM timestamps
            if timestamp(1,1)>0
                 error(' Time lag is not correct, check the origin')
                 
            end
                
            listmode_resp_AMP_5bin(:,1)=int64(listmode_resp_AMP_5bin(:,1))-timeshift;
            RespGating_listmode_AMP = listmode_resp_AMP_5bin;
            RespLM_AMP=table(RespGating_listmode_AMP);
            
            listmode_resp_Phase_5bin(:,1)=int64(listmode_resp_Phase_5bin(:,1))-timeshift;
            RespGating_listmode_PHASE= listmode_resp_Phase_5bin;
            RespLM_PHASE=table(RespGating_listmode_PHASE);
            
            rpm_listmode_resp_AMP_5bin(:,1)=int64(rpm_listmode_resp_AMP_5bin(:,1))-timeshift;
            RPM_listmode_AMP = rpm_listmode_resp_AMP_5bin;
            RPM_RespLM_AMP=table(RPM_listmode_AMP);
            
            rpm_resp_Phase_5bin(:,1)=int64(rpm_resp_Phase_5bin(:,1))-timeshift;
            RPM_RespGating_listmode_PHASE= rpm_resp_Phase_5bin;
            RPM_RespLM_PHASE=table(RPM_RespGating_listmode_PHASE);

            writetable(RespLM_AMP,['LM_RespGatingLM_5Bins_AMP_' txtfiles(k).name '.txt'],'Delimiter','\t');
            writetable(RespLM_PHASE,['LM_RespGatingLM_5Bins_PHASE_' txtfiles(k).name '.txt'],'Delimiter','\t');
            writetable(LM_Cardiac_5bin,['LMcardiac_MEMS_5Bins_' txtfiles(k).name '.txt'],'Delimiter','\t')
            writetable(LM_Cardiac_3bin,['LMcardiac_MEMS_3Bins_' txtfiles(k).name '.txt'],'Delimiter','\t')
            writetable(RPM_RespLM_AMP,['RPM_LM_RespGatingLM_5Bins_AMP_' txtfiles(k).name '.txt'],'Delimiter','\t');
            writetable(RPM_RespLM_PHASE,['RPM_LM_RespGatingLM_5Bins_PHASE_' txtfiles(k).name '.txt'],'Delimiter','\t');


        end
    
%         pause
        close all
end