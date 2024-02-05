% close all
clear all

cd

%%%%% HEALTHY GROUP
% datapath = 'D:\BIOSIGNAL_PROCESSING\MEDICAL IMAGING PROJECT\PilotStudy_2015\RPM and MEMS data (HealthySubjects)\MEMS/'
% datapath_rpm = 'D:\BIOSIGNAL_PROCESSING\MEDICAL IMAGING PROJECT\PilotStudy_2015\RPM and MEMS data (HealthySubjects)\RPM/'


% %%%%%%%%%%%%%%%% Patient group
datapath = 'D:\BIOSIGNAL_PROCESSING\MEDICAL IMAGING PROJECT\MinMotion_Micromachine paper project\RPM and MEMS data (HealthySubjects)\PLAQUA Study patients\MEMS/'
datapath_rpm = 'D:\BIOSIGNAL_PROCESSING\MEDICAL IMAGING PROJECT\MinMotion_Micromachine paper project\RPM and MEMS data (HealthySubjects)\PLAQUA Study patients\RPM/'

filepattern=fullfile(datapath,'*.txt')
txtfiles=dir(filepattern)

filepattern_rpm=fullfile(datapath_rpm,'*.vxp')
vpxfiles=dir(filepattern_rpm)
locs_min_mems=[];locs_min_rpm=[];locs_max_mems=[];locs_max_rpm=[];
MCG_HR=[];MCG_STI=[];MCG_CyCPerc=[];HR_ECG=[];
offset_peaks=[];offset_peaks_std=[];
for k=1:length(txtfiles)
    %% loading RPM Signals
    
    [Data_layoutamplitude,phase,timestamp,validflag,ttlin,mark,ttlout] = importRPM([datapath_rpm   vpxfiles(k).name]);
    
    RPM=(Data_layoutamplitude-mean(Data_layoutamplitude))/std(Data_layoutamplitude);
    RPM=fft_filter(Data_layoutamplitude,25,0.1,12);
    %  RPM= resample(RPM,32,1);
    %  GDR=downsample(RotX,32);
    %  ADR=downsample(ThetaB,32);
    %%
    data=load2_txt([datapath   txtfiles(k).name])
    fs=800;
    
    %%%%%%%%%%% for respiratory gating
    st=10*fs;
    sp=numel(data.gyroY)-fs*10;
    gyroY=double(data.gyroY(st:sp))*(250/32767);
    gyroX=double(data.gyroX(st:sp))*(250/32767);
    gyroZ=double(data.gyroZ(st:sp))*(250/32767);
    
    AccX=double(data.accX(st:sp))*(2/32760);
    AccY=double(data.accY(st:sp))*(2/32760);
    AccZ=double(data.accZ(st:sp))*(2/32760);
    Hp=0.01;
    GyroX=baseline_removal_imu(gyroX,Hp,fs);
    GyroY=baseline_removal_imu(gyroY,Hp,fs);
    GyroZ=baseline_removal_imu(gyroZ,Hp,fs);
    
    
    %%%%%  cardiac gating
    
     [ LM_Cardiac_5bin, LM_Cardiac_3bin,MCG_HR(k,:),MCG_STI(k,:),MCG_CyCPerc(k,:),HR_ECG(k,:)]=peakdetector_cardiacgating(data,fs,1);

    pause
%     continue
    %%
    %%%% respiration extraction from ACCELEROMETER X and Y using ARCTAN
    thetaA=atan(AccX./(sqrt(AccY.^2+AccZ.^2)));
    thetaB=atan(AccY./(sqrt(AccX.^2+AccZ.^2)));
    
    hpf=0.1;
    lpf=1;
    ThetaA=fft_filter(thetaA,fs,hpf,lpf);
    ThetaB=fft_filter(thetaB,fs,hpf,lpf);
    
    %%%%%%%%%%%%%%%%%%%%%%%  respiration extraction from ACCELEROMETER X and Y
    %%%%%%%%%%%%%%%%%%%%%%%  WITHOUT ARCTAN (APPROXIMATION)
    thetaA_aprx=(AccX./(sqrt(AccY.^2+AccZ.^2)));
    thetaB_aprx=(AccY./(sqrt(AccX.^2+AccZ.^2)));
    ThetaA_aprx=fft_filter(thetaA_aprx,fs,hpf,lpf);
    ThetaB_aprx=fft_filter(thetaB_aprx,fs,hpf,lpf);
    %%%%%%%%%%%%%%%%%%%%%%%%% respiration extraction from ACCELEROMETER X and
    %%%%%%%%%%%%%%%%%%%%%%%%% Y using ARCSIN
    alphaA=asin(AccX);
    alphaB=asin(AccY);
    AlphaA_aprx=fft_filter(alphaA,fs,hpf,lpf);
    AlphaB_aprx=fft_filter(alphaB,fs,hpf,lpf);
    
    AlphaA=fft_filter(AccX,fs,hpf,lpf);
    AlphaB=fft_filter(AccY,fs,hpf,lpf);
    
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% respiration extraction from GYROSCOPE X and
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Y by INTEGRATING (convering angular velocity to
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% angular displacement)
    samples=1:numel(AccX);
    t=samples/fs;
    
    rotX=cumtrapz(t,GyroX);
    RotX=fft_filter(rotX*0.017,fs,hpf,lpf);
    
    rotY=cumtrapz(t,GyroY);
    RotY=fft_filter(rotY*0.017,fs,hpf,lpf);
    
    
    rotZ=cumtrapz(t,GyroZ);
    RotZ=fft_filter(rotZ*0.017,fs,hpf,lpf);
    
    
    %%%%%%%%%
    Gama_zaxis=fft_filter(AccZ,fs,hpf,lpf);
%     [adrY_corr,~]=corr(Gama_zaxis,ThetaB,'Type','Pearson')
%      [adrX_corr,~]=corr(Gama_zaxis,ThetaA,'Type','Pearson')
%         [gdrX_corr,~]=corr(Gama_zaxis,RotX,'Type','Pearson')
%             [gdrY_corr,~]=corr(Gama_zaxis,RotY,'Type','Pearson')
%                 [gdrZ_corr,~]=corr(Gama_zaxis,RotZ,'Type','Pearson')
%     corr_val=[adrX_corr,adrY_corr,gdrX_corr,gdrY_corr,gdrZ_corr];

   [adrY_corr,~]=corr(Gama_zaxis,ThetaB,'Type','Pearson')
     [adrX_corr,~]=corr(Gama_zaxis,ThetaA,'Type','Pearson')
        [gdrX_corr,~]=corr(Gama_zaxis,RotX,'Type','Pearson')
            [gdrY_corr,~]=corr(Gama_zaxis,RotY,'Type','Pearson')
%                 [gdrZ_corr,~]=corr(Gama_zaxis,RotZ,'Type','Pearson')
    
    corr_val=[adrX_corr,adrY_corr,gdrX_corr,gdrY_corr];
    sign_val=sign(corr_val);
        resp_signals=[ThetaA,ThetaB,RotX,RotY];
        sign_mat=repmat(sign_val,[size(resp_signals,1),1]);
        new_resp_sig=sign_mat.*resp_signals;
        
        ADRx=zscore(new_resp_sig(:,1));
        ADRy=zscore(new_resp_sig(:,2));
        GDRx=zscore(new_resp_sig(:,3));
        GDRy=zscore(new_resp_sig(:,4));
        GDRz=zscore(RotZ);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ADRx_ap1=zscore(ThetaA_aprx);
        ADRy_ap1=zscore(ThetaB_aprx);
        
        ADRx_ap2=zscore(AlphaA_aprx);
        ADRy_ap2=zscore(AlphaB_aprx);
%        figure
%       ax(1)= subplot(411);
%           plot(t,GDRx,'LineWidth',2,'Linestyle','-')
%         hold on
%         plot(t,GDRy,'LineWidth',2,'Linestyle','-')
%         hold on
%         plot(t,-GDRz,'LineWidth',2,'Linestyle','-')
%         
%        ax(2)=subplot(412);
%         plot(t,ADRx,'LineWidth',2,'Linestyle','-')
%         hold on
%         plot(t,ADRy,'LineWidth',2,'Linestyle','-')
%        ax(3)= subplot(413);
%         plot(t,ADRx_ap1,'LineWidth',2,'Linestyle','-')
%         hold on
%         plot(t,ADRy_ap1,'LineWidth',2,'Linestyle','-')
%         ax(4)=subplot(414);
%         plot(t,ADRx_ap2,'LineWidth',2,'Linestyle','-')
%         hold on
%         plot(t,ADRy_ap2,'LineWidth',2,'Linestyle','-')
%         linkaxes([ax(4) ax(3) ax(2) ax(1)],'x');

        
    %%%%% ADR signals from Acc X and Y
    figure
    title('Accelerometeric Derived Respiration Rotation')
    ax(1)=subplot(2,1,1);plot(ADRx)
    xlabel('X axis')
    ylabel('angular Displacement')
    ax(2)=subplot(2,1,2);plot(ADRy)
    ylabel('angular Displacement')
    xlabel('Y axis')
    linkaxes([ax(2) ax(1)],'x');
    
    
    
    
    
    %% Principal component analysis (PCA)
    temp=[(ADRx');(ADRy');(GDRx');(GDRy');(GDRz')];
%         temp=[RotZ';(RotY');(RotX');(ThetaA');(ThetaB')];

    [coef,resp_pca,vars]=pca(temp');
%      resp_pca=fft_filter(resp_pca(:,1),fs,hpf,lpf);
    
        resp_pca=resp_pca(:,1);

    
    %% signal analysis
    figure
    plot(t,ADRx,'LineWidth',2,'Linestyle','-')
    hold on
    plot(t,ADRy,'LineWidth',2,'Linestyle','-')
    hold on
    plot(t,GDRx,'LineWidth',2,'Linestyle','-')
    hold on
    plot(t,GDRy,'LineWidth',2,'Linestyle','-')
    hold on
    plot(t,GDRz,'LineWidth',2,'Linestyle','-')
    plot(t,resp_pca,'LineWidth',1.5,'Linestyle','-')
%     hold on
%     plot(t,fused_four,'LineWidth',1.5,'Linestyle','-')
    %     plot(t,R,'LineWidth',1.5,'Linestyle','-')
    ylabel('angular Displacement(rad)')
    
    legend('ADR (acc_Y)','ADR (acc_X)','GDR (gyro_X)','GDR (gyro_Y)', 'GDR (gyro_Z)','PCA');
    
    
    
    %% downsampling
    ADR_accX=downsample(ADRx,32);
    ADR_accY=downsample(ADRy,32);
    GDR_gyrX=downsample(GDRx,32);
    GDR_gyrY=downsample(GDRy,32);
    GDR_gyrz=downsample(GDRz,32);
    Resp_PCAFusion=downsample(resp_pca,32);
%     Resp_medianFusion=downsample(fused_four,32);
    samples=1:numel(RPM);
    trpm=samples/25;
    %%%%%%%%%%%% Synching MEMS derived Respiration with RPM
    RPM=-(RPM-mean(RPM))/std(RPM);
   
    s7=Resp_PCAFusion;
    s6=GDR_gyrz;
    s5=GDR_gyrY;
    s4=GDR_gyrX;
    s3=ADR_accY;
    s2=ADR_accX;
    s1=RPM;
    %     [RPM_syn,ADR_accX_syn] = alignsignals(s1,s2);
    %     [x13,ADR_accY_syn] = alignsignals(s1,s3);
    %     [x14,GDR_gyrX_syn] = alignsignals(s1,s4);
    %     [x15,GDR_gyrY_syn] = alignsignals(s1,s5);
    %     [x16,MedianFusion_syn] = alignsignals(s1,s6);
    
    [RPM_syn,PCAFusion_syn] = alignsignals(-s1,s7);
    shimmer_pca=find(diff(PCAFusion_syn)>0);
    jitter_pca=floor(min(length(PCAFusion_syn),length(RPM_syn)));
    zero_padding=zeros(1,shimmer_pca(1));
    ADR_accX_syn=[zero_padding';  zscore(s2)];
    ADR_accY_syn=[zero_padding';  zscore(s3)];
    GDR_gyrX_syn=[zero_padding';  zscore(s4)];
    GDR_gyrY_syn=[zero_padding';  zscore(s5)];
    GDR_gyrZ_syn=[zero_padding';  zscore(s6)];

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
    MEC_OPT=PCADR_sync_cor+RPM_sync_cor;
%     
    
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
        plot(PCADR_sync_cor,'LineWidth',2)
        hold on
        plot(RPM_flag_cor)
        hold on
        plot(BeamOff,PCADR_sync_cor(BeamOff),'go','MarkerSize',3,...
    'MarkerEdgeColor','red',...
    'MarkerFaceColor',[1 .6 .6])
        legend('PCA derived respiration','Navigator Pulse','Irregular Respiration (Rejected) ');
        ylabel('Respiratory cycles')
        linkaxes([ax(2) ax(1)],'x');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT Sychnced signals
    
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
    
    
    [ minloc_pca, maxloc_pca, minloc_rpm, maxloc_rpm ] = find_peaks_rpm_guided( PCADR_sync_cor,RPM_sync_cor, 1 )
    offset_peaks(k,:)= mean(abs(maxloc_rpm-maxloc_pca))/25;
    offset_peaks_std(k,:)=std(abs(maxloc_rpm-maxloc_pca))/25;
    
    count(k,:)=k;
    
    %% rejection mode OFF
    [Sim_Score_ADRx(k,:),Pval_ADRx(k,:)]=corr(ADRx_sync_cor,RPM_sync_cor,'Type','Pearson') 
    [Sim_Score_ADRy(k,:),Pval_ADRy(k,:)]=corr(ADRy_sync_cor,RPM_sync_cor,'Type','Pearson')
    
    [Sim_Score_GDRx(k,:),Pval_GDRx(k,:)]=corr(GDRx_sync_cor,RPM_sync_cor,'Type','Pearson')
    
    [Sim_Score_GDRy(k,:),Pval_GDRy(k,:)]=corr(GDRy_sync_cor,RPM_sync_cor,'Type','Pearson')
    
    [Sim_Score_GDRz(k,:),Pval_pca(k,:)]=corr(GDRz_sync_cor,RPM_sync_cor,'Type','Pearson')

    [Sim_Score_pca(k,:),Pval_pca(k,:)]=corr(PCADR_sync_cor,RPM_sync_cor,'Type','Pearson')
    
    if corr(PCADR_sync_cor,RPM_sync_cor,'Type','Pearson')
    [Sim_Score_MO(k,:),Pval_pca(k,:)]=corr(MEC_OPT,RPM_sync_cor,'Type','Pearson')
    else
    [Sim_Score_MO(k,:),Pval_pca(k,:)]=corr(-MEC_OPT,RPM_sync_cor,'Type','Pearson')
    end

    
    


    
%     maxloc_rpm = ampd(RPM_sync_cor);
%     minloc_rpm = ampd(-RPM_sync_cor);
%     
%        [ minlocs,maxlocs ] = MinMaxLocFinder( ADRx_sync_cor,ADRy_sync_cor,GDRx_sync_cor,GDRy_sync_cor, PCADR_sync_cor );
%        
%        locs_min_mems=[ locs_min_mems, minlocs];
%        locs_min_rpm=[locs_min_rpm, minloc_rpm];
%        locs_max_mems=[ locs_max_mems, maxlocs];
%        locs_max_rpm=[ locs_max_rpm, maxloc_rpm];
% %         [ new_locs_rpm ] = match_resp_locs( locs_max_rpm,locs_max_mems );
% 
% RR_int_rpm=diff(locs_max_rpm)/25;
% RR_rpm_permin=(60./RR_int_rpm);
% 
% [a,b]=find(RR_rpm_permin>45 | RR_rpm_permin<0)
% 
% RR_rpm_permin (b)=[];
% 
% RR_int_mems=diff(locs_max_mems)/25;
% RR_mems_permin=(60./RR_int_mems);
% 
% [a,b]=find(RR_mems_permin>45 | RR_mems_permin<0)
% 
% RR_mems_permin (b)=[];
% 
% x1 = rand(5,1);
% x = [x1;RR_mems_permin'; RR_rpm_permin'];
% g1 = repmat({'First'},5,1);
% 
% g2 = repmat({'MEMS'},numel(RR_mems_permin),1);
% g3 = repmat({'RPM'},numel(RR_rpm_permin),1);
% g = [g1; g2;g3];
% figure;boxplot(x,g);
%  [p,h] = ranksum(RR_mems_permin,RR_rpm_permin)


%         [ listmode_resp_Phase_5bin ] = PhaseBasedGating(minlocs,PCADR_sync_cor);
%         [ listmode_resp_AMP_5bin ] = MEMS_Amplitude_GATER( PCADR_sync_cor,minlocs,maxlocs);
    %     read=0;
%       [Sim_Score_ADRx(k,:),Pval_ADRx(k,:)]=corr(normalize(ADR_accX_syn(shimmer_pca(1):jitter_1),1),normalize(RPM_syn(shimmer_1(1):jitter_1),1),'Type','Pearson')
%     
%       [Sim_Score_ADRy(k,:),Pval_ADRy(k,:)]=corr(normalize(ADR_accY_syn(shimmer_pca(1):jitter_2),1),normalize(RPM_syn(shimmer_2(1):jitter_2),1),'Type','Pearson')
%     
%       [Sim_Score_GDRx(k,:),Pval_GDRx(k,:)]=corr(normalize(GDR_gyrX_syn(shimmer_pca(1):jitter_3),1),normalize(RPM_syn(shimmer_3(1):jitter_3),1),'Type','Pearson')
%     
%       [Sim_Score_GDRy(k,:),Pval_GDRy(k,:)]=corr(normalize(GDR_gyrY_syn(shimmer_pca(1):jitter_4),1),normalize(RPM_syn(shimmer_4(1):jitter_4),1),'Type','Pearson')
%     
%       [Sim_Score_median(k,:),Pval_median(k,:)]=corr(normalize(MedianFusion_syn(shimmer_pca(1):jitter_5),1),normalize(RPM_syn(shimmer_5(1):jitter_5),1),'Type','Pearson')
%     
%       [Sim_Score_pca(k,:),Pval_pca(k,:)]=corr(normalize(PCAFusion_syn(shimmer_pca(1):jitter_6),1),normalize(RPM_syn(shimmer_6(1):jitter_6),1),'Type','Pearson')
%     
    
%       RespirationCurves=struct('RPM',RPM_syn(shimmer_1(1):jitter_1) ,'ADRx',ADR_accX_syn(shimmer_1(1):jitter_1),'ADRy',ADR_accY_syn(shimmer_2(1):jitter_2),'GDRx',GDR_gyrX_syn(shimmer_3(1):jitter_3),'GDRy',GDR_gyrY_syn(shimmer_4(1):jitter_4),'PCA',PCAFusion_syn(shimmer_6(1):jitter_6));
%     
%     
%                 [Sim_Score_ADRx(k,:),Pval_ADRx(k,:)]=corr(ADR_accX_syn(shimmer_pca(1):jitter_pca),RPM_syn(shimmer_pca(1):jitter_pca),'Type','Pearson')
%     
%             [Sim_Score_ADRy(k,:),Pval_ADRy(k,:)]=corr(ADR_accY_syn(shimmer_pca(1):jitter_pca),RPM_syn(shimmer_pca(1):jitter_pca),'Type','Pearson')
%     
%                     [Sim_Score_GDRx(k,:),Pval_GDRx(k,:)]=corr(GDR_gyrX_syn(shimmer_pca(1):jitter_pca),RPM_syn(shimmer_pca(1):jitter_pca),'Type','Pearson')
%     
%                         [Sim_Score_GDRy(k,:),Pval_GDRy(k,:)]=corr(GDR_gyrY_syn(shimmer_pca(1):jitter_pca),RPM_syn(shimmer_pca(1):jitter_pca),'Type','Pearson')
%     
%                                     [Sim_Score_pca(k,:),Pval_pca(k,:)]=corr(PCAFusion_syn(shimmer_pca(1):jitter_pca),RPM_syn(shimmer_pca(1):jitter_pca),'Type','Pearson')
    
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
    
    
%     RespirationCurves=struct('RPM',RPM_syn(shimmer_pca(1):jitter_pca),'ADRx',ADR_accX_syn(shimmer_pca(1):jitter_pca)*sign(Sim_Score_ADRx(k)),'ADRy',ADR_accY_syn(shimmer_pca(1):jitter_pca)*sign(Sim_Score_ADRy(k,:)),...
%         'GDRx',GDR_gyrX_syn(shimmer_pca(1):jitter_pca)*sign(Sim_Score_GDRx(k)),'GDRy',GDR_gyrY_syn(shimmer_pca(1):jitter_pca)*sign(Sim_Score_GDRy(k)),'PCA',PCAFusion_syn(shimmer_pca(1):jitter_pca)*sign(Sim_Score_pca(k)));
%     save(['RespirationCurves_Subject_ID' num2str(k) '.mat'],'RespirationCurves')
%     
%     
    
    T_respiration = table(count,Sim_Score_ADRx,Sim_Score_ADRy,Sim_Score_GDRx,Sim_Score_GDRy,Sim_Score_GDRz,Sim_Score_pca,offset_peaks,offset_peaks_std,...
        'VariableNames',{'Subject','ADRx','ADRy ','GDRx','GDRy','GDRz','FusionPCA','Offset','offsetSTD'})
    
%     Pvals_respiration = table(count,Pval_ADRx,Pval_ADRy,Pval_GDRx,Pval_GDRy,Pval_pca,...
%         'VariableNames',{'Subject','ADRx','ADRy ','GDRx','GDRy','FusionIPCA'})
    
   pause
%     
    close all
    
    
    correlations=table2array(T_respiration(:,2:end));
    col=@(x)reshape(x,numel(x),1);
    boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),...
        cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)),varargin{:},...
        'Labels',{'ADR_X','ADR_Y','GDR_x','GDR_y','GDR_z','PCA'},'Whisker',1);
    figure
    boxplot2({abs(correlations(:,1)),abs(correlations(:,2)),abs(correlations(:,3)),abs(correlations(:,4)),abs(correlations(:,5)),abs(correlations(:,6))});
 
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