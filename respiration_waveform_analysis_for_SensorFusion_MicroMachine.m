close all
clear all
cd
datapath = 'D:\BIOSIGNAL_PROCESSING\MEDICAL IMAGING PROJECT\PilotStudy_2015\RPM and MEMS data (HealthySubjects)\MEMS/'
datapath_rpm = 'D:\BIOSIGNAL_PROCESSING\MEDICAL IMAGING PROJECT\PilotStudy_2015\RPM and MEMS data (HealthySubjects)\RPM/'

% datapath = 'D:\BIOSIGNAL_PROCESSING\MEDICAL IMAGING PROJECT\MIC2018 project\Micromachines project/'
%  datapath='D:\BIOSIGNAL_PROCESSING\MEDICAL IMAGING PROJECT\PilotStudy_2015\RPM and MEMS data (HealthySubjects)/'
%datapath='D:\BIOSIGNAL_PROCESSING\MEDICAL IMAGING PROJECT\PilotStudy_2015\HI BMI subjects\data/'
filepattern=fullfile(datapath,'*.txt')
txtfiles=dir(filepattern)

filepattern_rpm=fullfile(datapath_rpm,'*.vxp')
vpxfiles=dir(filepattern_rpm)


for k=1:length(txtfiles)
      %% loading RPM Signals
    
%     load( 'plaque2_sub2_23032016' )
    
%      [Data_layoutamplitude,phase,timestamp,validflag,ttlin,mark,ttlout] = importRPM('patient_01_06_2016_4DPET.vxp');
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

    
    %%%%% For cardiac gating
    
    SCGX=double(data.accX)*(2/32760);
    SCGY=double(data.accY)*(2/32760);
    SCGZ=double(data.accZ)*(2/32760);
    
    GCGY=double(data.gyroY)*(250/32767);
    GCGX=double(data.gyroX)*(250/32767);
    GCGZ=double(data.gyroZ)*(250/32767);
    
    hp_A=4;
    lp_A=40;
    
    SCGx=fft_filter(SCGX,800,hp_A,lp_A);
    SCGy=fft_filter(SCGY,800,hp_A,lp_A);
    SCGz=fft_filter(SCGZ,800,hp_A,lp_A);
    % %
    hp_G=1;
    lp_G=20;
    GCGx=fft_filter(GCGX,800,hp_G,lp_G);
    GCGy=fft_filter(GCGY,800,hp_G,lp_G);
    GCGz=fft_filter(GCGZ,800,hp_G,lp_G);
    
    
    
    %%  %%%%%%%%%%%%%%%%%%%%%%%%FFT ANALYSIS
    FFTanalysis=0
    if FFTanalysis
        SCG_X=AccX(10000:140000); %%%% TAKE A SAMPLE PORTION OF THE SIGNALS
        SCG_Y=AccY(10000:140000);
        SCG_Z=AccZ(10000:140000);
        
        
        GCG_X=GyroX(10000:140000);
        GCG_Y=GyroY(10000:140000);
        GCG_Z=GyroZ(10000:140000);
        
        m = length(GCG_X);       % original sample length
        n = pow2(nextpow2(m));  % transform length
        y = fft(GCG_X,n);
        f = (0:n-1)*(fs/n);
        power = abs(y).^2/n;
        figure;
        ax(1)=subplot(321);plot(f(1:floor(n/2)),power(1:floor(n/2)))
        title('Spectral Analysis of GCG X axis')
        xlabel('Frequency')
        ylabel('Power')
        m = length(GCG_Y);       % original sample length
        n = pow2(nextpow2(m));  % transform length
        y = fft(GCG_Y,n);
        f = (0:n-1)*(fs/n);
        power = abs(y).^2/n;
        ax(2)=subplot(323);plot(f(1:floor(n/2)),power(1:floor(n/2)))
        title('Spectral Analysis of GCG Y axis')
        xlabel('Frequency')
        ylabel('Power')
        m = length(GCG_Z);       % original sample length
        n = pow2(nextpow2(m));  % transform length
        y = fft(GCG_Z,n);
        f = (0:n-1)*(fs/n);
        power = abs(y).^2/n;
        ax(3)=subplot(325);plot(f(1:floor(n/2)),power(1:floor(n/2)))
        title('Spectral Analysis of GCG Y axis')
        xlabel('Frequency')
        ylabel('Power')
        
        
        
        m = length(SCG_X);       % original sample length
        n = pow2(nextpow2(m));  % transform length
        y = fft(SCG_X,n);
        f = (0:n-1)*(fs/n);
        power = abs(y).^2/n;
        ax(4)=subplot(322);plot(f(1:floor(n/2)),power(1:floor(n/2)))
        title('Spectral Analysis of SCG X axis')
        xlabel('Frequency')
        ylabel('Power')
        
        m = length(SCG_Y);       % original sample length
        n = pow2(nextpow2(m));  % transform length
        y = fft(SCG_Y,n);
        f = (0:n-1)*(fs/n);
        power = abs(y).^2/n;
        ax(5)=subplot(324);plot(f(1:floor(n/2)),power(1:floor(n/2)))
        title('Spectral Analysis of SCG Y axis')
        xlabel('Frequency')
        ylabel('Power')
        
        m = length(SCG_Z);       % original sample length
        n = pow2(nextpow2(m));  % transform length
        y = fft(SCG_Z,n);
        f = (0:n-1)*(fs/n);
        power = abs(y).^2/n;
        ax(6)=subplot(326);plot(f(1:floor(n/2)),power(1:floor(n/2)))
        title('Spectral Analysis of SCG Z axis')
        xlabel('Frequency')
        ylabel('Power')
        linkaxes([ax(6) ax(5) ax(4) ax(3) ax(2) ax(1)],'x');
        
        
        
        gr=1; %%% PLOTING FLAG
        if gr
            figure
            
            ax(1)=subplot(3,2,1);plot(GyroX)
            xlabel('X axis')
            ylabel('angular Velocity')
            ax(2)=subplot(3,2,3);plot(GyroY)
            xlabel('Y axis')
            ylabel('angular Velocity')
            ax(3)=subplot(3,2,5);plot(GyroZ)
            xlabel('Z axis')
            ylabel('angular Velocity')
            ax(4)=subplot(3,2,2);plot(AccX)
            xlabel('X axis')
            ylabel('accelerations')
            ax(5)=subplot(3,2,4);plot(AccY)
            xlabel('Y axis')
            ylabel('accelerations')
            ax(6)=subplot(3,2,6);plot(AccZ)
            xlabel('Z axis')
            ylabel('accelerations')
            linkaxes([ax(6) ax(5) ax(4) ax(3) ax(2) ax(1)],'x');
        end
    end
    %%
    %%%% respiration extraction from ACCELEROMETER X and Y using ARCTAN
    thetaA=atan(AccY./(sqrt(AccX.^2+AccZ.^2)));
    thetaB=atan(AccX./(sqrt(AccY.^2+AccZ.^2)));
    
    hpf=0.1;
    lpf=1;
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
    %%%%% ADR signals from Acc X and Y
    figure
    title('Accelerometeric Derived Respiration Rotation')
    ax(1)=subplot(2,1,1);plot(ThetaA)
    xlabel('X axis')
    ylabel('angular Displacement')
    ax(2)=subplot(2,1,2);plot(ThetaB)
    ylabel('angular Displacement')
    xlabel('Y axis')
    linkaxes([ax(2) ax(1)],'x');

    
    %% Principal component analysis (PCA)
    temp=[RotY';RotX';ThetaA';ThetaB'];
    [coef,resp_pca,vars]=pca(temp');
    
    resp_pca=fft_filter(resp_pca(:,1),fs,hpf,lpf);
    
    figure
    title('MEMS Derived Respirations')
    ax(1)=subplot(6,1,1);plot(RotY)
    ylabel('Angular Disp (rad)')
    xlabel('gyroscope')
    ax(2)=subplot(6,1,2);plot(ThetaA)
    ylabel('Angular Disp (rad) ')
    xlabel('Theta ADR (arctang)')
    ax(3)=subplot(6,1,3);plot(ThetaA_aprx)
    ylabel('Angular Disp (rad)')
    xlabel('Theta ADR (approximation)')
    ax(4)=subplot(6,1,4);plot(AlphaB)
    ylabel('Angular Dispt (rad)')
    xlabel('Theta ADR (arcsin)')
    ax(5)=subplot(6,1,5);plot(AlphaB_aprx)
    ylabel('Acceleration (milli-G)')
    xlabel('Theta ADR (approximation)')
    ax(6)=subplot(6,1,6);plot(resp_pca(:,1))
    ylabel('Angular Dispt (rad)')
    xlabel('Respiration PCA')
    linkaxes([ax(6) ax(5) ax(4) ax(3) ax(2) ax(1)],'x');
    
    figure
    title('Accelerometeric Derived Respiration Rotation')
    ax(1)=subplot(2,1,1);plot(AlphaA)
    xlabel('Alpha A Y axis')
    ylabel('angular Displacement')
    ax(2)=subplot(2,1,2);plot(AlphaB)
    ylabel('angular Displacement')
    xlabel('Alpha B X axis')
    linkaxes([ax(2) ax(1)],'x');
    
    figure
    title('Gyroscopic Derived Respiration Rotation')
    ax(1)=subplot(2,1,1);plot(RotX)
    xlabel('X axis')
    ylabel('angular Displacement')
    ax(2)=subplot(2,1,2);plot(RotY)
    ylabel('angular Displacement')
    xlabel('Y axis')
    linkaxes([ax(2) ax(1)],'x');
    
    %% Median Data Fusion
    fused_four=median([ThetaA ThetaB RotX RotY],2);
    Ax=ThetaA-mean(ThetaA);
    Ay=ThetaB-mean(ThetaB);
    Gx=RotX-mean(RotX);
    Gy=RotY-mean(RotY);
%     R=sqrt((Ax.^2+Ay.^2+Gx.^2+Gy.^2)/4).*sign(mode([Ax Ay Gx Gy],2));
    
    %% signal analysis
    figure
    plot(t,ThetaA,'LineWidth',2,'Linestyle','-')
    hold on
    plot(t,ThetaB,'LineWidth',2,'Linestyle','-')
    hold on
    plot(t,RotX,'LineWidth',2,'Linestyle','-')
    hold on
    plot(t,RotY,'LineWidth',2,'Linestyle','-')
    hold on
    plot(t,resp_pca(:,1),'LineWidth',1.5,'Linestyle','-')
    hold on
    plot(t,fused_four,'LineWidth',1.5,'Linestyle','-')
%     plot(t,R,'LineWidth',1.5,'Linestyle','-')
    ylabel('angular Displacement(rad)')
    
    legend('ADR (acc_Y)','ADR (acc_X)','GDR (gyro_X)','GDR (gyro_Y)','FUSION I (PCA)', 'FUSION II (median)');
 
 

    %% downsampling
    ADR_accX=downsample(ThetaB,32);
    ADR_accY=downsample(ThetaA,32);
    GDR_gyrX=downsample(RotX,32);
    GDR_gyrY=downsample(RotY,32);
    Resp_PCAFusion=downsample(resp_pca(:,1),32);
    Resp_medianFusion=downsample(fused_four,32);
    samples=1:numel(RPM);
    trpm=samples/25;
    %%%%%%%%%%%% Synching MEMS derived Respiration with RPM
    RPM=(RPM-mean(RPM))/std(RPM);
    s7=Resp_PCAFusion;
    s6=Resp_medianFusion;
    s5=-GDR_gyrY;
    s4=GDR_gyrX;
    s3=ADR_accY;
    s2=ADR_accX;
    s1=RPM;
    [RPM_syn,ADR_accX_syn] = alignsignals(s1,s2);
    [x13,ADR_accY_syn] = alignsignals(s1,s3);
    [x14,GDR_gyrX_syn] = alignsignals(s1,s4);
    [x15,GDR_gyrY_syn] = alignsignals(s1,s5);
    [x16,MedianFusion_syn] = alignsignals(s1,s6);
    [x17,PCAFusion_syn] = alignsignals(s1,s7);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT Sychnced signals   
    
    figure
    plot(normalize(RPM_syn,1),'LineWidth',2,'Linestyle','-')
    hold on
    plot(normalize(ADR_accX_syn,1),'LineWidth',2,'Linestyle','-')
    hold on
    plot(normalize(ADR_accY_syn,1),'LineWidth',2,'Linestyle','-')
    hold on
    plot(normalize(GDR_gyrX_syn,1),'LineWidth',2,'Linestyle','-')
    hold on
    plot(normalize(GDR_gyrY_syn,1),'LineWidth',1.5,'Linestyle','-')
    hold on
    plot(normalize(PCAFusion_syn,1),'LineWidth',1.5,'Linestyle','-')
    hold on
    plot(normalize(MedianFusion_syn,1),'LineWidth',1.5,'Linestyle','-')

  ylabel('angular Displacement(rad)')
    
    legend('RPM','ADR (acc_X)','ADR (acc_Y)','GDR (gyro_X)','GDR (gyro_Y)','FUSION I (PCA)', 'FUSION II (Median)');
    
    
    jitter_1=floor(min(length(ADR_accX_syn),length(RPM_syn)));
    jitter_2=floor(min(length(ADR_accY_syn),length(RPM_syn)));
    jitter_3=floor(min(length(GDR_gyrX_syn),length(RPM_syn)));
    jitter_4=floor(min(length(GDR_gyrY_syn),length(RPM_syn)));
    jitter_5=floor(min(length(MedianFusion_syn),length(RPM_syn)));
    jitter_6=floor(min(length(PCAFusion_syn),length(RPM_syn)));

    shimmer_1=find(diff(ADR_accX_syn)>0);
    shimmer_2=find(diff(ADR_accY_syn)>0);
    shimmer_3=find(diff(GDR_gyrX_syn)>0);
    shimmer_4=find(diff(GDR_gyrY_syn)>0);
    shimmer_5=find(diff(MedianFusion_syn)>0);
    shimmer_6=find(diff(PCAFusion_syn)>0);
    
    
  count(k,:)=k;  
    [Sim_Score_ADRx(k,:),Pval_ADRx(k,:)]=corr(normalize(ADR_accX_syn(shimmer_1(1):jitter_1),1),normalize(RPM_syn(shimmer_1(1):jitter_1),1),'Type','Pearson')
       
        [Sim_Score_ADRy(k,:),Pval_ADRy(k,:)]=corr(normalize(ADR_accY_syn(shimmer_2(1):jitter_2),1),normalize(RPM_syn(shimmer_2(1):jitter_2),1),'Type','Pearson')

                [Sim_Score_GDRx(k,:),Pval_GDRx(k,:)]=corr(normalize(GDR_gyrX_syn(shimmer_3(1):jitter_3),1),normalize(RPM_syn(shimmer_3(1):jitter_3),1),'Type','Pearson')

                    [Sim_Score_GDRy(k,:),Pval_GDRy(k,:)]=corr(normalize(GDR_gyrY_syn(shimmer_4(1):jitter_4),1),normalize(RPM_syn(shimmer_4(1):jitter_4),1),'Type','Pearson')

                            [Sim_Score_median(k,:),Pval_median(k,:)]=corr(normalize(MedianFusion_syn(shimmer_5(1):jitter_5),1),normalize(RPM_syn(shimmer_5(1):jitter_5),1),'Type','Pearson')

                                [Sim_Score_pca(k,:),Pval_pca(k,:)]=corr(normalize(PCAFusion_syn(shimmer_6(1):jitter_6),1),normalize(RPM_syn(shimmer_6(1):jitter_6),1),'Type','Pearson')

    
    T_respiration = table(count,Sim_Score_ADRx,Sim_Score_ADRy,Sim_Score_GDRx,Sim_Score_GDRy,Sim_Score_pca,Sim_Score_median,...
    'VariableNames',{'Subject','ADRx','ADRy ','GDRx','GDRy','FusionIPCA','FusionIIMedian'})

 Pvals_respiration = table(count,Pval_ADRx,Pval_ADRy,Pval_GDRx,Pval_GDRy,Pval_pca,Pval_median,...
    'VariableNames',{'Subject','ADRx','ADRy ','GDRx','GDRy','FusionIPCA','FusionIIMedian'})


            
            
              [ LM_Cardiac_5bin, LM_Cardiac_3bin]=peakdetector_cardiacgating(data,fs,1);
%   pause
 
 close all
  
 
correlations=table2array(T_respiration(:,2:end));
 col=@(x)reshape(x,numel(x),1);
            boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),...
                cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)),varargin{:},...
                'Labels',{'ADR_X','ADR_Y','GDR_x','GDR_y','PCA','Median'},'Whisker',1);
            figure
            boxplot2({abs(correlations(:,1)),abs(correlations(:,2)),abs(correlations(:,3)),abs(correlations(:,4)),abs(correlations(:,5)),abs(correlations(:,6)),});
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