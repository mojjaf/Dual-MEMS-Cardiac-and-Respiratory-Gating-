% close all
clear all
cd
% datapath = 'D:\BIOSIGNAL_PROCESSING\Physiological Recordings\CoronaryDiseasedPatients/'
datapath = 'D:\BIOSIGNAL_PROCESSING\MEDICAL IMAGING PROJECT\MIC2018 project\Micromachines project/'
% datapath='D:\BIOSIGNAL_PROCESSING\MEDICAL IMAGING PROJECT\PilotStudy_2015\RPM and MEMS data (HealthySubjects)/'
% datapath='D:\BIOSIGNAL_PROCESSING\Atrial Fibrillation\Clinical Data_Confidential\AF measurements may-june 2014/'
filepattern=fullfile(datapath,'*.txt')
txtfiles=dir(filepattern)



for k=2:length(txtfiles)
    
    data=load2_txt([datapath   txtfiles(k).name])
    fs=800;
    
    
    
    %%%%%%%%%%% for respiratory gating
    gyroY=double(data.gyroY)*(250/32767);
    gyroX=double(data.gyroX)*(250/32767);
    gyroZ=double(data.gyroZ)*(250/32767);
    
    
    
    AccX=double(data.accX)*(2/32760);
    AccY=double(data.accY)*(2/32760);
    AccZ=double(data.accZ)*(2/32760);
    Hp=0.01;
    GyroX=baseline_removal_imu(gyroX,Hp,fs);
    GyroY=baseline_removal_imu(gyroY,Hp,fs);
    GyroZ=baseline_removal_imu(gyroZ,Hp,fs);
    % AccX=baseline_removal_imu(accX,Hp,fs);
    % AccY=baseline_removal_imu(accY,Hp,fs);
    % AccZ=baseline_removal_imu(accZ,Hp,fs);
    
    
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
    FFTanalysis=1
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
    
    ThetaA=fft_filter(thetaA,fs,.1,2);
    ThetaB=fft_filter(thetaB,fs,.1,2);
    %   ThetaB=medfilt1(ThetaB,7);
    
    %%%%%%%%%%%%%%%%%%%%%%%  respiration extraction from ACCELEROMETER X and Y
    %%%%%%%%%%%%%%%%%%%%%%%  WITHOUT ARCTAN (APPROXIMATION)
    thetaA_aprx=(AccY./(sqrt(AccX.^2+AccZ.^2)));
    thetaB_aprx=(AccX./(sqrt(AccY.^2+AccZ.^2)));
    ThetaA_aprx=fft_filter(thetaA_aprx,fs,.1,2);
    ThetaB_aprx=fft_filter(thetaB_aprx,fs,.1,2);
    %%%%%%%%%%%%%%%%%%%%%%%%% respiration extraction from ACCELEROMETER X and
    %%%%%%%%%%%%%%%%%%%%%%%%% Y using ARCSIN
    alphaA=asin(AccY);
    alphaB=asin(AccX);
    AlphaA=fft_filter(alphaA,fs,.1,2);
    AlphaB=fft_filter(alphaB,fs,.1,2);
    
    AlphaB_aprx=fft_filter(AccX,fs,.1,2);
    AlphaA_aprx=fft_filter(AccY,fs,.1,2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% respiration extraction from GYROSCOPE X and
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Y by INTEGRATING (convering angular velocity to
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% angular displacement)
    samples=1:numel(AccX);
    t=samples/fs;
    rotX=cumtrapz(t,GyroX);
    RotX=fft_filter(rotX*0.017,fs,.1,2);
    
    rotY=cumtrapz(t,GyroY);
    RotY=fft_filter(rotY*0.017,fs,.1,2);
    
    
    rotZ=cumtrapz(t,GyroZ);
    RotZ=fft_filter(rotZ*0.017,fs,.1,2);
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
    temp=[RotY,RotX,ThetaA,ThetaA_aprx,ThetaB_aprx,ThetaB,AlphaA,AlphaA_aprx,AlphaB,AlphaB_aprx];
    resp_pca=pca(temp');
    resp_pca=fft_filter(resp_pca,fs,.1,2);
    
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
    
    %% Sensor Fusion
    
    fused_four=median([ThetaA ThetaB RotX -RotY],2);
    temp=[-RotY,RotX,ThetaA,ThetaB];
    respfused_pca=pca(temp');
    Kalman_fusion=fusion_4signals_2acc_2gyro (ThetaA',ThetaB',RotX',-RotY');
    
    Ax=ThetaA-mean(ThetaA);
    Ay=ThetaB-mean(ThetaB);
    Gx=RotX-mean(RotX);
    Gy=-RotY-mean(-RotY);
    R=sqrt((Ax.^2+Ay.^2+Gx.^2+Gy.^2)/4).*sign(mode([Ax Ay Gx Gy],2));
    
    %% signal analysis
    
    figure
    plot(t,ThetaA,'LineWidth',2,'Linestyle','-')
    hold on
    plot(t,ThetaB,'LineWidth',2,'Linestyle','-')
    hold on
    plot(t,RotX,'LineWidth',2,'Linestyle','-')
    hold on
    plot(t,-RotY,'LineWidth',2,'Linestyle','-')
    hold on
    plot(t,-respfused_pca(:,1),'LineWidth',1.5,'Linestyle','-')
    hold on
    plot(t,Kalman_fusion-mean(Kalman_fusion),'LineWidth',1.5,'Linestyle','-')
    hold on
    plot(t,fused_four,'LineWidth',1.5,'Linestyle','-')
    plot(t,R,'LineWidth',1.5,'Linestyle','-')
    ylabel('angular Displacement(rad)')
    
    legend('ADR (acc_Y)','ADR (acc_X)','GDR (gyro_X)','GDR (gyro_Y)','FUSION I (PCA)','FUSION II (Kalman)', 'FUSION III (median)', 'Fusion IV ()');
    
    %% loading RPM Signals
%     load( 'plaque2_sub2_23032016' )
    
     [Data_layoutamplitude,phase,timestamp,validflag,ttlin,mark,ttlout] = importRPM('patient_01_06_2016_4DPET.vxp');
%          [Data_layoutamplitude,phase,timestamp,validflag,ttlin,mark,ttlout] = importRPM('0211_2015_1234_1234.vxp');

    RPM=(Data_layoutamplitude-mean(Data_layoutamplitude))/std(Data_layoutamplitude);
    RPM=fft_filter(Data_layoutamplitude,25,0.1,12);
    %  RPM= resample(RPM,32,1);
    %  GDR=downsample(RotX,32);
    %  ADR=downsample(ThetaB,32);
    %% downsampling
    ADR_accX=downsample(ThetaB,32);
    ADR_accY=downsample(ThetaA,32);
    GDR_gyrX=downsample(RotX,32);
    GDR_gyrY=downsample(-RotY,32);
    Resp_KalmanFusion=downsample(-respfused_pca(:,1),32);
    Resp_medianFusion=downsample(Kalman_fusion-mean(Kalman_fusion),32);
    Resp_PCAFusion=downsample(fused_four,32);
    
    %%%%%%%%%%%% Synching MEMS derived Respiration with RPM
    RPM=(RPM-mean(RPM))/std(RPM);
    
    s6=Resp_medianFusion;
    s5=GDR_gyrY;
    s4=GDR_gyrX;
    s3=ADR_accY;
    s2=ADR_accX;
    s1=RPM;
    [RPM_syn,ADR_accX_syn] = alignsignals(s1,s2);
    [x13,ADR_accY_syn] = alignsignals(s1,s3);
    [x14,GDR_gyrX_syn] = alignsignals(s1,s4);
    [x15,GDR_gyrY_syn] = alignsignals(s1,s5);
    [x16,MedianFusion_syn] = alignsignals(s1,s6);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT Sychnced signals
    figure
    plot(normalize(RPM_syn,1))
    hold on
    plot(normalize(ADR_accX_syn,1))
    hold on
    plot(normalize(ADR_accY_syn,1))
    hold on
    plot(normalize(GDR_gyrX_syn,1))
    hold on
    plot(normalize(GDR_gyrY_syn,1))
    hold on
    plot(normalize(MedianFusion_syn,1))
    hold on
    ylabel('Arbitrary Units (norm)')
    
    legend('RPM','ADR (X)','ADR(Y)','GDR(X)','GDR(Y)','Fusion(Acc+Gyr)');
    
    figure
    plot(ADR_accX_syn)
    hold on
    plot(ADR_accY_syn)
    hold on
    plot(GDR_gyrX_syn)
    hold on
    plot(GDR_gyrY_syn)
    hold on
    plot(MedianFusion_syn)
    ylabel('angular Displacement(rad)')
    
    legend('ADR (X)','ADR(Y)','GDR(X)','GDR(Y)','Fusion(Acc+Gyr)');
    
    %% Gating Data Generation
    [ minlocs,maxlocs ] = MinMaxLocFinder( ADR_accX_syn,ADR_accY_syn,GDR_gyrX_syn,GDR_gyrY_syn,Resp_KalmanFusion,MedianFusion_syn, Resp_PCAFusion );
    [ listmode_resp_Phase_5bin ] = PhaseBasedGating(minlocs,MedianFusion_syn);
    [ listmode_resp_AMP_5bin ] = MEMS_Amplitude_GATER( MedianFusion_syn,minlocs,maxlocs);
    read=0;
    if read
        savdir1 ='D:\BIOSIGNAL_PROCESSING\MEDICAL IMAGING PROJECT\MIC2018 project\Gating Files' ;
        RespGating_listmode_AMP = int64(listmode_resp_AMP_5bin);
        RespLM_gyroBindataT=table(RespGating_listmode_AMP);
        RespGating_listmode_PHASE= int64(listmode_resp_Phase_5bin);
        RespLM_PHASE=table(RespGating_listmode_PHASE);
%         writetable(RespLM_gyroBindataT,'LM_RespGatingLM_5Bins_withoutAMP_01062016.txt','Delimiter','\t');
%         writetable(RespLM_PHASE,'LM_RespGatingLM_5Bins_PHASE_01062016.txt','Delimiter','\t');
        writetable(RespLM_gyroBindataT,'LM_RespGatingLM_5Bins_withoutAMP_23032016.txt','Delimiter','\t');
        writetable(RespLM_PHASE,'LM_RespGatingLM_5Bins_PHASE_23032016.txt','Delimiter','\t');
    end
    
    
    %% cardiac gating Generation
    
    [ LM_Cardiac_5bin, LM_Cardiac_3bin]=peakdetector_cardiacgating(data,fs,1);
    read=0
    if read
        
%         savdir1 ='D:\BIOSIGNAL_PROCESSING\MEDICAL IMAGING PROJECT\MIC2018 project\Gating Files' ;
        %  gyro_gatingbins=(listmode_data);
%           save(fullfile(savdir1,sprintf('GyroGateBins_sub_%02d.txt',2)),'listmode_data','-ascii','-tabs');
        
        writetable(LM_Cardiac_5bin,'LMcardiac_MEMS_5Bins_01062016.txt','Delimiter','\t')
        writetable(LM_Cardiac_3bin,'LMcardiac_MEMS_3Bins_01062016.txt','Delimiter','\t')
    end
    
    %  [ LM_CardiacData ] = cardiac_gate_generator( signal, fd_points,gr )
    
    % [peaks_scg,fd_gyro]=peakdetector_dualgating(data,fs,detrend(tmptrigy),0,1);
    
    pause
    close all
end