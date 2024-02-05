function [ RespirationCurves ] = respiration_extractor_from_imu( data, fs )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
 %%  PREPROCESS for respiratory gating
%     st=10*fs;
%       sp=numel(data.gyroY)-fs*10;

    st=1;
    sp=numel(data.gyroY);
    gyroY=double(data.gyroY(st:sp))*(250/32767);% 250 corresponds to gyroscope range (dps) and 32767=(2^(16-1)-1) 
    gyroX=double(data.gyroX(st:sp))*(250/32767);
    gyroZ=double(data.gyroZ(st:sp))*(250/32767);
    
    AccX=double(data.accX(st:sp))*(2/32767); % 2 corresponds to the accelerometer dynamic range (in g= gravity)
    AccY=double(data.accY(st:sp))*(2/32767);
    AccZ=double(data.accZ(st:sp))*(2/32767);
    Hp=0.01;
    GyroX=baseline_removal_imu(gyroX,Hp,fs);
    GyroY=baseline_removal_imu(gyroY,Hp,fs);
    GyroZ=baseline_removal_imu(gyroZ,Hp,fs);
    
    %% EXTRACTING RESPIRARION TRACES FROM  ACCELEROMETER X and Y
    %%%% respiration extraction from ACCELEROMETER X and Y using ARCTAN
    thetaA=atan(AccX./(sqrt(AccY.^2+AccZ.^2)));% check how it looks if we filter before angle calculations% should we remove atan??
    thetaB=atan(AccY./(sqrt(AccX.^2+AccZ.^2)));
    
    hpf=0.1;
    lpf=1; % assuming that respiration rate is lower than 1 Hz
    ThetaA=fft_filter(thetaA,fs,hpf,lpf);
    ThetaB=fft_filter(thetaB,fs,hpf,lpf);
    
  %% EXTRACTING RESPIRARION TRACES FROM  GYROSCOPE X, Y , Z
 
    samples=1:numel(AccX);
    t=samples/fs;
    
    rotX=cumtrapz(t,GyroX);
    RotX=fft_filter(rotX*0.017,fs,hpf,lpf); % 0.017 is the degree to radian conversion constant.
    
    rotY=cumtrapz(t,GyroY);
    RotY=fft_filter(rotY*0.017,fs,hpf,lpf);
    
    rotZ=cumtrapz(t,GyroZ);
    RotZ=fft_filter(rotZ*0.017,fs,hpf,lpf);
    
    %% CHECK SIGNAL POLARITY (compared to RPM) THIS is OUTDATED
    %%%%%%%%%
    Gama_zaxis=fft_filter(AccZ,fs,hpf,lpf); % we need to figure out a new MEMS standalone way for polarity check

   [adrY_corr,~]=corr(Gama_zaxis,ThetaB,'Type','Pearson');
   [adrX_corr,~]=corr(Gama_zaxis,ThetaA,'Type','Pearson')
   [gdrX_corr,~]=corr(Gama_zaxis,RotX,'Type','Pearson')
   [gdrY_corr,~]=corr(Gama_zaxis,RotY,'Type','Pearson')
%  [gdrZ_corr,~]=corr(Gama_zaxis,RotZ,'Type','Pearson')
    
    corr_val=[adrX_corr,adrY_corr,gdrX_corr,gdrY_corr];
    sign_val=sign(corr_val);
        resp_signals=[ThetaA,ThetaB,RotX,RotY];
        sign_mat=repmat(sign_val,[size(resp_signals,1),1]);
        new_resp_sig=sign_mat.*resp_signals;
        
        %%%%% Standard Normalization of ADR AND GDR Signals
        ADRx=zscore(new_resp_sig(:,1));
        ADRy=zscore(new_resp_sig(:,2));
        GDRx=zscore(new_resp_sig(:,3));
        GDRy=zscore(new_resp_sig(:,4));
        GDRz=zscore(RotZ);
        
    %% Principal component analysis (PCA)
    temp=[(ADRx');(ADRy');(GDRx');(GDRy');(GDRz')];
    [coef,resp_pca,vars]=pca(temp');   %% perhaps GDRz is not neccessary!
%      resp_pca=fft_filter(resp_pca(:,1),fs,hpf,lpf);
    
     resp_pca=resp_pca(:,1);

    
    %% signal visualization
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
    ylabel('angular Displacement(rad)')
    legend('ADR (acc_Y)','ADR (acc_X)','GDR (gyro_X)','GDR (gyro_Y)', 'GDR (gyro_Z)','PCA');

%%
 RespirationCurves=struct('ADRx',ADRx,'ADRy',ADRy,...
    'GDRx',GDRx,'GDRy',GDRy,'GDRz',GDRz,'PCA',resp_pca);
end

