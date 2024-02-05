function [ LM_Cardiac_5bin, LM_Cardiac_3bin]=peakdetector_cardiacgating_tmp(data,fs,gr)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% Author : Mojtaba Jafari Tadi
% Department of Future Technology, University of Turku
% email : mojjaf@utu.fi
% 

% Any direct or indirect use of this code should be referenced 
% Copyright OCT 2018
%%
if ~isstruct(data)
  error('data must be a structure based input');
end
Fs=fs;

if nargin < 4
    gr = 1;   % on default the function always plots
end


gyrY=medfilt1(double(data.gyroY),7);
gyrX=medfilt1(double(data.gyroX),7);
gyrZ=medfilt1(double(data.gyroZ),7);
accX=medfilt1(double(data.accX),7);
accY=medfilt1(double(data.accY),7);
accZ=medfilt1(double(data.accZ),7);

data_ds=struct('accX',downsample(accX,8),'accY',downsample(accY,8),'accZ',downsample(accZ,8),'gyroX',downsample(gyrX,8),'gyroY',downsample(gyrY,8),'gyroZ',downsample(gyrZ,8));
SCGX=(data_ds.accX)*(2/32760);
SCGY=(data_ds.accY)*(2/32760);
SCGZ=(data_ds.accZ)*(2/32760);

GCGY=(data_ds.gyroY)*(250/32767);
GCGX=(data_ds.gyroX)*(250/32767);
GCGZ=(data_ds.gyroZ)*(250/32767);
FS_ds=100;
fs=FS_ds;


% %%%%%%%%%%%%%%%%%%%FFT Filter
% hp_A=3;
% lp_A=35;
%   
% SCGx=fft_filter(SCGX,100,hp_A,lp_A);
% SCGy=fft_filter(SCGY,100,hp_A,lp_A);
% SCGz=fft_filter(SCGZ,100,hp_A,lp_A);
% % % 
% hp_G=1;
%   lp_G=20;
% GCGx=fft_filter(GCGX,100,hp_G,lp_G);
% GCGy=fft_filter(GCGY,100,hp_G,lp_G);
% GCGz=fft_filter(GCGZ,100,hp_G,lp_G);


%% %%%%%%%%%%%%%%%% Butterworth filter
f1=1; %cuttoff low frequency to get rid of baseline wander
f2=20; %cuttoff frequency to discard high frequency noise
Wn=[f1 f2]*2/fs; % cutt off based on fs
N = 4; % order of 3 less processing
[a,b] = butter(N,Wn); %bandpass filtering
indvec_GCG = validseismodata(SCGX, fs);
gcg_x = filtfilt(a,b,GCGX);
GCGx = gcg_x/ median( abs(gcg_x(indvec_GCG)));
gcg_y = filtfilt(a,b,GCGY);
GCGy = gcg_y/ median( abs(gcg_y(indvec_GCG)));
gcg_z = filtfilt(a,b,GCGZ);
GCGz = gcg_z/ median( abs(gcg_z(indvec_GCG)));

f1=3; %cuttoff low frequency to get rid of baseline wander
f2=35; %cuttoff frequency to discard high frequency noise
Wn=[f1 f2]*2/fs; % cutt off based on fs
N = 4; % order of 3 less processing
[a,b] = butter(N,Wn);
 indvec_SCG = validseismodata(SCGX, fs);

scg_hx = filtfilt(a,b,SCGX);
SCGx = scg_hx/ median( abs(scg_hx(indvec_SCG)));
scg_hy = filtfilt(a,b,SCGY);
SCGy = scg_hy/ median( abs(scg_hy(indvec_SCG)));
scg_hz = filtfilt(a,b,SCGZ);
SCGz = scg_hz/ median( abs(scg_hz(indvec_SCG)));


%%  PCA Data Fusion

[coeffs,SCG_pca,variances] = pca([SCGx'; SCGy'; SCGz']');
[coeffs,GCG_pca,variances] = pca([GCGx';GCGy']');
SCG_pca_norm=zscore(SCG_pca(:,1));
GCG_pca_norm=zscore(GCG_pca(:,1));
[coeffs,signals,variances] = pca([SCG_pca_norm';GCG_pca_norm']');


Sensor_fusion_AccGyr=signals(:,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dt=1/fs;
% %Kalman_acc_fusion = kalman([SCGx'; SCGy'; SCGz'], 1, [0.95;0.05;0.95],diag([1 1 1]), 0.05);
% A=[1 dt dt.^2/2; 0 1 dt ; 0 0 1];
% Kalman_acc_fusion = kalman([SCGx'; SCGy'; SCGz'], A, [0.95 0 0;0 0.05 0; 0 0 0.95],diag([0.1 0.1 0.1].^2),ones(size(A))*0.01.^2);
% 
% % Kalman_acc_fusion = kalman([SCGx'; SCGy'; SCGz'], 1, [0.95;0.05;0.95],diag([0.1 0.1 0.1].^2),  0.01.^2);
% % Kalman_gyr_fusion = kalman([GCGx'; SCGy'], 1, [0.95;0.95], diag([1 1]), 0.05);
% % Kalman_gyr_fusion = kalman([GCGx'; GCGy'], 1, [0.95;0.95], diag([1 1]), 0.05);
% % Kalman_gyr_fusion = kalman([GCGx'; GCGy'], 1, [0.95;0.95], diag([1*pi/180 1*pi/180].^2), 0.1^2);
% Kalman_gyr_fusion = kalman([GCGx'; GCGy'; GCGz'], A, [1 0 0;0 1 0; 0 0 0.1], diag([1*pi/180 1*pi/180 1*pi/180].^2), ones(size(A))*0.1^2);
% 
% % Kalman_fusion_AccGyr = kalman([Kalman_acc_fusion(1:end-1); diff(Kalman_gyr_fusion)], 1, [1.5;1.5], diag([1 1]), 0.005);
% 
% % Kalman_fusion_AccGyr = kalman([Kalman_acc_fusion(1,1:end-1); diff(Kalman_gyr_fusion(1,:))], 1, [0.25;0.75], diag([0.01 1*pi/180].^2), 0.01^2);
% Kalman_fusion_AccGyr = kalman([Kalman_acc_fusion(1,1:end-1); diff(Kalman_gyr_fusion(1,:))], [1 dt;0 1], [1 0;0 1], diag([0.01 1*pi/180].^2), [ 0.01^2 0.01^2; 0.01^2 0.01^2]);
% Kalman_fusion_AccGyr=Kalman_fusion_AccGyr(1,:);

%%%%%%%%%%%%%%%%%%%%%%
%  [ accX_hf ] = HF_sphone(  SCGX,1 );
%  [ accY_hf ] = HF_sphone(  SCGY,1 );
%  [ accZ_hf ] = HF_sphone(  SCGZ,1 );
%  [ gyroX_hf ] = HF_sphone(  GCGX,2 );
%  [ gyroY_hf ] = HF_sphone(  GCGY,2 );
%  [ gyroZ_hf ] = HF_sphone(  GCGZ,2 );
 
 [ Fusion_hf ] = HF_sphone(  Sensor_fusion_AccGyr,2 );
%  Kalman_fusion = kalman([accX_hf'; accY_hf'; accZ_hf';gyroX_hf';gyroY_hf';gyroZ_hf'], 1, [0.5;0.5;0.5;0.5;0.5;0.5], diag([1 1 1 1 1 1]), 0.05);

%  [coeff, env_fuse,variance] = pca([zscore(accX_hf'); zscore(accY_hf'); zscore(accZ_hf');zscore(gyroX_hf');zscore(gyroY_hf');zscore(gyroZ_hf')]');
 
f1=1; %cuttoff low frequency to get rid of baseline wander
f2=4; %cuttoff frequency to discard high frequency noise
Wn=[f1 f2]*2/fs; % cutt off based on fs
N = 3; % order of 3 less processing

[a,b] = butter(N,Wn); %bandpass filtering
%  temp_pca = filtfilt(a,b,env_fuse(:,1));
 temp_pca = filtfilt(a,b,Fusion_hf);

 PCA_fusion_filt = temp_pca/ max( abs(temp_pca));
filteredsignal = (PCA_fusion_filt(:))' + (1:numel(PCA_fusion_filt))/numel(PCA_fusion_filt)/100;
filteredsignal = filteredsignal/max(filteredsignal);

% MPH = max(filteredsignal(FS_ds*10:FS_ds*20))*0.3;
MPH = max(filteredsignal(indvec_SCG))*0.3;

MPD = round(fs/2);
% [pval, locsKalEnv] = findpeaksold(filteredsignal,'MINPEAKHEIGHT',MPH,'MINPEAKDISTANCE',MPD);

 % [pksKalEnv,locsKalEnv]=findpeaks(Kalman_fusion_filt,'MinPeakHeight',(max(Kalman_fusion_filt(FS_ds*10:FS_ds*20))*0.35), 'MinPeakDistance',FS_ds*0.35);
%%
 %%
 Gcomdiv=gcd(numel(filteredsignal),fs*20);
 signal_segmented = reshape(filteredsignal, [], Gcomdiv);
 
 %                    locs_temp=[];
 locs_all_max=[];
 locs_all_min=[];
 
 [row,col]=size(signal_segmented);
 for i=1:col
     locs_max = ampd(signal_segmented(:,i));
     locs_min = ampd(-signal_segmented(:,i));
     
     if i<2
         locs_all_max=[locs_all_max, locs_max];
         locs_all_min=[locs_all_min, locs_min];
         
     else
         locs_all_max=[locs_all_max, locs_max+numel(signal_segmented(:,2:i))];
         locs_all_min=[locs_all_min, locs_min+numel(signal_segmented(:,2:i))];
         
     end
     
 end
%%
                   PPint=diff(locs_all_max);
                  thrs=median(PPint);
                  
                  
                  
                  [a,b]=find(PPint< 0.35*thrs | PPint> 1.85*thrs);
                  
                  for i=1:(numel(b))
                      start = locs_all_max(b(i))+35; %in samples
                      stop = locs_all_max(b(i)+1)-35;
                      sample = filteredsignal(start:stop);
                      
                      for j=1:numel(sample)
                          if sample(j)==max(filteredsignal(start:stop))
                              Missed_location(i)=start+j-1;
                          end
                      end
                      
                  end
                  
                  locs_all_max_full=[locs_all_max,Missed_location];
                  locs_all_max2=sort(locs_all_max_full);
                  
                   PPint2=diff(locs_all_max2);
                  thrs=median(PPint2);
                   [a,b]=find(PPint2< 0.35*thrs);
                  
                   locs_all_max2(b+1)=[];

        [ LM_Cardiac_5bin, LM_Cardiac_3bin ] = cardiac_gate_generator( Sensor_fusion_AccGyr, (locs_all_max2-40)',fs,1 );
        
        
        
        
        if gr
            
            figure
plot(PCA_fusion_filt)
hold on
plot(locs_all_max2,PCA_fusion_filt(locs_all_max2),'r*')


% figure
% plot(PCA_fusion_filt)
% hold on
% plot(locsKalEnv,PCA_fusion_filt(locsKalEnv),'r+')
        end

%%
%  EnvGX=zscore(gyroX_hf);
% EnvGY=zscore(gyroY_hf);
% EnvAZ=zscore(accZ_hf);
% EnvAX=zscore(accX_hf);
% Fused_Env=median([EnvGX,EnvGY,EnvAZ,EnvAX],2);
% [data_ssa,data_hf_ssa]= ssa_Analysis( data_ds,2:3,7,0 );
% 
% EnvGX=zscore(data_hf_ssa.gyroX);
% EnvGY=zscore(data_hf_ssa.gyroY);
% EnvAZ=zscore(data_hf_ssa.accZ);
% EnvAX=zscore(data_hf_ssa.accX);


% [pksGx,locsGX]=findpeaks(EnvGX,'MinPeakHeight',(max(EnvGX(FS_ds*10:FS_ds*20))*0.35), 'MinPeakDistance',FS_ds*0.35);
% [pksGy,locsGY]=findpeaks(EnvGY,'MinPeakHeight',(max(EnvGY(FS_ds*10:FS_ds*20))*0.35), 'MinPeakDistance',FS_ds*0.35);
% [pksAz,locsAZ]=findpeaks(EnvAZ,'MinPeakHeight',(max(EnvAZ(FS_ds*10:FS_ds*20))*0.35), 'MinPeakDistance',FS_ds*0.35);
% [pksAx,locsAX]=findpeaks(EnvAX,'MinPeakHeight',(max(EnvAX(FS_ds*10:FS_ds*20))*0.35), 'MinPeakDistance',FS_ds*0.35);

% [kal_ssa_acc,kk,jj] = ssa(Kalman_acc_fusion,10,2,gr);
% 
% [kal_ssa_acc,kk,jj] = ssa(Kalman_acc_fusion,20,1:3,gr);
% [kal_ssa_gyr,kk,jj] = ssa(Kalman_gyr_fusion,10,2:3,gr);
% comb_kal_sigs=kal_ssa_acc+kal_ssa_gyr;


       
%         [ locs_scg_Fuse ] = adaptivepeakfinder_dualgate( Kalman_fusion_AccGyr,Kalman_fusion_filt,locsKalEnv-3,FS_ds,gr );
%                 [ locs_scgZ ] = adaptivepeakfinder_dualgate( SCGz,KalEnv,locsKalEnv,FS_ds,gr )


%         [ locs_gcg_Fuse ] = adaptivepeakfinder_dualgate( Kalman_gyr_fusion,KalEnv,locsKalEnv,FS_ds,gr )
%          [ locs_gcg_Fuse ] = adaptivepeakfinder_dualgate( kal_ssa_gyr,Kalman_fusion_filt,locsKalEnv-3,FS_ds,gr );

%                   [ locs_comb_sigs ] = adaptivepeakfinder_dualgate( comb_kal_sigs,Kalman_fusion_filt,locsKalEnv-5,FS_ds,gr );

                  %%

                                   
                  

% figure
% plot(SCGz)
% hold on
% % plot(SCGx)
% % plot(locs_scgX,SCGx(locs_scgX),'c*')
% plot(locs_scgZ,SCGz(locs_scgZ),'r*')
% plot(GCGy)
% plot(locs_gcgY,GCGy(locs_gcgY),'g*')
% plot(GCGx)
% plot(locs_gcgX,GCGx(locs_gcgX),'k*')

% 
% 
% 
% 
% all_maxlocs={locs_scg_Fuse',locs_gcgy'};
% % %%
% all_maxlocs_vec=cell2mat(all_maxlocs);
% temp=all_maxlocs_vec;
% max_indx=[];
% j=1;
% i=1;
%  while (numel(temp)>2)
%                   distances= abs(all_maxlocs_vec(i)-temp);
%                   [flag,indx]=find(distances<20);
%                   if numel(flag>2)
%                    max_indx(j)=median(temp(indx));
%                    temp(indx)=[];
%                    j=j+1;
%                   end
%                   i=i+1;
%           
%  end
%  
%  maxlocs=int32(max_indx);
% % 
% % [minGx,locsGX_min]=findpeaks(-EnvGX,'MinPeakHeight',(max(-EnvGX(FS_ds*10:FS_ds*20))*0.35), 'MinPeakDistance',FS_ds*0.5);
% % [minGy,locsGY_min]=findpeaks(-EnvGY,'MinPeakHeight',(max(-EnvGY(FS_ds*10:FS_ds*20))*0.35), 'MinPeakDistance',FS_ds*0.5);
% % [minAz,locsAZ_min]=findpeaks(-EnvAZ,'MinPeakHeight',(max(-EnvAZ(FS_ds*10:FS_ds*20))*0.35), 'MinPeakDistance',FS_ds*0.5);
% % [minAx,locsAX_min]=findpeaks(-EnvAX,'MinPeakHeight',(max(EnvAX(FS_ds*10:FS_ds*20))*0.35), 'MinPeakDistance',FS_ds*0.35);
% % 
% % 
%  figure
%  plot(GCGy)
% hold on
% plot(maxlocs,GCGy(maxlocs),'r*')
% plot(locsGX_min,EnvGX(locsGX_min),'r+')
% 
% plot(EnvGY)
% plot(locsGY,EnvGY(locsGY),'g*')
% plot(locsGY_min,EnvGY(locsGY_min),'g+')
% 
% plot(EnvAZ)
% plot(locsAZ,EnvAZ(locsAZ),'k*')
% plot(locsAZ_min,EnvAZ(locsAZ_min),'k+')
% 
% plot(EnvAX)
% plot(locsAX,EnvAX(locsAX),'c*')
% plot(locsAX_min,EnvAX(locsAX_min),'c+')
% 
% all_minlocs={locsGX_min',locsGY_min',locsAZ_min',locsAX_min'};
% 
% all_minlocs_vec=cell2mat(all_minlocs);
% temp_min=all_minlocs_vec;
% min_indx=[];
% jj=1;
% ii=1;
%  while (numel(temp_min)>2)
%                   distances= abs(all_minlocs_vec(ii)-temp_min);
%                   [flag,indx]=find(distances<100);
%                   if numel(flag>2)
%                    min_indx(jj)=median(temp_min(indx));
%                    temp_min(indx)=[];
%                    jj=jj+1;
%                   end
%                   ii=ii+1;
%           
%  end
%  
%         
%        minlocs=int32(min_indx);
% 
% 
% figure
% plot(EnvGY)
% hold on
% plot(maxlocs,EnvGY(maxlocs),'go')
% plot(EnvGY)
% hold on
% plot(minlocs,EnvGY(minlocs),'co')
% 
% 
% 
% 
% if motionArtifact
% 
% indvec=validseismodata(accZ, Fs);
% 
% 
% 
% readingJitter = Fs*2;
% start=indvec(readingJitter);
% stop=indvec(end-readingJitter);
% samples=1:numel(gyrY(start:stop));
% t=samples/Fs;
% 
% GyrY=fft_filter(gyrY(start:stop),fs,1,20);
% GyrX=fft_filter(gyrX(start:stop),fs,1,20);
% GyrZ=fft_filter(gyrZ(start:stop),fs,1,20);
% 
% AccZ=fft_filter(accZ(start:stop),fs,4,40);
% AccY=fft_filter(accY(start:stop),fs,4,40);
% AccX=fft_filter(accX(start:stop),fs,4,40);
% 
% % ECG=ECG(start:stop);
% % sig_acc_z=fft_filter(AccZ,800,4,40);
% % sig_gyr_y=fft_filter(GyrY,800,1,20);
% SCG=accZ(start:stop);
% GCG=gyrX(start:stop);
% else
% samples=1:numel(gyrY);
% t=samples/Fs;
% % 
% GyrY=fft_filter(gyrY,fs,1,20);
% GyrX=fft_filter(gyrX,fs,1,20);
% GyrZ=fft_filter(gyrZ,fs,1,20);
% 
% AccZ=fft_filter(accZ,fs,4,40);
% AccY=fft_filter(accY,fs,4,40);
% AccX=fft_filter(accX,fs,4,40); 
% 
% 
% 
% 
%  SCG=AccZ;
%  GCG=GyrY;
% % data_maf=struct('accX',accX,'accY',accY,'accZ',accZ,'gyroX',gyrX,'gyroY',gyrY,'gyroZ',gyrZ);
% %   data_ssa=ssa_Analysis( data_maf,1,12 );
% % SCG_Z=data_ssa.accZ;
% % GCG_Y=data_ssa.gyroY;
% end 
% 
% 
% % CardiacFusion=median([GyrY GyrX AccZ AccX],2);
% % temp_car=[GyrY GyrX AccZ AccX];
% % cardFus_pca=pca(temp_car');
% %  [tmptrigy1 ] = scg_filt(cardFus_pca(:,1),1 );
% %  [tmptrigy2 ] = scg_filt(cardFus_pca(:,2),2 );
% % tmptrigy= tmptrigy1+tmptrigy2;
% % figure
% % ax(1)=subplot(311);plot(detrend(tmptrigy1))
% % title('PC 1 Envelope')
% % 
% % ax(2)=subplot(312);plot(detrend(tmptrigy2))
% % title('PC 2Envelope')
% % ax(3)=subplot(313);plot(detrend(tmptrigy))
% % title('PC 1 + PC 2 Envelope')
% % linkaxes([ ax(3) ax(2) ax(1)],'x');
% 
% 
% % ECG=filter([1/15, 1/15, 1/15],1,ECG);

%% 
% sig_comb_m=waveletfilt;
%% 





end




