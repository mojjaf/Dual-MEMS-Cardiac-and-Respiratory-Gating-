function [ LM_Cardiac_5bin, LM_Cardiac_3bin,MCG_HR,MCG_STI,MCG_CyCPerc,HR_ECG]=peakdetector_cardiacgating(data,fs,gr)
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

%%%%%%%%%%%%%%%%%%%% removing sharp and fast spikes due to mis sampling 
gyrY=medfilt1(double(data.gyroY),7);
gyrX=medfilt1(double(data.gyroX),7);
gyrZ=medfilt1(double(data.gyroZ),7);
accX=medfilt1(double(data.accX),7);
accY=medfilt1(double(data.accY),7);
accZ=medfilt1(double(data.accZ),7);
ECG=medfilt1(double(data.ECG),7);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Downsampling the signals to 100 Hz
data_ds=struct('accX',downsample(accX,8),'accY',downsample(accY,8),'accZ',downsample(accZ,8),'gyroX',downsample(gyrX,8),'gyroY',downsample(gyrY,8),'gyroZ',downsample(gyrZ,8));

SCGX=(data_ds.accX)*(2/32760);
SCGY=(data_ds.accY)*(2/32760);
SCGZ=(data_ds.accZ)*(2/32760);

GCGY=(data_ds.gyroY)*(250/32767);
GCGX=(data_ds.gyroX)*(250/32767);
GCGZ=(data_ds.gyroZ)*(250/32767);
ECG_ds=downsample(ECG,8);
FS_ds=100;
fs=FS_ds;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Spike/Motion artifact removal
SCGX=schmidt_spike_removal(SCGX,FS_ds);
SCGY=schmidt_spike_removal(SCGY,FS_ds);
SCGZ=schmidt_spike_removal(SCGZ,FS_ds);
GCGX=schmidt_spike_removal(GCGX,FS_ds);
GCGY=schmidt_spike_removal(GCGY,FS_ds);
GCGZ=schmidt_spike_removal(GCGZ,FS_ds);

%% %%%%%%%%%%%%%%%% Butterworth filtering to get rid of noise
f1=2; %cuttoff low frequency to get rid of baseline wander
f2=20; %cuttoff frequency to discard high frequency noise
Wn=[f1 f2]*2/fs; % cutt off based on fs
N = 4; % order of 3 less processing
[a,b] = butter(N,Wn); %bandpass filtering
indvec_GCG = validseismodata(SCGX, fs);
% gcg_x = filtfilt(a,b,GCGX);
% GCGx = gcg_x/ median( abs(gcg_x(indvec_GCG)));
% gcg_y = filtfilt(a,b,GCGY);
% GCGy = gcg_y/ median( abs(gcg_y(indvec_GCG)));
% gcg_z = filtfilt(a,b,GCGZ);
% GCGz = gcg_z/ median( abs(gcg_z(indvec_GCG)));
% 
% f1=3; %cuttoff low frequency to get rid of baseline wander
% f2=40; %cuttoff frequency to discard high frequency noise
% Wn=[f1 f2]*2/fs; % cutt off based on fs
% N = 4; % order of 3 less processing
% [a,b] = butter(N,Wn);
indvec_SCG = validseismodata(SCGX, fs);
% 
% scg_hx = filtfilt(a,b,SCGX);
% SCGx = scg_hx/ median( abs(scg_hx(indvec_SCG)));
% scg_hy = filtfilt(a,b,SCGY);
% SCGy = scg_hy/ median( abs(scg_hy(indvec_SCG)));
% scg_hz = filtfilt(a,b,SCGZ);
% SCGz = scg_hz/ median( abs(scg_hz(indvec_SCG)));

ekg = filtfilt(a,b,ECG_ds);
ECG_filt = ekg/ median( abs(ekg(indvec_SCG)));
ECG_z=zscore(ECG_filt);
%%%%%%%%%%%%%%% SSA filter
%  winlen=floor(fs/5);
%         order=3;
%         gr_in=0;
% [GCGx,kk,jj] = ssa(GCGX,winlen,order,gr_in);
% [GCGy,kk,jj]= ssa(GCGY,winlen,order,gr_in);
% [GCGz,kk,jj] = ssa(GCGZ,winlen,order,gr_in);
% [SCGx,kk,jj]= ssa(SCGX,winlen,order,gr_in);
% [SCGy,kk,jj]= ssa(SCGY,winlen,order,gr_in);
% [SCGz,kk,jj]= ssa(SCGZ,winlen,order,gr_in);
hp_A=2;
lp_A=50;
hp_G=2;
lp_G=25;
    SCGx=fft_filter(SCGX,fs,hp_A,lp_A);
    SCGy=fft_filter(SCGY,fs,hp_A,lp_A);
    SCGz=fft_filter(SCGZ,fs,hp_A,lp_A);
    GCGx=fft_filter(GCGX,fs,hp_G,lp_G);
    GCGy=fft_filter(GCGY,fs,hp_G,lp_G);
    GCGz=fft_filter(GCGZ,fs,hp_G,lp_G);
    % %
    % %
    
%% Perform PCA Data Fusion (converting five axis signals into only one chanel signal)

% [coeffs_acc,SCG_pca,variances] = pca([zscore(SCGx)'; zscore(SCGy)'; zscore(SCGz)']');
% [coeffs_gyr,GCG_pca,variances] = pca([zscore(GCGx)';zscore(GCGy)']');
% SCG_pca_norm=zscore(SCG_pca(:,1));
% GCG_pca_norm=zscore(GCG_pca(:,1));
% [coeffs,signals,variances] = pca([SCG_pca_norm';GCG_pca_norm']');

% [coeffs,signals,variances] = pca([SCG_pca_norm';GCG_pca_norm']');

% [coeffs_comb,SCG_GCG_pca,variances] = pca([SCGx'; SCGy'; SCGz'; GCGx';GCGy']');

% Sensor_fusion_AccGyr_sep=signals(:,1);
% Sensor_fusion_AccGyr_all=SCG_GCG_pca(:,1);
%% Perform ICA 
r=2;
% SCG_ica = fastICA([zscore(SCGx)'; zscore(SCGy)'; zscore(SCGz)'],r);
% GCG_ica = fastICA([zscore(GCGx)';zscore(GCGy)';zscore(SCGz)'],r);

Zfica = fastICA([zscore(SCGz)';zscore(GCGy)'],r);
Sensor_fusion_AccGyr_ica=Zfica(1,:);
 %%%%%%%%%%%%%%%%%%%%%%%%%
%  [ env_SCGx ] = Envelope_Filter(SCGX,2 );
%  [ env_SCGy ] = Envelope_Filter(SCGY,2 );
% [ env_SCGz ] = Envelope_Filter(SCGZ,2 );
% [ env_GCGx ] = Envelope_Filter(GCGX,2 );
% [ env_GCGy ] = Envelope_Filter(GCGY,2 );
% figure
% subplot(511);plot(env_SCGx)
% subplot(512);plot(env_SCGy)
% subplot(513);plot(env_SCGz)
% subplot(514);plot(env_GCGx)
% subplot(515);plot(env_GCGy)
%  [coeffs_combenv,SCG_GCG_Envpca,variances] = pca([env_SCGx'; env_SCGy'; env_SCGz'; env_GCGx';env_GCGy']');
%     env_SCG_GCG=baseline_removal_imu(SCG_GCG_Envpca(:,1),0.5,fs);

%%%%%%%%%%%%%%%%%%%%%% Envelope of PCA signal

[ env_pca ] = Envelope_Filter(Sensor_fusion_AccGyr_ica,2 );
   env_pca_filt=baseline_removal_imu(env_pca,0.5,fs);
   env_SCG_GCG=SMQT(env_pca_filt,1,8);

   
%    [ env_SCG ] = Envelope_Filter(SCG_pca_norm,2 );
% [ env_GCG ] = Envelope_Filter(GCG_pca_norm,2 );
% figure
% subplot(311);plot(env_SCG)
% subplot(312);plot(env_GCG)
% subplot(313);plot(env_SCG_GCG)

%% Finding heartbeats on the envelope signal using a rolling window scheme
% Gcomdiv=gcd(numel(env_SCG_GCG),fs*20);
% signal_segmented = reshape(env_SCG_GCG, [], Gcomdiv);
% 
% locs_all_max=[];
% locs_all_min=[];
% 
% [row,col]=size(signal_segmented);
% for i=1:col
%     locs_max = ampd(signal_segmented(:,i));
%     locs_min = ampd(-signal_segmented(:,i));
%     
%     if i<2
%         locs_all_max=[locs_all_max, locs_max];
%         locs_all_min=[locs_all_min, locs_min];
%         
%     else
%         locs_all_max=[locs_all_max, locs_max+numel(signal_segmented(:,2:i))];
%         locs_all_min=[locs_all_min, locs_min+numel(signal_segmented(:,2:i))];
%         
%     end
%     
% end
% %% Filtering false detected peaks 
% 
% PPint=diff(locs_all_max);
% thrs=median(PPint);
% 
% [a,b]=find(PPint< 0.35*thrs | PPint> 1.85*thrs);
% 
% for i=1:(numel(b))
%     start = locs_all_max(b(i))+35; %in samples
%     stop = locs_all_max(b(i)+1)-35;
%     sample = env_SCG_GCG(start:stop);
%     
%     for j=1:numel(sample)
%         if sample(j)==max(env_SCG_GCG(start:stop))
%             Missed_location(i)=start+j-1;
%         end
%     end
%     
% end
% 
% locs_all_max_full=[locs_all_max,Missed_location];
% locs_all_max_tuned=sort(locs_all_max_full);
% 
% PPint2=diff(locs_all_max_tuned);
% thrs=median(PPint2);
% [a,b]=find(PPint2< 0.35*thrs);
% 
% locs_all_max_tuned(b+1)=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Adaptive peak detection and
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% segmentation
FusedICA_sig=schmidt_spike_removal(Sensor_fusion_AccGyr_ica',FS_ds);

FusedPCA_env=schmidt_spike_removal(env_SCG_GCG',FS_ds);


%%%%%%%%%%%%%%%%%MATCH FILTER

[AO_amp,AO_locs]=findpeaks(env_SCG_GCG,'MinPeakHeight',(max(env_SCG_GCG(1:fs*10))*0.5), 'MinPeakDistance',round(fs/2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[AO_amp_raw,AO_i_raw,heartrate_median]=AdaptivePeakDetector((FusedICA_sig),FusedPCA_env,AO_locs,fs,gr);


% [ LM_Cardiac_5bin, LM_Cardiac_3bin ] = cardiac_gate_generator( FusedPCA_sig, (locs_all_max2-40)',fs,1 );

% [ mat_SCG ] = ensemble_average(SCG_pca_norm,AO_i_raw',fs,10,99,0);
% [ mat_GCG ] = ensemble_average(GCG_pca_norm,AO_i_raw',fs,10,99,0);
%  [ mat_SCG_GCG] = ensemble_average(Sensor_fusion_AccGyr_ica,AO_i_raw',fs,10,99,0  );
%  [ mat_envSCG_GCG] = ensemble_average(env_pca_filt,AO_i_raw',fs,10,99,0  );

% [ mat_ECG] = ensemble_average(ECG_z,AO_i_raw',fs,10,99,0  );
% 
% figure
% subplot(131);plot(median(mat_SCG))
% subplot(132);plot(median(mat_GCG))
% subplot(133);plot(median(mat_SCG_GCG))
% 
% y = filter( median(mat_SCG_GCG), 1, Sensor_fusion_AccGyr_ica );

% 
% 
% [ mat_SCGx ] = ensemble_average(SCGx,AO_i_raw',fs,10,99,1);
% [ mat_SCGy ] = ensemble_average(SCGy,AO_i_raw',fs,10,99,1);
% [ mat_SCGz] = ensemble_average(SCGz,AO_i_raw',fs,10,99,1  );
% 
% figure
% subplot(311);plot(median(mat_SCGx))
% subplot(312);plot(median(mat_SCGy))
% subplot(313);plot(median(mat_SCGz))
% 
% 
% [ mat_GCGx ] = ensemble_average(GCGx,AO_i_raw',fs,10,99,1);
% [ mat_GCGy ] = ensemble_average(GCGy,AO_i_raw',fs,10,99,1);
% [ mat_GCGz] = ensemble_average(GCGz,AO_i_raw',fs,10,99,1  );
% 
% figure
% subplot(311);plot(median(mat_GCGx))
% subplot(312);plot(median(mat_GCGy))
% subplot(313);plot(median(mat_GCGz))

%  figure
%  subplot(411);plot(median(mat_SCGz))
%  subplot(412);plot(median(mat_GCGy))
%  subplot(413);plot(median(mat_SCG_GCG))
%  subplot(414);plot(median(mat_envSCG_GCG))
%%
 %%%% Generating cardiac gating files
% [ LM_Cardiac_5bin, LM_Cardiac_3bin ] = cardiac_gate_generator(FusedPCA_sig, AO_i_raw',fs,1 );

[ LM_Cardiac_5bin, LM_Cardiac_3bin ] = cardiac_gate_generator_dynamic(FusedICA_sig, AO_i_raw',fs,1);

if nargout > 2
        disp('Calculating heartRate and systolic time intervals ...')
   
[heartRate,systolicTimeInterval] = getSTIs(FusedICA_sig,20*fs,fs);
MCG_HR=median(heartRate);
MCG_STI=median(systolicTimeInterval);
Cycle_percentages=systolicTimeInterval./(heartRate./60);
MCG_CyCPerc=median(Cycle_percentages);

[HR_ECG, ~] = getHeartRateSchmidt(ECG_ds, fs, 0);
 end
% [heartRateECG,systolicTimeIntervalECG] = getSTIs(ECG_ds,20*fs,fs);



end




