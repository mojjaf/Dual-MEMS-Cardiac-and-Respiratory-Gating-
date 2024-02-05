function [ LM_Cardiac_5bin, LM_Cardiac_3bin]=peakdetector_cardiacgating_sync(data,fs,shifter, gr)
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

%%%%%%%%%%%%%%%%%%%% removing sharp and fast spikes due to mis-sampling 
gyrY=medfilt1(double(data.gyroY),7);
gyrX=medfilt1(double(data.gyroX),7);
gyrZ=medfilt1(double(data.gyroZ),7);
accX=medfilt1(double(data.accX),7);
accY=medfilt1(double(data.accY),7);
accZ=medfilt1(double(data.accZ),7);
ECG=medfilt1(double(data.ECG),7);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Downsampling the signals to 100 Hz
data_ds=struct('accX',downsample(accX,8),'accY',downsample(accY,8),'accZ',downsample(accZ,8),'gyroX',downsample(gyrX,8),'gyroY',downsample(gyrY,8),'gyroZ',downsample(gyrZ,8));

SCGX=(data_ds.accX)*(2/32760);% should fix it.
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

indvec_SCG = validseismodata(SCGX, fs);
ekg = filtfilt(a,b,ECG_ds);
ECG_filt = ekg/ median( abs(ekg(indvec_SCG)));
ECG_z=zscore(ECG_filt);
%%%%%%%%%%%%%%% BP filter

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
%% Perform ICA 
r=2;
% SCG_ica = fastICA([zscore(SCGx)'; zscore(SCGy)'; zscore(SCGz)'],r);
% GCG_ica = fastICA([zscore(GCGx)';zscore(GCGy)';zscore(SCGz)'],r);

Zfica = fastICA([zscore(SCGz)';zscore(GCGy)'],r);
Sensor_fusion_AccGyr_ica=Zfica(1,:);
 %%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%% Envelope of ICA signal

[ env_ica ] = Envelope_Filter(Sensor_fusion_AccGyr_ica,2 );
   env_ica_filt=baseline_removal_imu(env_ica,0.5,fs);
   env_SCG_GCG=SMQT(env_ica_filt,1,8);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Adaptive peak detection and
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% segmentation
FusedICA_sig=schmidt_spike_removal(Sensor_fusion_AccGyr_ica',FS_ds);

FusedICA_env=schmidt_spike_removal(env_SCG_GCG',FS_ds);


%%%%%%%%%%%%%%%%%MATCH FILTER

[AO_amp,AO_locs]=findpeaks(env_SCG_GCG,'MinPeakHeight',(max(env_SCG_GCG(1:fs*10))*0.5), 'MinPeakDistance',round(fs/2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[AO_amp_raw,AO_i_raw,heartrate_median]=AdaptivePeakDetector((FusedICA_sig),FusedICA_env,AO_locs,fs,gr);


%%
 %%%% Generating cardiac gating files
% [ LM_Cardiac_5bin, LM_Cardiac_3bin ] = cardiac_gate_generator(FusedPCA_sig, AO_i_raw',fs,1 );
input_cardiac_gated_sig=[zeros(1,round(shifter))';FusedICA_sig];
fiducials_shifted=round(shifter)+AO_i_raw';
[ LM_Cardiac_5bin, LM_Cardiac_3bin ] = cardiac_gate_generator_dynamic(input_cardiac_gated_sig,fiducials_shifted ,fs,1);

if nargout > 2
        disp('Calculating heartRate and systolic time intervals ...')
   
[heartRate,systolicTimeInterval] = getSTIs(FusedICA_sig,20*fs,fs);
% MCG_HR=median(heartRate);
% MCG_STI=median(systolicTimeInterval);
% Cycle_percentages=systolicTimeInterval./(heartRate./60);
% MCG_CyCPerc=median(Cycle_percentages);
% 
% [HR_ECG, ~] = getHeartRateSchmidt(ECG_ds, fs, 0);
 end
% [heartRateECG,systolicTimeIntervalECG] = getSTIs(ECG_ds,20*fs,fs);



end




