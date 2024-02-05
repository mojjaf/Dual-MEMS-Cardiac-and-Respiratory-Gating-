function [AO_amp_raw,AO_i_raw,heartrate_median]=AdaptivePeakDetector(inputsig1,inputsig2,locs,fs,gr)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
% Complete implementation of Seismocardiogram and Gyrocardiogram peak detection algorithm

%% Inputs
% data : raw measured data including ECG, 3 axis Accelerometer(AccX, AccY and AccZ) and 3 Axes
% Gyroscope signals(GyrX,GyrY and GyrZ)
% fs : sampling frequency e.g. 800 Hz
% gr : flag to plot or not plot (set it 1 to have a plot or set it zero not
% to see any plots
% flag: 2 to have peak detection in Gyroscope signal(Y axis) and 1 to have peak
% detection in Accelerometer (Z axis)
%% Outputs
% AO_amp_raw : amplitude of AO waves amplitudes
% qrs_i_raw : index/ location of AO waves
% delay : number of samples which the signal is delayed due to the
% filtering
%% Method :

%% PreProcessing
% % 1) Signal is preprocessed, first the raw motion signals obtained from each axes (
% i.e.  X,Y and Z)of motion sensors were combined together by summing up the acceleration and rotational informa-
% tion  to  form  a  new  signal.   The  new  summed  signals  called  3D-SCG  and
% 3D-GCG  are  pre-processed  by  applying  a  moving  average  ?lter  as  well  as
% band-pass ?ltering to reduce noise,  baseline wander,  muscle noise and etc 
% (a combination of low pass and high pass filter 5-15 Hz)
% % to get rid of the baseline wander and muscle noise. 


%% Decision Rule 
% At this point in the algorithm, the preceding stages have produced a roughly pulse-shaped
% waveform at the output of the MWI . The determination as to whether this pulse
% corresponds to a QRS complex (as opposed to a high-sloped T-wave or a noise artefact) is
% performed with an adaptive thresholding operation and other decision
% rules outlined below;

% a) FIDUCIAL MARK - The waveform is first processed to produce a set of weighted unit
% samples at the location of the MWI maxima. This is done in order to localize the QRS
% complex to a single instant of time. The w[k] weighting is the maxima value.

% b) THRESHOLDING - When analyzing the amplitude of the MWI output, the algorithm uses
% two threshold values (THR_SIG and THR_NOISE, appropriately initialized during a brief
% 2 second training phase) that continuously adapt to changing ECG signal quality. The
% first pass through y[n] uses these thresholds to classify the each non-zero sample
% (CURRENTPEAK) as either signal or noise:
% If CURRENTPEAK > THR_SIG, that location is identified as a “QRS complex
% candidate” and the signal level (SIG_LEV) is updated:
% SIG _ LEV = 0.125 ×CURRENTPEAK + 0.875× SIG _ LEV

% If THR_NOISE < CURRENTPEAK < THR_SIG, then that location is identified as a
% “noise peak” and the noise level (NOISE_LEV) is updated:
% NOISE _ LEV = 0.125×CURRENTPEAK + 0.875× NOISE _ LEV
% Based on new estimates of the signal and noise levels (SIG_LEV and NOISE_LEV,
% respectively) at that point in the ECG, the thresholds are adjusted as follows:
% THR _ SIG = NOISE _ LEV + 0.25 × (SIG _ LEV ? NOISE _ LEV )
% THR _ NOISE = 0.5× (THR _ SIG)
% These adjustments lower the threshold gradually in signal segments that are deemed to
% be of poorer quality.


% c) SEARCHBACK FOR MISSED QRS COMPLEXES - In the thresholding step above, if
% CURRENTPEAK < THR_SIG, the peak is deemed not to have resulted from a QRS
% complex. If however, an unreasonably long period has expired without an abovethreshold
% peak, the algorithm will assume a QRS has been missed and perform a
% searchback. This limits the number of false negatives. The minimum time used to trigger
% a searchback is 1.66 times the current R peak to R peak time period (called the RR
% interval). This value has a physiological origin - the time value between adjacent
% heartbeats cannot change more quickly than this. The missed QRS complex is assumed
% to occur at the location of the highest peak in the interval that lies between THR_SIG and
% THR_NOISE. In this algorithm, two average RR intervals are stored,the first RR interval is 
% calculated as an average of the last eight QRS locations in order to adapt to changing heart 
% rate and the second RR interval mean is the mean 
% of the most regular RR intervals . The threshold is lowered if the heart rate is not regular 
% to improve detection.

% d) ELIMINATION OF MULTIPLE DETECTIONS WITHIN REFRACTORY PERIOD - It is
% impossible for a legitimate QRS complex to occur if it lies within 200ms after a previously
% detected one. This constraint is a physiological one – due to the refractory period during
% which ventricular depolarization cannot occur despite a stimulus[1]. As QRS complex
% candidates are generated, the algorithm eliminates such physically impossible events,
% thereby reducing false positives.

% e) T WAVE DISCRIMINATION - Finally, if a QRS candidate occurs after the 200ms
% refractory period but within 360ms of the previous QRS, the algorithm determines
% whether this is a genuine QRS complex of the next heartbeat or an abnormally prominent
% T wave. This decision is based on the mean slope of the waveform at that position. A slope of
% less than one half that of the previous QRS complex is consistent with the slower
% changing behaviour of a T wave – otherwise, it becomes a QRS detection.
% Extra concept : beside the points mentioned in the paper, this code also
% checks if the occured peak which is less than 360 msec latency has also a
% latency less than 0,5*mean_RR if yes this is counted as noise

% f) In the final stage , the output of R waves detected in smoothed signal is analyzed and double
% checked with the help of the output of the bandpass signal to improve
% detection and find the original index of the real R waves on the raw ecg
% signal



%% Author : Mojtaba Jafari Tadi
% Technology Research Center, University of Turku
% email : mojjaf@utu.fi
% 

% Any direct or indirect use of this code should be referenced 
% Copyright JUNE 2015
%%
% if ~isstruct(data)
%   error('data must be a structure based input');
% end
Fs=fs;

if nargin < 4
    gr = 1;   % on default the function always plots
end
% if nargin<5
%     coef=[1 1];
% end
% if isfield(data,'EKG')
%   ECG=double(data.EKG);
%     
%  elseif  isfield(data,'ECG2')
% 
%  ECG=double(data.ECG2);
% 
% elseif isfield(data,'ECG')
% ECG=double(data.ECG);
% else
%     ECG=double(data.accX);
% 
% end
% ECG=medfilt1(ECG,3);
% 
% GyrY=double(data.gyroY)*(2000/32767)/256;
% GyrX=double(data.gyroX)*(2000/32767)/256;
% GyrZ=double(data.gyroZ)*(2000/32767)/256;
% AccX=double(data.accX)*(2000/32760);
% AccY=double(data.accY)*(2000/32760);
% AccZ=double(data.accZ)*(2000/32760);
% ECG=fft_filter(ECG,800,5,40);
% ECG = detrend((ECG-mean(ECG))/std(ECG));
% tt=1:length(GyrY);
% R=[AccX';AccY';AccZ'];
% R2=[GyrX';GyrY';GyrZ'];
% comb_signal_acc = max_ratio_combination3(tt,R,3,2400);
% comb_signal_gyr = max_ratio_combination3(tt,R2,2,2400);
% 
% comb_signal_gyr=fft_filter(comb_signal_gyr,800,1,20);
% comb_signal_acc=fft_filter(comb_signal_acc,800,4,40);
% indvec=validseismodata(GyrY, Fs);
% 
%         
%     
%     
% start=indvec(14000);
% stop=indvec(end-4000);
% samples=1:numel(GyrY(start:stop));
% t=samples/Fs;
% 
% 
% 
% 
% GyrY=fft_filter(GyrY(start:stop),800,1,20);
% GyrX=fft_filter(GyrX(start:stop),800,1,20);
% GyrZ=fft_filter(GyrZ(start:stop),800,1,20);
% 
% AccZ=fft_filter(AccZ(start:stop),800,4,40);
% AccY=fft_filter(AccY(start:stop),800,4,40);
% AccX=fft_filter(AccX(start:stop),800,4,40);
% 
% ECG=ECG(start:stop);
% 
% tt=1:length(GyrY);
% R=[AccX';AccY';AccZ'];
% R2=[GyrX';GyrY';GyrZ'];
% comb_signal_acc = comb_signal_acc(start:stop);
% comb_signal_gyr = comb_signal_gyr(start:stop);
% 
% 
% % sig_acc_z=fft_filter(AccZ,800,4,40);
% % sig_gyr_y=fft_filter(GyrY,800,1,20);
% sig_acc_z=AccZ;
% sig_gyr_y=GyrY;
% ECG=filter([1/15, 1/15, 1/15],1,ECG);
%% Initialize
AO_c =[]; %amplitude of R
AO_i =[]; %index
SIG_LEV = 0; 
nois_c =[];
nois_i =[];
delay = 0;
skip = 0; % becomes one when a T wave is detected
not_nois = 0; % it is not noise when not_nois = 1
selected_RR =[]; % Selected RR intervals
m_selected_AOAO = 0;
mean_AO = 0;
AO_i_raw =[];
AO_amp_raw=[];
ser_back = 0; 
test_m = 0;
SIGL_buf = [];
NOISL_buf = [];
THRS_buf = [];
SIGL_buf1 = [];
NOISL_buf1 = [];
THRS_buf1 = [];

%% Plot ECG 
% if gr
% figure
% ax(1)=subplot(4,2,[1 2]);plot(ECG)
% title('ECG signal')
% xlabel('ECG')
% ylabel('Amplitude')
% ax(2)=subplot(4,2,3);plot(GyrX)
% xlabel('X axis')
% ylabel('angular Velocity')
% ax(3)=subplot(4,2,5);plot(GyrY)
% xlabel('Y axis')
% ylabel('angular Velocity')
% ax(4)=subplot(4,2,7);plot(GyrZ)
% xlabel('Z axis')
% ylabel('angular Velocity')
% ax(5)=subplot(4,2,4);plot(AccX)
% xlabel('X axis')
% ylabel('accelerations')
% ax(6)=subplot(4,2,6);plot(AccY)
% xlabel('Y axis')
% ylabel('accelerations')
% ax(7)=subplot(4,2,8);plot(AccZ)
% xlabel('Z axis')
% ylabel('accelerations')
% linkaxes([ax(7) ax(6) ax(5) ax(4) ax(3) ax(2) ax(1)],'x');
% end
% if gr
%  if fs == 800
%   figure,  ax(1)=subplot(321);plot(ECG);axis tight;title('Raw ECG Signal');
%  else
%   figure,  ax(1)=subplot(3,2,[1 2]);plot(ECG);axis tight;title('Raw ECG Signal');
%  end
% end    
% %% Summing up motion signals
% if fs == 800
% % sig_acc=detrend((AccZ+AccX+AccY));
% sig_acc=comb_signal_acc; %combined signal MRC
% % sig_acc=detrend((AccZ+AccY));
% 
% % sig_acc_z=sig_acc;
% if gr
% ax(2)=subplot(322);plot(sig_acc);axis tight;title('Raw combined 3D SCG');
% end
% 
% % sig_gyr=(GyrY+GyrX+GyrZ);
% sig_gyr=comb_signal_gyr; %combined signal MRC
% % sig_gyr=(GyrY+GyrX);
% % sig_gyr_y=sig_gyr;
% if gr
% ax(3)=subplot(323);plot(sig_gyr);axis tight;title('Raw combined 3D GCG');
% end
% else
% 
% end
% %% Hilbert Transform and Noise cancelation(Filtering) % Moving average and Bandpass Filters 
% % (MV= 5 samples, cutt of freq for BP Filter 4-40 Hz for SCG and 1-15 Hz for GCG)
% 
% sig_acc_filt= Transform( sig_acc,Fs,1 );
% 
% 
% if gr
% ax(4)=subplot(324);plot(sig_acc_filt);axis tight;title('Transformed 3D SCG');
% end
% 
% sig_gyr_filt= Transform( sig_gyr,Fs ,2);
% 
% if gr
% ax(5)=subplot(325);plot(sig_gyr_filt);axis tight;title('Transformed 3D GCG');
% end
% 
% 
% 
% %% Combining transfomed version of SCG and GCG by summing up 3D-SCG and 3D-GCG
% acc_int=cumtrapz(tt,sig_acc_filt);
% gyr_int=cumtrapz(tt,sig_gyr_filt);
% % sig_comb=sig_acc_filt+sig_gyr_filt*10;
% % [ acc_int ] = normalize( acc_int );
% % [ gyr_int ] = normalize( sig_gyr_filt );
% % coef=[0 1];
% j=coef(1);
% k=coef(2);
% 
% sig_comb_int=acc_int*j+gyr_int*k;
% 
% % sig_comb_int=cumtrapz(t,sig_comb);
% sig_comb_int=filter(0.1,1,sig_comb_int);
% % sig_comb_int=filter([1/25, 1/25, 1/25],1,sig_comb_int);
% sig_comb_int=filter([1/3, 1/3, 1/3],1,sig_comb_int);
% 
% sig_comb_int = conv(sig_comb_int ,ones(1 ,round(0.150*fs))/round(0.150*fs));
% % delay = delay + 15;
% envelope_of_fusedsig=sig_comb_int;
% if gr
% ax(6)=subplot(326);plot(sig_comb_int);axis tight;title('Averaged with 15 samples length,Black noise,Green Adaptive Threshold,RED Sig Level,Red circles QRS adaptive threshold');
% linkaxes([ax(6) ax(5) ax(4) ax(3) ax(2) ax(1)],'x');
% end

%% Navigator signal for Peak Detection, i.e. Integrated version of Combined 3D SCG and GCG 
% Note : a minimum distance of 350 samples is considered between each heart beat wave
% since in physiological point of view no RR wave can occur in less than
% 200-300 msec distance
% wander=abs(max(sig_comb_int)-min(sig_comb_int))/2;% 
% envelope_of_fusedsig=sig_comb_int+wander;
% % [~,locs_comb]=findpeaks(sig_comb,'MinPeakHeight',(max(sig_comb)*0.3), 'MinPeakDistance',300);
% [pks,locs]=findpeaks(envelope_of_fusedsig,'MinPeakHeight',(max(envelope_of_fusedsig(1:fs*10))*0.25), 'MinPeakDistance',350);



fusedsignal=inputsig1;
envelope_of_fusedsig=inputsig2;

    

%% initialize the training phase (2 seconds of the signal) to determine the THR_SIG and THR_NOISE
THR_SIG = max(envelope_of_fusedsig(1:2*fs))*1/2; % 0.25 of the max amplitude 
THR_NOISE = mean(envelope_of_fusedsig(1:2*fs))*1/2; % 0.5 of the mean signal is considered to be noise**2
SIG_LEV= THR_SIG;
NOISE_LEV = THR_NOISE;


%% Initialize bandpath filter threshold(2 seconds of the bandpass signal)
THR_SIG1 = max(fusedsignal(5:10*fs))*1/2; % 0.25 of the max amplitude 
THR_NOISE1 = mean(fusedsignal(5:10*fs))*1/2; %
SIG_LEV1 = THR_SIG1; % Signal level in Bandpassed filter
NOISE_LEV1 = THR_NOISE1; % Noise level in Bandpassed filter
%% Thresholding and online desicion rule
pks=envelope_of_fusedsig(locs);
for i = 1 : length(pks)
    
   %% locate the corresponding peak in the filtered signal 
    if locs(i)-round(0.70*fs)>= 1 && locs(i)<= length(fusedsignal)
          [y_i x_i] = max(fusedsignal(locs(i)-round(0.7*fs):locs(i)));
       else
          if i == 1
            [y_i x_i] = max(fusedsignal(1:locs(i)));
            ser_back = 1;
          elseif locs(i)>= length(fusedsignal)
            [y_i x_i] = max(fusedsignal(locs(i)-round(0.70*fs):end));
          end
        
     end
    
    
  %% update the heart_rate (Two heart rate means one the moste recent and the other selected)
    if length(AO_c) >= 9 
        
        diffRR = diff(AO_i(end-8:end)); %calculate AO-AO interval
        mean_AO = mean(diffRR); % calculate the mean of 8 previous AO waves interval
        comp =AO_i(end)-AO_i(end-1); %latest AO-AO
        if comp <= 0.92*mean_AO || comp >= 1.16*mean_AO
            % lower down thresholds to detect better in MVI
                THR_SIG = 0.5*(THR_SIG);
                %THR_NOISE = 0.5*(THR_SIG);  
               % lower down thresholds to detect better in Bandpass filtered 
                THR_SIG1 = 0.5*(THR_SIG1);
                %THR_NOISE1 = 0.5*(THR_SIG1); 
                
        else
            m_selected_AOAO = mean_AO; %the latest regular beats mean
        end 
          
    end
    
      %% calculate the mean of the last 8 AO waves to make sure that AO is not
       % missing(If no AO detected , trigger a search back) 1.66*mean
       
       if m_selected_AOAO
           test_m = m_selected_AOAO; %if the regular AO-AO availabe use it   
       elseif mean_AO && m_selected_AOAO == 0
           test_m = mean_AO;   
       else
           test_m = 0;
       end
        
    if test_m
          if (locs(i) - AO_i(end)) >= round(1.66*test_m)% it shows a AO is missed 
              [pks_temp,locs_temp] = max(envelope_of_fusedsig(AO_i(end)+ round(0.250*fs):locs(i)-round(0.250*fs))); % search back and locate the max in this interval
              locs_temp = AO_i(end)+ round(0.250*fs) + locs_temp -1; %location 
             
              if pks_temp > THR_NOISE
               AO_c = [AO_c pks_temp];
               AO_i = [AO_i locs_temp];
              
               % find the location in filtered sig
               if locs_temp <= length(fusedsignal)
                [y_i_t x_i_t] = max(fusedsignal(locs_temp-round(0.70*fs):locs_temp));
               else
                [y_i_t x_i_t] = max(fusedsignal(locs_temp-round(0.70*fs):end));
               end
               % take care of bandpass signal threshold
               if y_i_t > THR_NOISE1 
                        
                      AO_i_raw = [AO_i_raw locs_temp-round(0.30*fs)+ (x_i_t - 1)];% save index of bandpass 
                      AO_amp_raw =[AO_amp_raw y_i_t]; %save amplitude of bandpass 
                      SIG_LEV1 = 0.25*y_i_t + 0.75*SIG_LEV1; %when found with the second thres 
               end
               
               not_nois = 1;
               SIG_LEV = 0.25*pks_temp + 0.75*SIG_LEV ;  %when found with the second threshold             
             end 
              
          else
              not_nois = 0;
              
          end
    end
      
    
    
    
    %%  find noise and QRS peaks
    if pks(i) >= THR_SIG
        
                 % if a AO candidate occurs within 360ms of the previous QRS
                 % ,the algorithm determines if its AC wave or AO
                 if length(AO_c) >= 3
                      if (locs(i)-AO_i(end)) <= round(0.3600*fs)
                        Slope1 = mean(diff(envelope_of_fusedsig(locs(i)-round(0.075*fs):locs(i)))); %mean slope of the waveform at that position
                        Slope2 = mean(diff(envelope_of_fusedsig(AO_i(end)-round(0.075*fs):AO_i(end)))); %mean slope of previous AO wave
                             if abs(Slope1) <= abs(0.5*(Slope2)) || (locs(i)-AO_i(end)) <= round(0.4*test_m)  % slope less then 0.5 of previous AO
                                 nois_c = [nois_c pks(i)];
                                 nois_i = [nois_i locs(i)];
                                 skip = 1; % T wave identification
                                 % adjust noise level in both filtered and
                                 % MVI
                                 NOISE_LEV1 = 0.125*y_i + 0.875*NOISE_LEV1;
                                 NOISE_LEV = 0.125*pks(i) + 0.875*NOISE_LEV; 
                             else
                                 skip = 0;
                             end
            
                      end
                 end
        
        if skip == 0  % skip is 1 when a AC wave is detected       
        AO_c = [AO_c pks(i)];
        AO_i = [AO_i locs(i)];
        
        % bandpass filter check threshold
         if y_i >= THR_SIG1
                        if ser_back 
                           AO_i_raw = [AO_i_raw x_i];  % save index of bandpass 
                        else
                           AO_i_raw = [AO_i_raw locs(i)-round(0.70*fs)+ (x_i - 1)];% save index of bandpass 
                        end
                           AO_amp_raw =[AO_amp_raw y_i];% save amplitude of bandpass 
          SIG_LEV1 = 0.125*y_i + 0.875*SIG_LEV1;% adjust threshold for bandpass filtered sig
         end
         
        % adjust Signal level
        SIG_LEV = 0.125*pks(i) + 0.875*SIG_LEV ;
        end
        
        
    elseif THR_NOISE <= pks(i) && pks(i)<THR_SIG
        
         %adjust Noise level in filtered sig
         NOISE_LEV1 = 0.125*y_i + 0.875*NOISE_LEV1;
         %adjust Noise level in MVI
         NOISE_LEV = 0.125*pks(i) + 0.875*NOISE_LEV; 
        
        
      
    elseif pks(i) < THR_NOISE
        nois_c = [nois_c pks(i)];
        nois_i = [nois_i locs(i)];
        
        % noise level in filtered signal
        NOISE_LEV1 = 0.125*y_i + 0.875*NOISE_LEV1;
        %end
        
         %adjust Noise level in MVI
        NOISE_LEV = 0.125*pks(i) + 0.875*NOISE_LEV;  
        
           
    end
    
    
    
 
    
    %% adjust the threshold with SNR
    if NOISE_LEV ~= 0 || SIG_LEV ~= 0
        THR_SIG = NOISE_LEV + 0.25*(abs(SIG_LEV - NOISE_LEV));
        THR_NOISE = 0.5*(THR_SIG);
    end
    
    % adjust the threshold with SNR for bandpassed signal
    if NOISE_LEV1 ~= 0 || SIG_LEV1 ~= 0
        THR_SIG1 = NOISE_LEV1 + 0.25*(abs(SIG_LEV1 - NOISE_LEV1));
        THR_NOISE1 = 0.5*(THR_SIG1);
    end
    
    
% take a track of thresholds of smoothed signal
SIGL_buf = [SIGL_buf SIG_LEV];
NOISL_buf = [NOISL_buf NOISE_LEV];
THRS_buf = [THRS_buf THR_SIG];

% take a track of thresholds of filtered signal
SIGL_buf1 = [SIGL_buf1 SIG_LEV1];
NOISL_buf1 = [NOISL_buf1 NOISE_LEV1];
THRS_buf1 = [THRS_buf1 THR_SIG1];



    
 skip = 0; %reset parameters
 not_nois = 0; %reset parameters
 ser_back = 0;  %reset bandpass param   
end

if gr
hold on,scatter(AO_i,AO_c,'m');
hold on,plot(locs,NOISL_buf,'--k','LineWidth',2);
hold on,plot(locs,SIGL_buf,'--r','LineWidth',2);
hold on,plot(locs,THRS_buf,'--g','LineWidth',2);
zoom on;

end

%%%%%%%%%%%%%%%%%%%
% [pk,ekg_locs] = findpeaks(ECG,'MinPeakHeight',(max(ECG)*0.3), 'MinPeakDistance',Fs*0.3);
% % [pk,ekg_locs,delay]=pan_tompkin(ECG,Fs,1)
% % ekg_locs=(ekg_locs);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% kick out the outlier picks
% [a,b]=find(AO_amp_raw<(mean(AO_amp_raw)*0.35));
% [c,d]=find(pk<(mean(pk)*0.35));
% 
% AO_amp_raw(b)=[];
% AO_i_raw(b)=[];
% pk(c)=[];
% ekg_locs(c)=[];
% [aa,bb]=find(diff(AO_i_raw)<=300);
% AO_amp_raw(bb)=[];
% AO_i_raw(bb)=[];
%% overlay on the signals
if gr
figure,az(1)=subplot(311);plot(fusedsignal);

title('Cardiac Pulsataion Detection in Fusion Signal');
axis tight;
hold on,scatter(AO_i_raw,AO_amp_raw,'r');
hold on,plot(locs,NOISL_buf1,'LineWidth',2,'Linestyle','--','color','k');
hold on,plot(locs,SIGL_buf1,'LineWidth',2,'Linestyle','-.','color','g');
hold on,plot(locs,THRS_buf1,'LineWidth',2,'Linestyle','-.','color','c');
az(2)=subplot(312);plot(envelope_of_fusedsig);title('AO on MVI signal and Noise level(black),Signal Level (g) and Adaptive Threshold(cyan)');axis tight;
hold on,scatter(AO_i,AO_c,'r');
hold on,plot(locs,NOISL_buf,'LineWidth',2,'Linestyle','--','color','k');
hold on,plot(locs,SIGL_buf,'LineWidth',2,'Linestyle','-.','color','g');
hold on,plot(locs,THRS_buf,'LineWidth',2,'Linestyle','-.','color','c');
az(3)=subplot(313);plot(fusedsignal-mean(fusedsignal));title('Pulse train of the found AO on MEMS signal');axis tight;
line(repmat(AO_i_raw,[2 1]),repmat([min(fusedsignal-mean(fusedsignal))/2; max(fusedsignal-mean(fusedsignal))/2],size(AO_i_raw)),'LineWidth',2.5,'LineStyle','-.','Color','r');
linkaxes(az,'x');
zoom on;
end
% [CTI_values ] = STIPD( fusedsignal,AO_i_raw,ekg_locs,flag,gr);
[ heartrate_median,~ ] = heartrate( fusedsignal,AO_i_raw,fs);
% [ HR_ekg_avg ] = heartrate( ECG,ekg_locs,fs)

end
 













