clear all, close all;
data_folder = 'D:\Seafile\Minmotion\Data\DDG_example\SSRB_C_unlister_01_June_2016\SSRB_C_unlister_01_June_2016/'; %SSRB_Eero, SSRB_22_Nov_2016
tvec_struct = load(strcat(data_folder,'tvec_60s.mat')); %60 seconds from the offset of PET imaging
tvec = tvec_struct.tvec;
cvec_struct = load(strcat(data_folder,'cvec_60s.mat')); %60 seconds from the offset of PET imaging
cvec_all = cvec_struct.cvec;

no_slices = size(cvec_all,1);

cvec = zeros(1,size(cvec_all,2)); %this will contain the center-of-mass values

k = 1; %exponent of signal (k = 2 corresponds to signal power)
p = 0.001; %attenuation coefficient

%In the following, the center of gravity of slice activities is computed
%and used as the DDG signal. To mitigate border effects, the end-slices are
%attenuated w.r.t. middle slices.

non_attenuated_slice_no = 10; %First slice that has no attenuated effect for the DDG signal
alpha = -log(p)/(non_attenuated_slice_no^2); %attenuation exponent parameter
total_counts = zeros(size(cvec));
for i = 1:no_slices
    attenuation_coeff = 1.;
    if i < non_attenuated_slice_no
        attenuation_distance = non_attenuated_slice_no-i;
        attenuation_coeff = exp(-alpha*attenuation_distance^2);
    end
    if i > no_slices-non_attenuated_slice_no+1
        attenuation_distance = i-(no_slices-non_attenuated_slice_no+1);
        attenuation_coeff = exp(-alpha*attenuation_distance^2);
    end
    
    total_counts = total_counts+attenuation_coeff*(cvec_all(i,:).^k);
    cvec = cvec+i*attenuation_coeff*(cvec_all(i,:).^k);
end
cvec = -cvec./total_counts;

load(strcat(data_folder,'RPM_reference'));

dt = tvec(2)-tvec(1);

%FFT:
cvec_mean = mean(cvec);
cvec_zero_mean = cvec-cvec_mean;


Fs = 1/dt; %Sampling frequency
L = length(tvec); %Length of signal
L2 = floor(L/2); %Halflength of signal
Y = fft(cvec_zero_mean);
P2 = abs(Y/L); %Two-sided spectrum P2
P1 = P2(1:L2+1);
P1(2:end-1) = 2*P1(2:end-1); %Single-sided spectrum P1

figure()
plot(tvec,total_counts)
title('Total counts')
xlabel('t (s)')
ylabel('Events per 50 ms')

f = Fs*(0:L2)/L;
figure()
plot(f,P1)
title('Single-Sided Amplitude Spectrum of counts(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

%Smooth FFT:
f0 = 0.5;
f1 = 1.0;
p1 = 0.001;
[slp,shp] = smooth_fft_filter(tvec,cvec,f0,f1-f0,p1);
title_txt = strcat(['Filtered signal, with k = ', num2str(k), ' and cut-off freq of ', num2str(f0), '-', num2str(f1), ' Hz']);

figure()
plot(tvec,slp)
title(title_txt)
xlabel('t (s)')
ylabel('Center of mass (negative)')

figure()
plot(tvec,shp)
title('High frequency part of the Fourier filtered signal')
xlabel('t (s)')
ylabel('Center of mass (negative)')


%RPM: find the starting point of PET
t_offset = 50*1e3; %50 seconds from the beginning of PET
start_of_PET=min(find(RPMDataAfterSync.time>=t_offset)); % find first index
time_limit=60; %s 
samples_RPM=25*time_limit; % 25 Hz sampling

tvec_RPM = (RPMDataAfterSync.time(start_of_PET:start_of_PET+samples_RPM)-t_offset)/1e3;
RPM_signal = -RPMDataAfterSync.amp(start_of_PET:start_of_PET+samples_RPM);
RPM_signal_normalized = (RPM_signal-mean(RPM_signal))/std(RPM_signal); %normalized to unit variance

%For the smooth Fourier filtering:
slp_normalized = (slp-cvec_mean)/std(slp); %normalized to unit variance
sgol_normalized = DDG_smoother( total_counts,Fs, 4, 71 );
figure()
hold all
plot(tvec,slp_normalized)
plot(tvec,sgol_normalized)
plot(tvec_RPM,RPM_signal_normalized)
hold off
title('Fourier filtered PET center of mass vs. RPM')
xlabel('t (s)')
ylabel('Arbitrary units')
legend('DDG FFT filter','DDG Savitzky-Golay', 'RPM')