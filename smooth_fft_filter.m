function [slp,shp] = smooth_fft_filter(t,s,fmax,fbw,p) %t = time vector, s = signal
    %fmax = maximum frequency without attenuation
    %fbw = bandwidth where the frequencies are attenuated to p
    %p = smallest attenuation before zero attenuation
    %slp = low-pass filtered signal, shp = high-pass filtered signal

    s_mean = mean(s);
    s_zero_mean = s-s_mean; %zero-mean signal
    dt = t(2)-t(1); %assuming constant sampling rate

    Fs = 1/dt; %Sampling frequency
    L = length(t); %Length of signal
    L2 = floor(L/2); %Halflength of signal
    Y = fft(s_zero_mean);
    f = Fs*(0:L2)/L;

    fmax_ind = find(f>fmax,1); %finds the smallest index for which f(fmax_ind)>fmax
    Y_filtered = Y;

    alpha = -log(p)/(fbw^2); %attenuation coefficient for Gaussian-shaped attenuation
    fbw_ind = find(f>fmax+fbw,1);
    f_attenuated = f(fmax_ind:fbw_ind)-fmax;
    attenuation_factors = exp(-alpha*f_attenuated.^2); %Gaussian-shaped attenuation between fmax and fmax+fbw

    Y_filtered(fmax_ind:fbw_ind) = attenuation_factors.*Y_filtered(fmax_ind:fbw_ind);
    Y_filtered(L-fbw_ind+1:L-fmax_ind+1) = fliplr(attenuation_factors).*Y_filtered(L-fbw_ind+1:L-fmax_ind+1); %+1 due to Matlab indexing
    Y_filtered(fbw_ind:L-fbw_ind) = 0.;

    slp = real(ifft(Y_filtered))+s_mean;
    shp = s-slp;
end