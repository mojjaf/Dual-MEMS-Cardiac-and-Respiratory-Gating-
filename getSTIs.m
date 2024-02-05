function [heartRate,systolicTimeInterval] = getSTIs(original_signal,windowsize,Fs)



homomorphic_envelope = Homomorphic_Envelope_with_Hilbert(original_signal, Fs);

% windowsize=Fs*20;
% 

%% Find any samples outside of a integer number of windows:
WindowNumber = mod(length(homomorphic_envelope), windowsize);

%% Reshape the signal into a number of windows:
sampleframes = reshape(homomorphic_envelope(1:end-WindowNumber), windowsize, []);


%% Find the autocorrelation:
% y=homomorphic_envelope-mean(homomorphic_envelope,2);
y=bsxfun(@minus, sampleframes, mean(sampleframes, 2));
a2=mat2cell(y',ones(1,size(y,2)),size(y,1));
b2=cellfun(@xcorr,a2,'uniformoutput',false);
b2=cell2mat(b2);
% [c] = xcorr(y,'coeff');
b2=b2';
signal_autocorrelation = b2(length(sampleframes)+1:end,:)';

min_index = 0.5*Fs;
max_index = 2*Fs;
signal_acc=signal_autocorrelation(:,min_index:max_index);
[~,index] = max(signal_acc,[],2);

true_index = index+min_index-1;

heartRate = 60./(true_index/Fs);


%% Find the systolic time interval:
% From Schmidt: "The systolic duration is defined as the time from lag zero
% to the highest peak in the interval between 200 ms and half of the heart
% cycle duration"


max_sys_duration = round(((60./heartRate)*Fs)/2);
min_sys_duration = round(0.2*Fs);

[~,pos] = max(signal_acc(:,min_sys_duration:max_sys_duration),[],2);
systolicTimeInterval = (min_sys_duration+pos)/Fs;

figures=0;
if(figures)
    figure('Name', 'Heart rate calculation figure');
    plot(signal_autocorrelation);
    hold on;
    plot(true_index, signal_autocorrelation(true_index),'ro');
    plot((min_sys_duration+pos), signal_autocorrelation((min_sys_duration+pos)), 'mo');
    xlabel('Samples');
    legend('Autocorrelation', 'Position of max peak used to calculate HR', 'Position of max peak within systolic interval');

end

