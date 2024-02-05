
function [ HR_avg, HR_i_avg ] = heartrate( sig_in,locs,Fs)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
sig_len=numel(sig_in);

num=floor((sig_len-Fs)/Fs);
num_10=num/5;
% rr_int=diff(locs);
for ii=1:num_10-1
    rr_intval=Fs*5;
    indvec =( (1:rr_intval) + (ii-1)*Fs);
    HR(ii)=60*(1/(mean(diff(locs(locs>min(indvec)& locs<max(indvec))))/Fs));
end
HR_avg=median(HR);
% HR_std=std(HR);
HR_i_avg=HR;


end

