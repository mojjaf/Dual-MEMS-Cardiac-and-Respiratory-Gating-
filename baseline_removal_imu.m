function [ corrected ] = baseline_removal_imu( inputsig,lowban, fs )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
freslut=fft(inputsig);
freslut(1:round(length(freslut)*lowban/fs))=0;
freslut(end-round(length(freslut)*lowban/fs): end)=0;
corrected=real(ifft(freslut));

end

