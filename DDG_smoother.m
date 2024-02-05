function [ respiration ] = DDG_smoother( signal,fs,order, framelen )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

sgf = sgolayfilt(signal,order, framelen);
sgf_filt  = baseline_removal_imu(sgf,0.05,fs);
respiration=zscore(sgf_filt);



end

