close all
clear all

cd

%%%%% HEALTHY GROUP
% datapath = 'D:\BIOSIGNAL_PROCESSING\MEDICAL IMAGING PROJECT\PilotStudy_2015\RPM and MEMS data (HealthySubjects)\MEMS/'
% datapath_rpm = 'D:\BIOSIGNAL_PROCESSING\MEDICAL IMAGING PROJECT\PilotStudy_2015\RPM and MEMS data (HealthySubjects)\RPM/'


% %%%%%%%%%%%%%%%% Patient group
datapath = 'D:\BIOSIGNAL_PROCESSING\MEDICAL IMAGING PROJECT\MinMotion_Micromachine paper project\RPM and MEMS data (HealthySubjects)\PLAQUA Study patients\MEMS/'

filepattern=fullfile(datapath,'*.txt')
txtfiles=dir(filepattern)

for k=1:length(txtfiles)
    %% loading RPM Signals
    
    
    %%
    data=load2_txt([datapath   txtfiles(k).name])
    fs=800;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         ICA TRANSFORM --- SENSOR FUSION ---
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   gr=0
    [AO_amp_ica,AO_i_ica,ekg_locs]=peakdetector_ICA(data,fs,gr);
    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         HILBERT TRANSFORM --- BINARY COMBINATION ---
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    coef=[1 1];
    [AO_amp_acc,AO_i_acc,pk,ekg_locs_hil]=Adpeakdet_hilbert(data,fs,coef,1,gr,0);%accelerometer

    [AO_amp_gyr,AO_i_gyr,pk,ekg_locs_hil]=Adpeakdet_hilbert(data,fs,coef,2,gr,0); %gyroscope

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        STATISTICAL ANALYSIS --- ALL METHODS ---
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fs_ds=100;
    [Precision_ICA(k,:),Sensitivity_ICA(k,:),Fscore_ICA(k,:),RMSE_ica(k,:)] = statistical_evaluation( ekg_locs', AO_i_ica,fs_ds )
    [Precision_ACC(k,:),Sensitivity_ACC(k,:),Fscore_ACC(k,:),RMSE_acc(k,:)] = statistical_evaluation( ekg_locs_hil', AO_i_acc,fs )
    [Precision_GYR(k,:),Sensitivity_GYR(k,:),Fscore_GYR(k,:),RMSE_gyr(k,:)] = statistical_evaluation( ekg_locs_hil', AO_i_gyr,fs )
    
    count(k,:)=k;
    format short
    T_evaluations = table(count,Precision_ICA,Sensitivity_ICA,Fscore_ICA,RMSE_ica,...
        Precision_ACC,Sensitivity_ACC,Fscore_ACC,RMSE_acc,...
        Precision_GYR,Sensitivity_GYR,Fscore_GYR,RMSE_gyr,...
    'VariableNames',{'subjectID','PR_ica','SE_ica','F1_ica','RMSE_ica','PR_acc','SE_acc','F1_acc','RMSE_acc','PR_gyr','SE_gyr','F1_gyr','RMSE_gyr'})
    close all
    
    
    
end