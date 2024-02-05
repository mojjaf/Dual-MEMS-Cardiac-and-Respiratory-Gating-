close all
clear all
cd
% datapath = 'D:\BIOSIGNAL_PROCESSING\Physiological Recordings\CoronaryDiseasedPatients/'
datapath = 'D:\BIOSIGNAL_PROCESSING\MEDICAL IMAGING PROJECT\MIC2018 project\Micromachines project/'
% datapath='D:\BIOSIGNAL_PROCESSING\MEDICAL IMAGING PROJECT\PilotStudy_2015\RPM and MEMS data (HealthySubjects)/'
% datapath='D:\BIOSIGNAL_PROCESSING\Atrial Fibrillation\Clinical Data_Confidential\AF measurements may-june 2014/'
filepattern=fullfile(datapath,'*.txt')
txtfiles=dir(filepattern)
 


for k=1:length(txtfiles)

 
data=load2_txt([datapath   txtfiles(k).name])
fs=800;  
gyroY=double(data.gyroY)*(250/32767);
gyroX=double(data.gyroX)*(250/32767);
gyroZ=double(data.gyroZ)*(250/32767);
accX=double(data.accX)*(2/32760);
accY=double(data.accY)*(2/32760);
accZ=double(data.accZ)*(2/32760);


Hp=0.1;
GyroX=baseline_removal_imu(gyroX,Hp,fs);
GyroY=baseline_removal_imu(gyroY,Hp,fs);
GyroZ=baseline_removal_imu(gyroZ,Hp,fs);
AccX=baseline_removal_imu(accX,Hp,fs);
AccY=baseline_removal_imu(accY,Hp,fs);
AccZ=baseline_removal_imu(accZ,Hp,fs);
%%%%% IN CASE WE NEEDED TO BANDPASS SIGNAL   
%   hp_G=1;
%   lp_G=20;
% GyroY=fft_filter(gyroY,800,hp_G,lp_G);
% GyroX=fft_filter(gyroX,800,hp_G,lp_G);
% GyroZ=fft_filter(gyroZ,800,hp_G,lp_G);
% % 
% hp_A=1;
% lp_A=40;
% AccZ=fft_filter(accZ,800,hp_A,lp_A);
% AccY=fft_filter(accY,800,hp_A,lp_A);
% AccX=fft_filter(accX,800,hp_A,lp_A);
  %%%%%%%%%%%%%%%%%%%%%%%%FFT ANALYSIS


SCG_X=AccX(10000:140000); %%%% TAKE A SAMPLE PORTION OF THE SIGNALS
SCG_Y=AccY(10000:140000);
SCG_Z=AccZ(10000:140000);


GCG_X=GyroX(10000:140000);
GCG_Y=GyroY(10000:140000);
GCG_Z=GyroZ(10000:140000);

m = length(GCG_X);       % original sample length
n = pow2(nextpow2(m));  % transform length
y = fft(GCG_X,n);
f = (0:n-1)*(fs/n);
power = abs(y).^2/n;
figure;
ax(1)=subplot(321);plot(f(1:floor(n/2)),power(1:floor(n/2)))
title('Spectral Analysis of GCG X axis')
xlabel('Frequency')
ylabel('Power')
m = length(GCG_Y);       % original sample length
n = pow2(nextpow2(m));  % transform length
y = fft(GCG_Y,n);
f = (0:n-1)*(fs/n);
power = abs(y).^2/n;
ax(2)=subplot(323);plot(f(1:floor(n/2)),power(1:floor(n/2)))
title('Spectral Analysis of GCG Y axis')
xlabel('Frequency')
ylabel('Power')
m = length(GCG_Z);       % original sample length
n = pow2(nextpow2(m));  % transform length
y = fft(GCG_Z,n);
f = (0:n-1)*(fs/n);
power = abs(y).^2/n;
ax(3)=subplot(325);plot(f(1:floor(n/2)),power(1:floor(n/2)))
title('Spectral Analysis of GCG Y axis')
xlabel('Frequency')
ylabel('Power')



    m = length(SCG_X);       % original sample length
n = pow2(nextpow2(m));  % transform length
y = fft(SCG_X,n);
f = (0:n-1)*(fs/n);
power = abs(y).^2/n;      
ax(4)=subplot(322);plot(f(1:floor(n/2)),power(1:floor(n/2)))
title('Spectral Analysis of SCG X axis')
xlabel('Frequency')
ylabel('Power')

 m = length(SCG_Y);       % original sample length
n = pow2(nextpow2(m));  % transform length
y = fft(SCG_Y,n);
f = (0:n-1)*(fs/n);
power = abs(y).^2/n;      
ax(5)=subplot(324);plot(f(1:floor(n/2)),power(1:floor(n/2)))
title('Spectral Analysis of SCG Y axis')
xlabel('Frequency')
ylabel('Power')

 m = length(SCG_Z);       % original sample length
n = pow2(nextpow2(m));  % transform length
y = fft(SCG_Z,n);
f = (0:n-1)*(fs/n);
power = abs(y).^2/n;      
ax(6)=subplot(326);plot(f(1:floor(n/2)),power(1:floor(n/2)))
title('Spectral Analysis of SCG Z axis')
xlabel('Frequency')
ylabel('Power')
linkaxes([ax(6) ax(5) ax(4) ax(3) ax(2) ax(1)],'x');


  
   gr=1; %%% PLOTING FLAG
if gr
figure

ax(1)=subplot(3,2,1);plot(GyroX)
xlabel('X axis')
ylabel('angular Velocity')
ax(2)=subplot(3,2,3);plot(GyroY)
xlabel('Y axis')
ylabel('angular Velocity')
ax(3)=subplot(3,2,5);plot(GyroZ)
xlabel('Z axis')
ylabel('angular Velocity')
ax(4)=subplot(3,2,2);plot(AccX)
xlabel('X axis')
ylabel('accelerations')
ax(5)=subplot(3,2,4);plot(AccY)
xlabel('Y axis')
ylabel('accelerations')
ax(6)=subplot(3,2,6);plot(AccZ)
xlabel('Z axis')
ylabel('accelerations')
linkaxes([ax(6) ax(5) ax(4) ax(3) ax(2) ax(1)],'x');
end
%%
%%%% respiration extraction from ACCELEROMETER X and Y using ARCTAN
thetaA=atan(AccY./(sqrt(AccX.^2+AccZ.^2)));
thetaB=atan(AccX./(sqrt(AccY.^2+AccZ.^2)));

 ThetaA=fft_filter(thetaA,fs,.1,2);
 ThetaB=fft_filter(thetaB,fs,.1,2);
%   ThetaB=medfilt1(ThetaB,7);

 %%%%%%%%%%%%%%%%%%%%%%%  respiration extraction from ACCELEROMETER X and Y
 %%%%%%%%%%%%%%%%%%%%%%%  WITHOUT ARCTAN (APPROXIMATION)
  thetaA_aprx=(AccY./(sqrt(AccX.^2+AccZ.^2)));
  thetaB_aprx=(AccX./(sqrt(AccY.^2+AccZ.^2)));
  ThetaA_aprx=fft_filter(thetaA_aprx,fs,.1,2);
  ThetaB_aprx=fft_filter(thetaB_aprx,fs,.1,2);
 %%%%%%%%%%%%%%%%%%%%%%%%% respiration extraction from ACCELEROMETER X and
 %%%%%%%%%%%%%%%%%%%%%%%%% Y using ARCSIN
 alphaA=asin(AccY);
 alphaB=asin(AccX);
 AlphaA=fft_filter(alphaA,fs,.1,2);
 AlphaB=fft_filter(alphaB,fs,.1,2);
 
 AlphaB_aprx=fft_filter(AccX,fs,.1,2);
 AlphaA_aprx=fft_filter(AccY,fs,.1,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% respiration extraction from GYROSCOPE X and
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Y by INTEGRATING (convering angular velocity to
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% angular displacement)
samples=1:numel(AccX);
t=samples/fs;
rotX=cumtrapz(t,GyroX);
  RotX=fft_filter(rotX*0.017,fs,.1,2);

rotY=cumtrapz(t,GyroY);
  RotY=fft_filter(rotY*0.017,fs,.1,2);


  
 %%%%% ADR signals from Acc X and Y
 figure
 title('Accelerometeric Derived Respiration Rotation')
ax(1)=subplot(2,1,1);plot(ThetaA)
xlabel('X axis')
ylabel('angular Displacement')
ax(2)=subplot(2,1,2);plot(ThetaB)
ylabel('angular Displacement')
xlabel('Y axis')
linkaxes([ax(2) ax(1)],'x');


%% Principal component analysis (PCA) 
temp=[RotY,RotX,ThetaA,ThetaA_aprx,ThetaB_aprx,ThetaB,AlphaA,AlphaA_aprx,AlphaB,AlphaB_aprx];
 resp_pca=pca(temp');
  resp_pca=fft_filter(resp_pca,fs,.1,2);

figure
 title('MEMS Derived Respirations')
ax(1)=subplot(6,1,1);plot(RotY)
ylabel('Angular Disp (rad)')
xlabel('gyroscope')
ax(2)=subplot(6,1,2);plot(ThetaA)
ylabel('Angular Disp (rad) ')
xlabel('Theta ADR (arctang)')
ax(3)=subplot(6,1,3);plot(ThetaA_aprx)
ylabel('Angular Disp (rad)')
xlabel('Theta ADR (approximation)')
ax(4)=subplot(6,1,4);plot(AlphaB)
ylabel('Angular Dispt (rad)')
xlabel('Theta ADR (arcsin)')
ax(5)=subplot(6,1,5);plot(AlphaB_aprx)
ylabel('Acceleration (milli-G)')
xlabel('Theta ADR (approximation)')

ax(6)=subplot(6,1,6);plot(resp_pca(:,1))
ylabel('Angular Dispt (rad)')
xlabel('Respiration PCA')

linkaxes([ax(6) ax(5) ax(4) ax(3) ax(2) ax(1)],'x');




 figure
plot(RotY);
hold on
plot(ThetaB*100,'r')


figure
title('Accelerometeric Derived Respiration Rotation')
ax(1)=subplot(2,1,1);plot(AlphaA)
xlabel('Alpha A Y axis')
ylabel('angular Displacement')
ax(2)=subplot(2,1,2);plot(AlphaB)
ylabel('angular Displacement')
xlabel('Alpha B X axis')
linkaxes([ax(2) ax(1)],'x');
figure
 title('Gyroscopic Derived Respiration Rotation')

ax(1)=subplot(2,1,1);plot(RotX)
xlabel('X axis')
ylabel('angular Displacement')
ax(2)=subplot(2,1,2);plot(RotY)
ylabel('angular Displacement')
xlabel('Y axis')
linkaxes([ax(2) ax(1)],'x');



%% loading RPM Signals
 load( 'plaque2_sub2_23032016' )
 
%  [Data_layoutamplitude,phase,timestamp,validflag,ttlin,mark,ttlout] = importRPM('patient_01_06_2016_4DPET.vxp');

 RPM=(Data_layoutamplitude-mean(Data_layoutamplitude))/std(Data_layoutamplitude);
 RPM=fft_filter(Data_layoutamplitude,25,0.1,12);
 RPM= resample(RPM,32,1);
%  GDR=downsample(RotX,32);
%  ADR=downsample(ThetaB,32);


%%%%%%%%%%%% Synching MEMS derived Respiration with RPM
RPM=(RPM-mean(RPM))/std(RPM);
  GDR=-RotY;
  ADR=ThetaB;
  s7=-resp_pca(:,1);
  s4=ThetaB_aprx;
  s5=AlphaB;
  s6=AlphaB_aprx;
  s3=ADR;
  s2=GDR;
  s1=RPM;
 [x1,x2] = alignsignals(s1,s2);
  [x1,x3] = alignsignals(s1,s3);
  [x1,x4] = alignsignals(s1,s4);
  [x1,x5] = alignsignals(s1,s5);
  [x1,x6] = alignsignals(s1,s6);
  [x1,x7] = alignsignals(s1,s7);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT Sychnced signals
figure
plot(normalize(x1,1))
hold on
plot(normalize(x2,1))
hold on
plot(normalize(x3,1))
hold on
plot(normalize(x4,1))
hold on
plot(normalize(x5,1))
hold on
plot(normalize(x6,1))
hold on
plot(normalize(x7,1))
ylabel('Arbitrary Units (norm)')

legend('RPM','GDR (theta)','ADR(arctan)','ADR(approximate)','ADR(arcsin)','ADR(approximate)','PCA');

figure
plot(x2)
hold on
plot(x3)
hold on
plot(x4)
hold on
plot(x5)
hold on
plot(x6)
hold on
plot(x7)
ylabel('angular Displacement(rad)')

legend('GDR (theta)','ADR(arctan)','ADR(approximate)','ADR(arcsin)','ADR(approximate)','PCA');


pause
close all
end