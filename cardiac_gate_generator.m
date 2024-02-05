function [ LM_5bin,LM_3bin ] = cardiac_gate_generator( signal, fd_points,fs,gr )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

cycle_length=[];
MEMS_bin=[];
format short
for i=1:numel(fd_points)-1
cycle_length(i)=fd_points(i+1)-fd_points(i);
MEMS_bin(i,1)=fd_points(i);
bin_num(i,1)=1;
for j=2:5
MEMS_bin(i,j)=MEMS_bin(i,j-1)+cycle_length(i)/10;

bin_num(i,j)=j;
end
gating_bins = reshape( MEMS_bin.' ,1,numel(MEMS_bin));
bin_nums = reshape( bin_num.' ,1,numel(bin_num));
listmode_data=[round(1000*gating_bins')/fs bin_nums'];
end
listmode_data_5bin=listmode_data;
listmode_data(listmode_data(:,:)==4)=3;
listmode_data(listmode_data(:,:)==5)=3;
listmode_data= int64(listmode_data);
listmode_data_5bin=int64(listmode_data_5bin);
read=0;
 LM_5bin=table(listmode_data_5bin);
 LM_3bin=table(listmode_data);


if gr
figure
plot(signal-mean(signal));title('Pulse train of the found AO on MEMS signal');axis tight;
line(repmat(gating_bins,[2 1]),repmat([min(signal-mean(signal))/2; max(signal-mean(signal))/2],size(gating_bins)),'LineWidth',2.5,'LineStyle','-.','Color','g');
line(repmat(fd_points',[2 1]),repmat([min(signal-mean(signal))/2; max(signal-mean(signal))/2],size(fd_points')),'LineWidth',2.5,'LineStyle','-.','Color','r');
end

% [ ADR,GDR ] = MEMS_ADR_GDR( data)
% [bin, bin_num,Resp_listmode_data ] = respiratory_gatingbins(GDR,5);
% Resp_listmode_data= int64(Resp_listmode_data);
% read=0;
% if read
% 
% savdir1 ='D:\BIOSIGNAL PROCESSING\MEDICAL IMAGING PROJECT\PET_PLAQUA2 project\dualgating_bins' ;
% %  gyro_gatingbins=(listmode_data);
% %  save(fullfile(savdir1,sprintf('GyroGateBins_sub_%02d.txt',2)),'listmode_data','-ascii','-tabs');
%  
%  RespLM_gyroBindataT=table(Resp_listmode_data);
% writetable(RespLM_gyroBindataT,'RespiratoryLM_GYRO_5Bins_PET.txt','Delimiter','\t')
% end 
LM_CardiacData=listmode_data;

end

