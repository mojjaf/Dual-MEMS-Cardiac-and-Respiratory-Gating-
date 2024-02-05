function [ listmode_data_5bin ] = PhaseBasedGating(minlocs,inputsignal)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%% PhaseBasedGating
% minlocs=sort(minlocs);
% maxlocs=[];
% for i=1:numel(minlocs)-1
%         start =  minlocs(i);
%         stop = minlocs(i+1);
%         sample = inputsignal(start:stop);
%         
%      for j=1:numel(sample)
%            if sample(j)==max(inputsignal(start:stop))
%            maxlocs(i)=start+j-1;
%            end
%      end
% end    
    fs=25;
  Resp_cycle_length=[];
MEMS_resp_bin=[];
bin_num=[];
format short
for i=1:numel(minlocs)-1
Resp_cycle_length(i)=minlocs(i+1)-minlocs(i);
MEMS_resp_bin(i,1)=minlocs(i);
bin_num(i,1)=1;
for j=2:5
MEMS_resp_bin(i,j)=MEMS_resp_bin(i,j-1)+Resp_cycle_length(i)/5;
bin_num(i,j)=j;
end
gating_bins = reshape( MEMS_resp_bin.' ,1,numel(MEMS_resp_bin));
bin_nums = reshape( bin_num.' ,1,numel(bin_num));
listmode_data=[round(1000*gating_bins')/fs bin_nums'];
end
% listmode_data_5bin=listmode_data;
% 
% listmode_data(listmode_data(:,:)==3)=2;
% listmode_data(listmode_data(:,:)==4)=3;
% listmode_data(listmode_data(:,:)==5)=3;
% listmode_data(listmode_data(:,:)==6)=4;
% listmode_data(listmode_data(:,:)==7)=4;
% listmode_data(listmode_data(:,:)==8)=5;
% listmode_data(listmode_data(:,:)==9)=5;
% listmode_data(listmode_data(:,:)==10)=1;

listmode_data= int64(listmode_data);
listmode_data_5bin=int64(listmode_data);


figure
plot(inputsignal-mean(inputsignal));title('End-Expiratory Cycle Segmentation');axis tight;
line(repmat(gating_bins,[2 1]),repmat([min(inputsignal-mean(inputsignal))/2; max(inputsignal-mean(inputsignal))/2],size(gating_bins)),'LineWidth',2.5,'LineStyle','-.','Color','g');
line(repmat(minlocs,[2 1]),repmat([min(inputsignal-mean(inputsignal))/2; max(inputsignal-mean(inputsignal))/2],size(minlocs)),'LineWidth',2.5,'LineStyle','-.','Color','r');


end
