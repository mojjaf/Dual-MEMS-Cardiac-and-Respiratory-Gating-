function [ listmode_data_phase,listmode_data_AMP ] = MEMS_GATER( inputsignal,minima,maxima,fs )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% inputsignal=SMQT(inputsignal,1,8);
M_ee=inputsignal(minima);
M_ei=inputsignal(maxima);
SD_ee=std(M_ee);
SD_ei=std(M_ei);
N=5;
UL=[];
LL=[];
fs=fs/25;
rmslen = round(fs/2);


rmsthresh = rms(inputsignal,rmslen);
mrval = median(rmsthresh);
% qtind3 = find(rmsaccZ<2*mrval);


badind = find(rmsthresh>1.25*mrval)






for i=1:N
    for j=1:numel(minima)
        UL(i,j)=(M_ee(j)-SD_ee)+(N-i+1)*(M_ei(j)+SD_ei-(M_ee(j)-SD_ee))/N;
        LL(i,j)=(M_ee(j)-SD_ee)+(N-i)*(M_ei(j)+SD_ei-(M_ee(j)-SD_ee))/N;
    end
end

len=numel(inputsignal);

figure
plot(inputsignal)
hold on,scatter(minima,inputsignal(minima),'m');
hold on,scatter(maxima,inputsignal(maxima),'c');
hold on;plot(maxima,UL(1,:),'LineWidth',2,'Linestyle','-.','color','g');
hold on;plot(maxima,UL(2,:),'LineWidth',2,'Linestyle','-.','color','g');
hold on;plot(maxima,UL(3,:),'LineWidth',2,'Linestyle','-.','color','g');
hold on;plot(maxima,UL(4,:),'LineWidth',2,'Linestyle','-.','color','g');
hold on;plot(maxima,UL(5,:),'LineWidth',2,'Linestyle','-.','color','g');
hold on;plot(maxima,LL(1,:),'LineWidth',2,'Linestyle','-.','color','r');
hold on;plot(maxima,LL(2,:),'LineWidth',2,'Linestyle','-.','color','r');
hold on;plot(maxima,LL(3,:),'LineWidth',2,'Linestyle','-.','color','r');
hold on;plot(maxima,LL(4,:),'LineWidth',2,'Linestyle','-.','color','r');
hold on;plot(minima,LL(5,:),'LineWidth',2,'Linestyle','-.','color','r');
hold on,scatter(maxima,UL(1,:),'k');
hold on,scatter(maxima,UL(2,:),'k');
hold on,scatter(maxima,UL(3,:),'k');
hold on,scatter(maxima,UL(4,:),'k');
hold on,scatter(maxima,UL(5,:),'k');
hold on,scatter(minima,LL(5,:),'k');

 line('XData', [0 len], 'YData', [mean(UL(1,:)) mean(UL(1,:))], 'LineStyle', '--', ...
    'LineWidth', 2, 'Color','m')
 line('XData', [0 len], 'YData', [mean(UL(2,:)) mean(UL(2,:))], 'LineStyle', '--', ...
    'LineWidth', 2, 'Color','m')
 line('XData', [0 len], 'YData', [mean(UL(3,:)) mean(UL(3,:))], 'LineStyle', '--', ...
    'LineWidth', 2, 'Color','m')
 line('XData', [0 len], 'YData', [mean(UL(4,:)) mean(UL(4,:))], 'LineStyle', '--', ...
    'LineWidth', 2, 'Color','m')
 line('XData', [0 len], 'YData', [mean(UL(5,:)) mean(UL(5,:))], 'LineStyle', '--', ...
    'LineWidth', 2, 'Color','m')
 line('XData', [0 len], 'YData', [mean(LL(5,:)) mean(LL(5,:))], 'LineStyle', '--', ...
    'LineWidth', 2, 'Color','m')

figure
plot(inputsignal)
hold on,scatter(minima,inputsignal(minima),'m');
hold on,scatter(maxima,inputsignal(maxima),'c');
line(repmat(maxima,[2 1]),repmat([min(inputsignal-mean(inputsignal))/2; max(inputsignal-mean(inputsignal))/2],size(maxima)),'LineWidth',2,'LineStyle','-.','Color','r');
line(repmat(minima,[2 1]),repmat([min(inputsignal-mean(inputsignal))/2; max(inputsignal-mean(inputsignal))/2],size(minima)),'LineWidth',2,'LineStyle','-.','Color','g');
bin_num=[];

for i=1:len;
if mean(UL(2,:))<=inputsignal(i) && inputsignal(i)<=mean(UL(1,:))
bin_num(i)=1;
elseif mean(UL(3,:))<=inputsignal(i) && inputsignal(i)<=mean(UL(2,:))
bin_num(i)=2;
elseif mean(UL(4,:))<=inputsignal(i) && inputsignal(i)<=mean(UL(3,:))
bin_num(i)=3;
elseif mean(UL(5,:))<=inputsignal(i) && inputsignal(i)<=mean(UL(4,:))
bin_num(i)=4;
% elseif mean(LL(5,:))<=GDR(i) && GDR(i)<=mean(UL(5,:))
elseif inputsignal(i)<=mean(UL(5,:))
bin_num(i)=5;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% interpolation 


bin_num=bin_num';

[p1 o]=find(bin_num==1);
[p2 o]=find(bin_num==2);
[p3 o]=find(bin_num==3);
[p4 o]=find(bin_num==4);
[p5 o]=find(bin_num==5);
pd1=numel(p1)/numel(bin_num);
pd2=numel(p2)/numel(bin_num);
pd3=numel(p3)/numel(bin_num);
pd4=numel(p4)/numel(bin_num);
pd5=numel(p5)/numel(bin_num);

gate_start=diff(bin_num);
gate_start(1)=bin_num(1);
 for ii=2:numel(gate_start)
if gate_start(ii)~=0
gate_start(ii)=bin_num(ii+1);
end
end
[kk jj]=find(gate_start>0);

start_points=kk+1;
listmode_data_AMP=[1000*start_points/25,gate_start(kk)];
cycle_length=[];
% MEMS_phasebin=[];
% phasebin_num=[];
% for i=1:numel(maxima)-1
% cycle_length(i)=maxima(i+1)-maxima(i);
% MEMS_phasebin(i,1)=maxima(i);
% phasebin_num(i,1)=1;
% for j=2:5
% MEMS_phasebin(i,j)=MEMS_phasebin(i,j-1)+cycle_length(i)/5;
% phasebin_num(i,j)=j;
% end
% end
% gating_bins = reshape( MEMS_phasebin.' ,1,numel(MEMS_phasebin));
% phasebin_num = reshape( phasebin_num.' ,1,numel(phasebin_num));
% 
% listmode_data_phase=[round(1000*gating_bins')/25 phasebin_num'];   
figure
plot(inputsignal)
hold on,scatter(minima,inputsignal(minima),'m');
hold on,scatter(maxima,inputsignal(maxima),'c');
line(repmat(gating_bins,[2 1]),repmat([min(inputsignal-mean(inputsignal)); max(inputsignal-mean(inputsignal))],size(gating_bins)),'LineWidth',1,'LineStyle','-.','Color','g');
line(repmat(maxima,[2 1]),repmat([min(inputsignal-mean(inputsignal)); max(inputsignal-mean(inputsignal))],size(maxima)),'LineWidth',1,'LineStyle','-.','Color','r');
  line('XData', [0 len], 'YData', [mean(UL(1,:)) mean(UL(1,:))], 'LineStyle', '--', ...
    'LineWidth', 2, 'Color','m')
 line('XData', [0 len], 'YData', [mean(UL(2,:)) mean(UL(2,:))], 'LineStyle', '--', ...
    'LineWidth', 2, 'Color','m')
 line('XData', [0 len], 'YData', [mean(UL(3,:)) mean(UL(3,:))], 'LineStyle', '--', ...
    'LineWidth', 2, 'Color','m')
 line('XData', [0 len], 'YData', [mean(UL(4,:)) mean(UL(4,:))], 'LineStyle', '--', ...
    'LineWidth', 2, 'Color','m')
 line('XData', [0 len], 'YData', [mean(UL(5,:)) mean(UL(5,:))], 'LineStyle', '--', ...
    'LineWidth', 2, 'Color','m')
 line('XData', [0 len], 'YData', [mean(LL(5,:)) mean(LL(5,:))], 'LineStyle', '--', ...
    'LineWidth', 2, 'Color','m')
end

