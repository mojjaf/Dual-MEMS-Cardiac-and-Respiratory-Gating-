function [ listmode_data_AMP ] = RPM_Amplitude_GATER( inputsignal,minima,maxima)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% inputsignal is the RPM signal and the minima and maxima are the locations
% of the local maximum and minimum points, respectively. M_ee and M_ei are
% the peak amplitudes at the end expiration and end inspiration.


%%  QualityIndex
 fs=25;
%  rmslen = round(fs/10);
% % 


%%
modality=0
[ M_ei,M_ee ] = gating_preprocessor( inputsignal,minima, maxima, modality );

% M_ee=inputsignal(minima);
% M_ei=inputsignal(maxima);
SD_ee=std(M_ee);
SD_ei=std(M_ei);
N=5;
UL=[];
LL=[];

numpeaks=min(numel(M_ee),numel(M_ei));

for i=1:N
    for j=1:numpeaks
        UL(i,j)=(M_ee(j)-SD_ee)+(N-i+1)*(M_ei(j)+SD_ei-(M_ee(j)-SD_ee))/N;
        LL(i,j)=(M_ee(j)-SD_ee)+(N-i)*(M_ei(j)+SD_ei-(M_ee(j)-SD_ee))/N;
    end
end




len=numel(inputsignal);


figure
plot(inputsignal)
hold on,scatter(minima,inputsignal(minima),'m');
hold on,scatter(maxima,inputsignal(maxima),'c');
%  line('XData', [0 len], 'YData', [mean(UL(1,:)) mean(UL(1,:))], 'LineStyle', '--', ...
%     'LineWidth', 2, 'Color','g')
 line('XData', [0 len], 'YData', [mean(UL(2,:)) mean(UL(2,:))], 'LineStyle', '--', ...
    'LineWidth', 2, 'Color','g')
text([0 len],[mean(UL(1,:)) mean(UL(1,:))],'BIN 1','Color','red','FontSize',14)
 line('XData', [0 len], 'YData', [mean(UL(3,:)) mean(UL(3,:))], 'LineStyle', '--', ...
    'LineWidth', 2, 'Color','m')
text([0 len],[mean(LL(2,:)) mean(LL(2,:))],'BIN 2','Color','red','FontSize',14)
 line('XData', [0 len], 'YData', [mean(UL(4,:)) mean(UL(4,:))], 'LineStyle', '--', ...
    'LineWidth', 2, 'Color','m')
text([0 len],[mean(LL(3,:)) mean(LL(3,:))],'BIN 3','Color','red','FontSize',14)
 line('XData', [0 len], 'YData', [mean(UL(5,:)) mean(UL(5,:))], 'LineStyle', '--', ...
    'LineWidth', 2, 'Color','g')
text([0 len],[mean(LL(4,:)) mean(LL(4,:))],'BIN 4','Color','red','FontSize',14)
%  line('XData', [0 len], 'YData', [mean(LL(5,:)) mean(LL(5,:))], 'LineStyle', '--', ...
%     'LineWidth', 2, 'Color','m')
text([0 len],[mean(LL(5,:)) mean(LL(5,:))],'BIN 5','Color','red','FontSize',14)
hold off


% figure (444)
% plot(inputsignal)
% hold on,scatter(minima,inputsignal(minima),'m');
%  hold on,scatter(maxima,inputsignal(maxima),'c');
% line(repmat(maxima,[2 1]),repmat([min(inputsignal-mean(inputsignal))/2; max(inputsignal-mean(inputsignal))/2],size(maxima)),'LineWidth',2,'LineStyle','-.','Color','r');
% line(repmat(minima,[2 1]),repmat([min(inputsignal-mean(inputsignal))/2; max(inputsignal-mean(inputsignal))/2],size(minima)),'LineWidth',2,'LineStyle','-.','Color','g');
% bin_num=[];



for i=1:len;
% if mean(UL(2,:))<=inputsignal(i) && inputsignal(i)<=mean(UL(1,:))
if mean(UL(2,:))<=inputsignal(i)
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
listmode_data_AMP=[1000*start_points/fs,gate_start(kk)];

end

