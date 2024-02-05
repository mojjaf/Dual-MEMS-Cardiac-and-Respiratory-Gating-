function [all_traces,avpeak,ploc] = averagedpeaks_ECG_REF(datain,fs,fignum, ploc)
% [avpeak,ploc] = averagedpeaks(signalin,fs,fignum)

if nargin<2
    fs = 200;
end
%UUOLLI
%   figure(fignum)
%UUOLLI

outputlen = 1000;

olli_avvi_median(1,1:outputlen)=0;

%ploc = findheartbeats(datain,fs);

num = numel(ploc)-1;

% figure(1)
% plot(datain)
% hold on
% plot(ploc,datain(ploc),'r*')
% hold off
% figure(2)
avpeak = zeros(1,outputlen);
for ii = 1 : num
%     tmp = datain(ploc(ii):ploc(ii+1));
tmp = datain(ploc(ii):ploc(ii+1));
    len = numel(tmp);
    indvec = ceil((1:outputlen)/outputlen*len);
%     keyboard 
    
    tmp2 = tmp(indvec);
%     plot(tmp2)
%     pause(.02)
%     keyboard
    avpeak = avpeak + tmp2(:)';
    
    %plot(tmp2(:)'); hold on;            
    
    olli_avvi_median(ii,1:outputlen)=tmp2(:)';
end
    

% if nargin>2
% %    figure(fignum)
% %    plot(avpeak)           
% 
% %plot(median(olli_avvi_median,1),'Color','b','LineWidth',3);
% 
% % hold off;
% end
all_traces=olli_avvi_median;
avpeak=(median(olli_avvi_median,1));
% avpeak=(olli_avvi_median);
