function [avpeak,ploc] = averagedpeaks(datain,fs,fignum)
% [avpeak,ploc] = averagedpeaks(signalin,fs,fignum)

if nargin<2
    fs = 200;
end


outputlen = 1000;
ploc = findheartbeats(datain,fs);
num = numel(ploc)-1;

%olli_avvi_median(1,1:outputlen)=0;

% figure(1)
% plot(datain)
% hold on
% plot(ploc,datain(ploc),'r*')
% hold off
% figure(2)
avpeak = zeros(1,outputlen);
for ii = 1 : num
    tmp = datain(ploc(ii):ploc(ii+1));
    len = numel(tmp);
    indvec = ceil((1:outputlen)/outputlen*len);
%     keyboard 
    
    tmp2 = tmp(indvec);
%     plot(tmp2)
%     pause(.02)
%     keyboard
    avpeak = avpeak + tmp2(:)';
    
%    olli_avvi_median(ii,1:outputlen)=tmp2(:)';
    
end
    

% if nargin>2
%     figure(fignum)
%     plot(avpeak)    
end


