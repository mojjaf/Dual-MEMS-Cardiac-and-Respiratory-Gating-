function [ploc,pval,filteredsignal] = findheartbeats(datain,fs,fignum)
% [ploc,pval,filteredsignal] = findheartbeats(datain,fs,fignum)
if nargin<2
    fs = 200 ;
end

filteredsignal = peakdetfilter(datain,fs);
filteredsignal = (filteredsignal(:))' + (1:numel(filteredsignal))/numel(filteredsignal)/100;
filteredsignal = filteredsignal/max(filteredsignal);

MPH = .3;
MPD = round(fs/3);
try
%[pval, ploc] = findpeaksold(filteredsignal,'MINPEAKHEIGHT',MPH,'MINPEAKDISTANCE',MPD);
[pval, ploc] = findpeaks(filteredsignal,'MINPEAKHEIGHT',MPH,'MINPEAKDISTANCE',MPD);
catch 
    keyboard
end

if nargin > 2
    figure(fignum);
    plot(filteredsignal)
    hold on
    plot(ploc,pval,'r*')
    hold off
end
