function [vecout, tmps] = peakdetfilter(vecin,fs,iirmode)
%vecout = peakdetfilter(vecin,fs)

if nargin<2;
    fs = 200;
    iirmode = 0;
end

if nargin<3;
    iirmode = 0;
end

medfiltlen = round(fs/4*3);
%
% vecin = simplefilter(vecin,fs,1);
% % vecin = simplefilter(vecin,fs,0);
% vecout = SMQT(vecin,1,8);
%
% vecout = (vecout-medfilt1(vecout,medfiltlen));
% vecout(vecout<0)=0;


Norig = numel(vecin);

tmpdata = vecin;

tmps.rawdata = tmpdata;

%SIMPLEFILTERBLOCK
if iirmode == 0
    nfir = round(fs/50);
    nfir = max(nfir,2);
    nfir2 = nfir*2;
    bhp = -ones(1,nfir)/nfir;
    bhp(1) = bhp(1)+1;
    blp = triang(nfir2);

    tmp = filter(bhp,1,tmpdata);
    tmp = filter(blp,1,tmp);
    tmp(1:10)=0;
    tmp(1:3)=[];
    tmp(end+1:end+3) = 0;
else
%     wbut = [0.05,0.3]*fs/200;
    wbut = [0.05,0.3];   %passband = 5 - 30 Hz
    [bbut,abut] = butter(2,wbut);
    tmp = filter(bbut,abut,tmpdata);
    tmp(1:10)=0;
    tmp(1:3)=[];
    tmp(end+1:end+3) = 0;
end
tmps.bpfiltdata = tmp;

lflen = round(fs/2+1);
% lflen = 40;
tmp(end+1:end+lflen)=0;
tmp = tmp.^2;
sfout = filter(triang(lflen),1,abs(tmp));
sfout(1:floor(lflen/2))=[];
tmps.triangfiltdata = sfout;

%SMQT

vecout = SMQT(sfout,1,8);

tmps.smqtdata = vecout;
%MEDIANREMOVAL
mftmp = medfilt1(vecout,medfiltlen+1);
%     mftmp(1) = []; mftmp(end+1) = 0;
vecout = (vecout-mftmp);
vecout(vecout<0)=0;

vecout = vecout(1:Norig);
%     figure(5)
%     plot(vecout)
tmps.finaldata = vecout;


%calculating peaklocs
try
    posind = zeros(size(vecout));
    posind(vecout>0)=1;
    pidiff = diff(posind);
    peakstart = find(pidiff==1);
    peakend = find(pidiff==-1);
    if posind(1)>0
        peakstart = [1; peakstart];
    end

    if posind(end)>0
        peakend(end+1) = numel(posind);
    end

    tmps.peakstart = peakstart;
    tmps.peakend = peakend;
catch tt
    %     keyboard
end
