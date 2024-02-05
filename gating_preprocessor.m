function [ M_ei,M_ee ] = gating_preprocessor( inputsignal,minima, maxima, modality )
%UNTITLED5 Summary of this function goes here
%   This function bounds the corresponding end expiration and end inspiration points and
%   estimates the missing ones. Modality 1 is MEMS and 0 is RPM.

%%%%%%%%%%%%%%%%%%%%%%
stps=numel(inputsignal)/numel(minima);
xq=1:stps:numel(inputsignal);

qx=double(xq);

V=inputsignal(maxima);
X=double(maxima);
vq2 = interp1(X,V,qx,'previous','extrap');

nanx = isnan(vq2);
vq2(nanx)=inputsignal(maxima(nanx));

for i=1:numel(vq2)
    maxlocs(i)=find(inputsignal==vq2(i));
end
%%%

maxima=maxlocs;
%%%
if modality

% 
P2P_thrsh=median(inputsignal(maxima)-inputsignal(minima));
badlocs_max=find(inputsignal(maxima)>2.5*P2P_thrsh);
badlocs_min=find(-inputsignal(minima)>2.5*P2P_thrsh);
maxima(badlocs_max)=[];
minima(badlocs_min)=[];
end

maxima_fixed=sort(maxima);
minima_fixed=sort(minima);

M_ei=inputsignal(maxima_fixed);
M_ee=inputsignal(minima_fixed);



end

