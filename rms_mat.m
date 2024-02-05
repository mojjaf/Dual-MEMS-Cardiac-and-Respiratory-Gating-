function vecout = rms_mat(vecin, len);
%function vecout = rms(vecin, len);
vecout=[];
%update 20170403:  changed filtering to conv, to get rid of filtering
%delay.
hlen = floor(len/2);
wlen = 2*hlen+1;
for i=1:size(vecin,2)
vecinsq = vecin(:,i).*vecin(:,i);
% vecinft = filter(ones(1,len),1,vecinsq)/len;
vecinft = conv(ones(1,wlen),vecinsq)/wlen;
vecinft(1:hlen) = [];
vecinft(end-hlen+1:end) = [];

vecout(:,i) = sqrt(vecinft);
end

end


