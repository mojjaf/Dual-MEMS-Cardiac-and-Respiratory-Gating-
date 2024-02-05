function M = SMQT(V, l, L)
if l>L 
    M = zeros(size(V), 'like', V);
    return;
end
meanV = nanmean(V(:));
D0 = V;
D1 = V;
if not(isnan(meanV)) 
    D0(D0 > meanV) = NaN;
    D1(D1 <= meanV) = NaN;
end
M = not(isnan(D1)) * (2^(L-l));
if l==L 
    return;
end
M0 = SMQT(D0, l+1, L);
M1 = SMQT(D1, l+1, L);
M = M + M0 + M1;
end