function [ minloc_pca, maxloc_pca, minloc_rpm, maxloc_rpm ] = find_peaks_rpm_guided( RPM, PCA, gr )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
% maxloc_pca = ampd(PCA);
% minloc_pca = ampd(-PCA);

maxloc_pca = [];
minloc_pca = [];

maxloc_rpm= ampd(RPM); % ampd is automated multiscale peak detection from Git
minloc_rpm= ampd(-RPM);

[resp_corr,~]=corr(RPM,PCA,'Type','Pearson')
if resp_corr<0
    PCA=-PCA;
end
    
for i=1:(numel(maxloc_rpm))
        peak_sample_n = maxloc_rpm(i); %in samples
        start = peak_sample_n-15; % 15 is the window size corresponding to 0.6 sec
        stop = peak_sample_n+15 ;
        sample_sig = PCA(start:stop);
        
     for j=1:numel(sample_sig)
         if sample_sig(j)==max(PCA(start:stop))
           maxloc_pca(i)=start+j-1;
         end
     end
       
end


for i=1:(numel(minloc_rpm))
        peak_sample_n = minloc_rpm(i); %in samples
        start = peak_sample_n-15;
        stop = peak_sample_n+15 ;
        sample_sig = PCA(start:stop);
        
     for j=1:numel(sample_sig)
         if sample_sig(j)==min(PCA(start:stop))
           minloc_pca(i)=start+j-1;
         end
     end
       
end


if gr
    
figure
ax(1)=subplot(211);
plot(RPM)
hold on
plot(maxloc_rpm,RPM(maxloc_rpm),'r*')
plot(minloc_rpm,RPM(minloc_rpm),'g*')

ax(2)=subplot(212)
plot(PCA)
hold on
plot(maxloc_pca,PCA(maxloc_pca),'r*')
plot(minloc_pca,PCA(minloc_pca),'g*')
linkaxes([ax(2) ax(1)],'x');
end


end

