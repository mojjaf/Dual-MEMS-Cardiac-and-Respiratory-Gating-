function [Precision,Sensitivity,Fscore,RMS_error] = statistical_evaluation( true_peak_indx, predicted_peak_indx,fs )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
truepos=[];
correct_peak=[];
for ii = 1 : numel(true_peak_indx)
     distances = abs(predicted_peak_indx - true_peak_indx(ii));
     mindist = min(distances);

     if mindist < 0.25*fs;
         truepos(ii) = 1;
        idx=find(distances==mindist)
         correct_peak(ii)=predicted_peak_indx(idx(1));
     else
         truepos(ii) = 0;
     end

end

truepos = mean(truepos);
falseneg=1-truepos;

for ii = 1 : numel(predicted_peak_indx)
     distances = abs(true_peak_indx - predicted_peak_indx(ii));
     mindist = min(distances);

     if mindist < 0.25*fs;
         truepos_b(ii) = 1;
                  

     else
         truepos_b(ii) = 0;
     end

end
truepos_b=mean(truepos_b);
falsepos=1-truepos_b;
if truepos==truepos_b
    
    Precision=(truepos/(truepos+falsepos));
    Sensitivity=(truepos/(truepos+falseneg));
    Fscore=2*truepos/(2*truepos+falsepos+falseneg);
    detection_E=(falseneg+falsepos)/numel(predicted_peak_indx);
    matched_peaks=correct_peak(correct_peak>0);
    RR_pred=diff(matched_peaks)/fs;
    RR_true=diff(true_peak_indx(correct_peak>0))/fs;
    RMS_error=rmse(RR_true,RR_pred)*1000;
else
    truepos=truepos_b;
    Precision=(truepos/(truepos+falsepos));
    Sensitivity=(truepos/(truepos+falseneg));
    Fscore=2*truepos/(2*truepos+falsepos+falseneg);
    detection_E=(falseneg+falsepos)/numel(predicted_peak_indx);
    matched_peaks=correct_peak(correct_peak>0);
    RR_pred=diff(matched_peaks)/fs;
    RR_true=diff(true_peak_indx(correct_peak>0))/fs;
    RMS_error=rmse(RR_true,RR_pred)*1000;
end


end









