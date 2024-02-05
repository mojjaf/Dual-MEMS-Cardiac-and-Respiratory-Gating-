function [ mat_MEMS ] = ensemble_average(signal,locs,fs,start,stop,gr  )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
n=numel(locs)-1;
mat_MEMS=zeros(n,(start+stop));
Jittar_on=start;
Jittar_off=stop;
for i=1:n
        peak_sample_n = locs(i); %in samples
        start = peak_sample_n-Jittar_on ;
        stop = locs(i)+Jittar_off;
        sample = signal(start:stop);
        
%       j=numel(sample);
          mat_MEMS(i,1:numel(sample))=sample;
          
    
end
     if gr
         figure
plot((mat_MEMS'))  
hold on
plot(mean((mat_MEMS)),'LineWidth',4.5,'LineStyle','-','Color','r')  

     end

end

