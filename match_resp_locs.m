function [ new_locs_rpm ] = match_resp_locs( locs_rpm,locs_mems )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

new_locs_rpm=[];
windowLength=5;
L=numel(locs_rpm);
step=1;
curPos=1;
numOfFrames = floor((L-windowLength)/step) + 1;

for i=1:numel(locs_mems)
 try  
       dist=locs_rpm(curPos:curPos+windowLength)-locs_mems(i);
       [a,b]=find(min(dist));
       temp=locs_rpm(curPos:curPos+windowLength);
       new_locs_rpm=[new_locs_rpm,temp(a)];
       curPos = curPos + step;
 catch
  end
end
end

