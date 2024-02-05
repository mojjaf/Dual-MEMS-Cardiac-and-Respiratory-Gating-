function [ listmode_data_5bin ] = MinMaxFinder( ADRx,ADRy,GDRx,GDRy,KalmanFusion,MedianFusion, PCAFusion )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

maxloc_adrX = ampd(ADRx);
minloc_adrX = ampd(-ADRx);


maxloc_adrY = ampd(ADRy);
minloc_adrY= ampd(-ADRy);

maxloc_gdrX = ampd(GDRx);
minloc_gdrX = ampd(-GDRx);

maxloc_gdrY = ampd(GDRy);
minloc_gdrY = ampd(-GDRy);

maxloc_klm = ampd(KalmanFusion);
minloc_klm = ampd(-KalmanFusion);

maxloc_med = ampd(MedianFusion);
minloc_med = ampd(-MedianFusion);

maxloc_pca= ampd(PCAFusion);
minloc_pac= ampd(-PCAFusion);

figure
plot(ADRx)
hold on
plot(maxloc_adrX,ADRx(maxloc_adrX),'r*')
plot(minloc_adrX,ADRx(minloc_adrX),'g*')

plot(ADRy)
plot(maxloc_adrY,ADRy(maxloc_adrY),'r+')
plot(minloc_adrY,ADRy(minloc_adrY),'g+')

plot(GDRx)
plot(maxloc_gdrX,GDRx(maxloc_gdrX),'k+')
plot(minloc_gdrX,GDRx(minloc_gdrX),'g+')

plot(GDRy)
plot(maxloc_gdrY,GDRy(maxloc_gdrY),'bo')
plot(minloc_gdrY,GDRy(minloc_gdrY),'co')

plot(MedianFusion)
plot(maxloc_med,MedianFusion(maxloc_med),'bo')
plot(minloc_med,MedianFusion(minloc_med),'co')

% plot(PCAFusion)
% plot(maxloc_pca,PCAFusion(maxloc_pca),'k*')
% plot(minloc_pac,PCAFusion(minloc_pac),'c*')


% all_maxlocs={maxloc_adrX,maxloc_adrY,maxloc_gdrX,maxloc_gdrY,maxloc_klm,maxloc_med,maxloc_pca};
% all_minlocs={minloc_adrX,minloc_adrY,minloc_gdrX,minloc_gdrY,minloc_klm,minloc_med,minloc_pac};
% 
all_maxlocs={maxloc_adrX,maxloc_adrY,maxloc_gdrX,maxloc_gdrY,maxloc_med};
all_minlocs={minloc_adrX,minloc_adrY,minloc_gdrX,minloc_gdrY,minloc_med};
%%
all_maxlocs_vec=cell2mat(all_maxlocs);
temp=all_maxlocs_vec;
max_indx=[];
j=1;
i=1;
 while (numel(temp)>2)
                  distances= abs(all_maxlocs_vec(i)-temp);
                  [flag,indx]=find(distances<100);
                  if numel(flag>2)
                   max_indx(j)=median(temp(indx));
                   temp(indx)=[];
                   j=j+1;
                  end
                  i=i+1;
          
 end
 
 
 maxlocs=int32(max_indx);
 
all_minlocs_vec=cell2mat(all_minlocs);
temp_min=all_minlocs_vec;
min_indx=[];
jj=1;
ii=1;
 while (numel(temp_min)>2)
                  distances= abs(all_minlocs_vec(ii)-temp_min);
                  [flag,indx]=find(distances<75);
                  if numel(flag>2)
                   min_indx(jj)=median(temp_min(indx));
                   temp_min(indx)=[];
                   jj=jj+1;
                  end
                  ii=ii+1;
          
 end
 
        
       minlocs=int32(min_indx);


%%
figure
plot(ADRx)
hold on
plot(maxlocs,ADRx(maxlocs),'r*')
plot(minlocs,ADRx(minlocs),'g*')

plot(ADRy)
plot(maxlocs,ADRy(maxlocs),'r+')
plot(minlocs,ADRy(minlocs),'g+')

plot(GDRx)
plot(maxlocs,GDRx(maxlocs),'k+')
plot(minlocs,GDRx(minlocs),'g+')

plot(GDRy)
plot(maxlocs,GDRy(maxlocs),'bo')
plot(minlocs,GDRy(minlocs),'co')

plot(MedianFusion,'LineWidth',2,'Linestyle','-','color','b')
plot(maxlocs,MedianFusion(maxlocs),'bo')
plot(minlocs,MedianFusion(minlocs),'co')
line(repmat(minlocs,[2 1]),repmat([min(MedianFusion-mean(MedianFusion)); max(MedianFusion-mean(MedianFusion))],size(minlocs)),'LineWidth',1,'LineStyle','-.','Color','r');
%%
minlocs=sort(minlocs);
fs=25;
  Resp_cycle_length=[];
MEMS_resp_bin=[];
bin_num=[];
format short
for i=1:numel(minlocs)-1
Resp_cycle_length(i)=minlocs(i+1)-minlocs(i);
MEMS_resp_bin(i,1)=minlocs(i);
bin_num(i,1)=1;
for j=2:5
MEMS_resp_bin(i,j)=MEMS_resp_bin(i,j-1)+Resp_cycle_length(i)/5;
bin_num(i,j)=j;
end
gating_bins = reshape( MEMS_resp_bin.' ,1,numel(MEMS_resp_bin));
bin_nums = reshape( bin_num.' ,1,numel(bin_num));
listmode_data=[round(1000*gating_bins')/fs bin_nums'];
end
listmode_data_5bin=listmode_data;


figure
plot(MedianFusion-mean(MedianFusion));title('End-Expiratory Cycle Segmentation');axis tight;
line(repmat(gating_bins,[2 1]),repmat([min(MedianFusion-mean(MedianFusion))/2; max(MedianFusion-mean(MedianFusion))/2],size(gating_bins)),'LineWidth',2.5,'LineStyle','-.','Color','g');
line(repmat(minlocs,[2 1]),repmat([min(MedianFusion-mean(MedianFusion))/2; max(MedianFusion-mean(MedianFusion))/2],size(minlocs)),'LineWidth',2.5,'LineStyle','-.','Color','r');


end
