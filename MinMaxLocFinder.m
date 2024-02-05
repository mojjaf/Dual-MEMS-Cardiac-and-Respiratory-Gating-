function [ minlocs,maxlocs ] = MinMaxLocFinder( ADRx,ADRy,GDRx,GDRy, PCAFusion )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


maxloc_adrX = ampd(ADRx);
minloc_adrX = ampd(-ADRx);


maxloc_adrY = ampd(ADRy);
minloc_adrY= ampd(-ADRy);

maxloc_gdrX = ampd(GDRx);
minloc_gdrX = ampd(-GDRx);

maxloc_gdrY = ampd(GDRy);
minloc_gdrY = ampd(-GDRy);

% maxloc_klm = ampd(KalmanFusion);
% minloc_klm = ampd(-KalmanFusion);

maxloc_pca = ampd(PCAFusion);
minloc_pca = ampd(-PCAFusion);
% 
% maxloc_pca= ampd(PCAFusion);
% minloc_pac= ampd(-PCAFusion);

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

plot(PCAFusion)
plot(maxloc_pca,PCAFusion(maxloc_pca),'bo')
plot(minloc_pca,PCAFusion(minloc_pca),'co')

% plot(PCAFusion)
% plot(maxloc_pca,PCAFusion(maxloc_pca),'k*')
% plot(minloc_pac,PCAFusion(minloc_pac),'c*')


% all_maxlocs={maxloc_adrX,maxloc_adrY,maxloc_gdrX,maxloc_gdrY,maxloc_klm,maxloc_med,maxloc_pca};
% all_minlocs={minloc_adrX,minloc_adrY,minloc_gdrX,minloc_gdrY,minloc_klm,minloc_med,minloc_pac};
% 
all_maxlocs={maxloc_adrX,maxloc_adrY,maxloc_gdrX,maxloc_gdrY,maxloc_pca};
all_minlocs={minloc_adrX,minloc_adrY,minloc_gdrX,minloc_gdrY,minloc_pca};
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
plot(ADRx,'LineWidth',1.5,'Linestyle','-')
hold on
plot(maxlocs,ADRx(maxlocs),'r+')
plot(minlocs,ADRx(minlocs),'g+')

plot(ADRy,'LineWidth',1.5,'Linestyle','-')
plot(maxlocs,ADRy(maxlocs),'r+')
plot(minlocs,ADRy(minlocs),'g+')

plot(GDRx,'LineWidth',1.5,'Linestyle','-')
plot(maxlocs,GDRx(maxlocs),'r+')
plot(minlocs,GDRx(minlocs),'g+')

plot(GDRy,'LineWidth',1.5,'Linestyle','-')
plot(maxlocs,GDRy(maxlocs),'r+')
plot(minlocs,GDRy(minlocs),'r+')

plot(PCAFusion,'LineWidth',2,'Linestyle','-.','color','r')
plot(maxlocs,PCAFusion(maxlocs),'bo')
plot(minlocs,PCAFusion(minlocs),'co')
% line(repmat(minlocs,[2 1]),repmat([min(PCAFusion-mean(PCAFusion)); max(PCAFusion-mean(PCAFusion))],size(minlocs)),'LineWidth',1,'LineStyle','-.','Color','r');
end

