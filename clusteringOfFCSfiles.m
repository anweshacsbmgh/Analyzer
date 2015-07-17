% clear
% close all
% clc
function [data,popDyn,popDynoverall] = clusteringOfFCSfiles (filename,data,popDyn,popDynoverall)
% close all
% Indices are  1-IAS, 2-ALL, 3-PSS, 4-DSS, 5-FL3, 6-time
addpath('/Users/anweshachaudhury/Desktop/Anwesha research/MATLABfuns/fcs_read')
% cd ('/Users/anweshachaudhury/Desktop/Anwesha research/FCSfiles/FCS_6992-7107_42973AZ.20121211.075949')
cd ('/Users/anweshachaudhury/Desktop/Anwesha research/FCSfiles/allFCS');
% filename='42318AZ96.20120622.190108.3945';
[FCSData, FCSHdrStd, FCSHdrAbt] = ReadFCS_abbott(filename);
AbtHdrTbl = struct2table(FCSHdrAbt);
fcsdat1 = FCSData(1).RawData;
fcsdat1(:,3)=[]; %remove ALL_W
if length(fcsdat1)<2000
    disp('Too few cells')
    data(end+1,:)=zeros(1,7);
    return
end

%% Following the protocol provided by Don for WBC classification
%Initial analyze and take only cells under FL3 and ALL threshold for
%further analysis
k=histogram(fcsdat1(:,5),256);
cutFL3 = k.BinEdges(100);
k1 = histogram(fcsdat1(:,2),256);
cutALL = k1.BinEdges(40);
initIdx = find(fcsdat1(:,5)<=cutFL3 | fcsdat1(:,2)<=cutALL);
fcsdat = fcsdat1(initIdx,:);
close
% Classifying mono-nuclear and poly-nuclear cells
[T,C] = kmeans(fcsdat(:,[1,3]),2); % finding mono and poly clusters in PSS vs IAS plot
% plot(fcsdat(T==1,1),fcsdat(T==1,3),'r*')
% hold on
% plot(fcsdat(T==2,1),fcsdat(T==2,3),'b*')
% xlabel('IAS')
% ylabel('PSS')
fidx(1,1)=max(fcsdat(T==1,3));
fidx(2,2)=min(fcsdat(T==2,3));
fidx(1,2)=max(fcsdat(T==2,3));
fidx(2,1)=min(fcsdat(T==1,3));
polyCluster = find(fidx(1,:)==max(max(fidx))); %identify the cluster number for the "polynucleated cells"
% poly=fcsdat(ismember(T,polyCluster),:); %polynucleated cells
% mono=fcsdat(~ismember(T,polyCluster),:); %mononucleated cells
k1=histogram(fcsdat(:,3),64);
dummyidx = find(k1.Values(10:25)==min(k1.Values(10:25)),1);
poly = fcsdat(find(fcsdat(:,3)>=k1.BinEdges(9+dummyidx)),:);
mono = fcsdat(find(fcsdat(:,3)<k1.BinEdges(9+dummyidx)),:);
close
figure
plot(poly(:,1),poly(:,3),'r.',mono(:,1),mono(:,3),'b.')
figure
h = histogram((mono(:,5)),128); %128 instead of 256
xlabel('FL3')
cutt = find(h.Values(40:70) == min(h.Values(40:70)),1);
cutOffFL3 = h.BinEdges(39+cutt); %defining the cutoff using FL3 histogram at 110 channel (used 55 instead of 110)
highFL3=mono(mono(:,5)>=cutOffFL3,:);
lowFL3=mono(mono(:,5)<cutOffFL3,:);

% manipulation using ALL (for NRBC)
% if length(highFL3)<20
%     n = 7;
%     dis=1;
% else
if length(highFL3)>50
if length(highFL3)<100
    n=50;
    dis=4;
else
    n=180;
    dis = 20;
end
figure
h1 = histogram(highFL3(:,2),n); %180 channels instead of 256
xlabel('ALL')
[pks,idx,w,p] = findpeaks(h1.Values,'MinPeakProminence',max(h1.Values)-max(h1.Values)/4,'MinPeakDistance',dis);
if isempty(pks)==1 || length(pks)<2
    disp('Two peaks not identified in ALL plot for highFL3 population, please modify the findpeaks parameter')
    % keyboard
    [pks,idx,w,p] = findpeaks(h1.Values,'MinPeakDistance',length(h1.Values)/4);
end
if length(pks)<2
    idx(end+1)=length(h1.Values)-ceil(n/3)+find(h1.Values(end-floor(n/3):end)==max(h1.Values(end-floor(n/3):end)));
end
isodd=mod(length(idx),2); % identifying number of peaks identified (even or odd)
centralPeak = h1.BinEdges(idx((length(idx)+1)/2*isodd+(length(idx))/2*~isodd)); % setting the middle peak as the cutoff
middle = idx((length(idx)+1)/2*isodd+(length(idx))/2*~isodd);
low = idx(max(1,(length(idx)+1)/2*isodd+(length(idx))/2*~isodd-1)); % finding index of peak previous to middle peak
upp = idx((length(idx)+1)/2*isodd+(length(idx))/2*~isodd+1);  % finding index of peak after middle peak

% tentatively labeling the NRBC population
lowBounddum=ceil(length(find(h1.Values(low:middle)==...
    min(h1.Values(low:middle))))/2);
lowBound=low+max(find(h1.Values(low:middle)==min(h1.Values(low:middle)),lowBounddum));
uppBounddum=ceil(length(find(h1.Values(middle:upp-1)==...
    min(h1.Values(middle:upp-1))))/2);
uppBound=middle+max(find(h1.Values(middle:upp-1)==min(h1.Values(middle:upp-1)),uppBounddum));
NRBC=highFL3(highFL3(:,2)>=h1.BinEdges(lowBound) & highFL3(:,2)<=h1.BinEdges(uppBound),:);
dead = highFL3(highFL3(:,2)<h1.BinEdges(lowBound),:);
mononuclear=highFL3(highFL3(:,2)>h1.BinEdges(uppBound),:);
remainingMono = [lowFL3;mononuclear];
iii = find(remainingMono(:,5)==0);
% remainingMono(iii,:) = []; 
else
    remainingMono = [lowFL3;highFL3];
    iii = find(remainingMono(:,5)==0);
% remainingMono(iii,:) = []; 
    mononuclear=zeros(1,6);
    NRBC = zeros(1,6);
    dead = zeros(1,6);
end
% further analysis using the poly nuclear cells/unknown cells
arctanVector = atand((poly(:,4)-5)./(poly(:,3)-25)); % tan inv of DSS/PSS
figure
h2 = histogram(arctanVector,50); %histogram of the arctan plot into 50 bins
xlabel('arctan(DSS/PSS)')

lowPop = h2.BinEdges(1:end-1)<22; % 22 is the default angle
highPop = h2.BinEdges(1:end-1)>22;
if length(h2.Values)-1<=length(find(lowPop))
    findPeak1 = max(h2.Values(find(lowPop(1:end-1))));
    findPeak2=findPeak1/1000;
    resp=1;
    disp('One peak')

else
    findPeak1 = max(h2.Values(find(lowPop)));
    findPeak2 = max(h2.Values(find(highPop(1:end-1))));
     resp=2;
    disp('Two peaks')
end

% if length(h2.Values)>length(find(lowPop))-1
%     resp=2;
%     disp('Two peaks')
% else
% resp=1;
%     disp('One peak')
% end
idxPeak1=find(h2.Values(find(lowPop))==findPeak1,1,'last');
if resp==2
idxPeak2=find(highPop,1)+find(h2.Values(find(highPop(1:end-1)))==findPeak2);
else
    idxPeak2=length(h2.Values);
end
if resp == 1
    peak=h2.BinEdges(find(h2.Values==max(findPeak1,findPeak2)));
    LineAngle = 22;
else
    peak = [h2.BinEdges(find(h2.Values==findPeak1)),h2.BinEdges(find(h2.Values==findPeak2,1,'last'))];
    valley = idxPeak1-1+h2.BinEdges(find(h2.Values(idxPeak1:idxPeak2)==min(h2.Values(idxPeak1:idxPeak2)),1));
    LineAngle = valley;
end

DSS0 = 5;
PSS0 = 25;
slope = tand(LineAngle);
DSS1 = 12; %one segment
PSS1= (DSS1-DSS0)/slope+PSS0;


DSS2=36;
PSS2= (DSS2-DSS0)/slope+PSS0;
slopeLine2 = slope*0.8;
yInterceptLine2 = DSS2 - slopeLine2*PSS2;
indices = (poly(:,4) > slopeLine2*poly(:,3)+yInterceptLine2) & poly(:,3)>PSS2;
xx=linspace(1,max(poly(:,3)),ceil(max(poly(:,3))/10));
y1 = 12.*(xx<=PSS1)+ ((PSS1-PSS0)*slope+DSS0).*(xx>PSS1 & xx<=PSS2)+(slopeLine2.*xx+yInterceptLine2).*(xx>PSS2);
figure
plot(poly(find(indices),3),poly(find(indices),4),'r.',poly(find(~indices),3),poly(find(~indices),4),'b.')
hold on
plot(xx,y1)
eosino = poly(find(indices),:);
neutro = poly(find(~indices),:);
% close all
% Classifying monocytes vs lymphocytes
figure
h3 = histogram((remainingMono(:,2)),70);
% figure
% hfit = histfit((remainingMono(:,2)),70,'kernel');
% xlabel('ALL')
% [pks1,idx11,w,p] = findpeaks(-hfit(2).YData,'MinPeakDistance',10);

[t1,c1]=kmeans(remainingMono(:,2),3,'Start',[0.4e4;1.3e4;2.1e4]);
% c2=sort(c1);
monocyteidx = find(t1==3);
stromaidx = find(t1==1);
stroma = remainingMono(stromaidx,:);
monocytes = remainingMono(monocyteidx,:);
remainingMononew = remainingMono(find(~ismember(1:length(remainingMono),[stromaidx;monocyteidx])),:);


if length(stroma)*100/length(fcsdat)>=20 %if stoma is more than 20%
    figure
    hfit = histfit((remainingMono(:,2)),100,'kernel');
    xlabel('ALL')
    [pks1,idx11,w,p] = findpeaks(hfit(2).YData,'MinPeakDistance',10);
        dummyI = find(hfit(2).YData == max(hfit(2).YData));

    if length(idx11)==2 
    lowbound1 = 2;
    uppbound1 = idx11(1)+find(hfit(2).YData(idx11(1):idx11(2))==min(hfit(2).YData(idx11(1):idx11(2))),1);

    else if dummyI==idx11(1)
         lowbound1 = idx11(1);  
         uppbound1 = idx11(2)+find(hfit(2).YData(idx11(1):idx11(2))==min(hfit(2).YData(idx11(1):idx11(2))));
        else
        dummyI = find(hfit(2).YData == max(hfit(2).YData));
        for i = 1:length(idx11)-1
            if idx11(i)<=dummyI && idx11(i+1)>dummyI
        uppbound1 = idx11(i)+find(hfit(2).YData(idx11(i):idx11(i+1))==min(hfit(2).YData(idx11(i):idx11(i+1))));
        lowbound1 =  idx11(i-1)+find(hfit(2).YData(idx11(i-1):idx11(i))==min(hfit(2).YData(idx11(i-1):idx11(i))));
            end
        end
        end
%             uppbound1 = idx11(2)+find(hfit(2).YData(idx11(2):idx11(3))==min(hfit(2).YData(idx11(2):idx11(3))),1);
%             lowbound1 =  idx11(1)+find(hfit(2).YData(idx11(1):idx11(2))==min(hfit(2).YData(idx11(1):idx11(2))),1);
end
    monocyteidx = find(remainingMono(:,2)>hfit(2).XData(uppbound1));
    stromaidx = find(remainingMono(:,2)<hfit(2).XData(lowbound1));
    monocytes = remainingMono(monocyteidx,:);
    stroma = remainingMono(stromaidx,:);
    remainingMononew = remainingMono(find(~ismember(1:length(remainingMono),[stromaidx;monocyteidx])),:);
end

% figure
% hist(monocytes(:,2))
% figure
% hist(remainingMononew(:,2))
% figure
% hist(stroma(:,2))
% if length(idx11)==2
% uppbound1=idx11(2);
% lowbound1=idx11(1);
% else if length(idx11)>2
%         dummyI = find(hfit(2).YData == max(hfit(2).YData));
%         for i = 1:length(idx11)-1
%             if idx11(i)<dummyI && idx11(i+1)>dummyI
%         uppbound1 = idx11(i+1);
%         lowbound1 =  idx11(i);
%             end
%         end
%     end
% end




% % idx1 = peakfinder(h3.Values,(max(h3.Values)-min(h3.Values))/4); %Finding the peaks of the histogram
% prom=3.5;
% [pks1,idx11,w,p] = findpeaks(hfit(2).YData,'MinPeakProminence',prom,'MinPeakDistance',10);
% [pks1,idx12,w,p] = findpeaks(h3.Values,'MinPeakProminence',prom,'MinPeakDistance',8);
% if length(idx11)>length(idx12)
%     idx1 = idx11;
% else
%     idx1 = idx12;
% end
% while length(idx1)<2
%     prom = prom-0.5;
%     [pks1,idx1,w,p] = findpeaks(hfit(2).YData,'MinPeakProminence',prom,'MinPeakDistance',10);
%     dummyindx = ones(1,length(idx1))*find(hfit(2).YData==max(hfit(2).YData));
% end
% 
% % if mod(length(idx1),2)>0
% %     mid1 = median(idx1);
% % else
% %     mid1=idx1(length(idx1)/2);
% % end
% if length(idx1)>=2
%     if length(idx1)==length(idx11)
%         mid1 = find(hfit(2).YData==max(hfit(2).YData),1);
%     else if length(idx1)==length(idx12)
%             mid1 = find(h3.Values==max(h3.Values),1);
%         else
%             mid1 = find(hfit(2).YData==max(hfit(2).YData),1);
%             
%         end
%     end
% end
% if mid1>idx1(1) & mid1<idx1(end)
%     low1 = idx1(find(idx1==mid1)-1);
%     upp1 = idx1(find(idx1==mid1)+1);
% else if mid1==idx1(1)
%         low1 = idx1(1);
%         upp1 = idx1(2);
%     else if mid1==idx1(end)
%             low1 = idx1(end-1);
%             upp1 = idx1(end);
%         end
%     end
% end
% 
% if length(idx1)>2
%     if length(idx1)==length(idx11)
%         if mid1>idx1(1) & mid1<idx1(end)
%             lowbound1 = low1+find(hfit(2).YData(low1:mid1)==min(hfit(2).YData(low1:mid1)),1);
%             uppbound1 = mid1+find(hfit(2).YData(mid1:upp1)==min(hfit(2).YData(mid1:upp1)),1);
%         else if    mid1==idx1(1)
%                 lowbound1 = 2+find(hfit(2).YData(2:idx1(1))==min(hfit(2).YData(2:idx1(1))),1);
%                 uppbound1 = idx1(1)+find(hfit(2).YData(idx1(1):idx1(2))==min(hfit(2).YData(idx1(1):idx1(2))),1);
%             else if mid1==idx1(end)
%                     uppbound1 = idx1(2)+find(hfit(2).YData(idx1(2):min(length(hfit(2).YData)-idx1(1),idx1(2)+30))==min(hfit(2).YData(idx1(2):min(length(hfit(2).YData)-idx1(1),idx1(2)+30))),1);
%                     lowbound1 = idx1(1)+find(hfit(2).YData(idx1(1):idx1(2))==min(hfit(2).YData(idx1(1):idx1(2))),1);
%                     
%                 end
%             end
%         end
%     else if length(idx1)==length(idx12)
%             if mid1>idx1(1) & mid1<idx1(end)
%                 lowbound1 = low1+find(h3.Values(low1:mid1)==min(h3.Values(low1:mid1)),1);
%                 uppbound1 = mid1+find(h3.Values(mid1:upp1)==min(h3.Values(mid1:upp1)),1);
%             else if mid1==idx1(1)
%                     lowbound1 = 2+find(h3.Values(2:idx1(1))==min(h3.Values(2:idx1(1))),1);
%                     uppbound1 = idx1(1)+find(h3.Values(idx1(1):idx1(2))==min(h3.Values(idx1(1):idx1(2))),1);
%                 else if mid1==idx1(end)
%                         uppbound1 = idx1(2)+find(h3.Values(idx1(2):min(length(h3.Values)-idx1(1),idx1(2)+30))==min(h3.Values(idx1(2):min(length(h3.Values)-idx1(1),idx1(2)+30))),1);
%                         lowbound1 = idx1(1)+find(h3.Values(idx1(1):idx1(2))==min(h3.Values(idx1(1):idx1(2))),1);
%                         
%                     end
%                 end
%             end
%         else if  mid1>idx1(1) & mid1<idx1(end)
%                 lowbound1 = low1+find(hfit(2).YData(low1:mid1)==min(hfit(2).YData(low1:mid1)),1);
%                 uppbound1 = mid1+find(hfit(2).YData(mid1:upp1)==min(hfit(2).YData(mid1:upp1)),1);
%                 
%                 
%             end
%         end
%     end
% end
% 
% if length(idx1)==2 & ismember(idx1,idx11)==[1 1]
%     dummyindx = find(hfit(2).YData==max(hfit(2).YData)).*(ismember(idx1,idx11))+find(h3.Values==max(h3.Values)).*(~ismember(idx1,idx11));
%     if dummyindx(1)==idx1(1)
%         lowbound1 = 2+find(hfit(2).YData(2:idx1(1))==min(hfit(2).YData(2:idx1(1))),1);
%         uppbound1 = idx1(1)+find(hfit(2).YData(idx1(1):idx1(2))==min(hfit(2).YData(idx1(1):idx1(2))),1);
%     else
%         uppbound1 = idx1(2)+find(hfit(2).YData(idx1(2):min(length(hfit(2).YData)-idx1(1),idx1(2)+30))==min(hfit(2).YData(idx1(2):min(length(hfit(2).YData)-idx1(1),idx1(2)+30))),1);
%         lowbound1 = idx1(1)+find(hfit(2).YData(idx1(1):idx1(2))==min(hfit(2).YData(idx1(1):idx1(2))),1);
%     end
% else if length(idx1)==2 & ismember(idx1,idx12)==[1 1]
%         dummyindx = find(h3.Values==max(h3.Values)).*(~ismember(idx1,idx11));
%         
%         if dummyindx(1)==idx1(1)
%             lowbound1 = 2+find(h3.Values(2:idx1(1))==min(h3.Values(2:idx1(1))),1);
%             uppbound1 = idx1(1)+find(h3.Values(idx1(1):idx1(2))==min(h3.Values(idx1(1):idx1(2))),1);
%         else
%             uppbound1 = idx1(2)+find(h3.Values(idx1(2):min(length(h3.Values),idx1(2)+30))==min(h3.Values(idx1(2):min(length(h3.Values),idx1(2)+30))),1);
%             lowbound1 = idx1(1)+find(h3.Values(idx1(1):idx1(2))==min(h3.Values(idx1(1):idx1(2))),1);
%         end
%     
%     
% else if length(idx1)>2 & max(idx1~=idx12)==1 & max(idx1~=idx11)==1
%         
%         if dummyindx(1)==idx1(1)
%             lowbound1 = 2+find(hfit(2).YData(2:idx1(1))==min(hfit(2).YData(2:idx1(1))),1);
%             uppbound1 = idx1(1)+find(hfit(2).YData(idx1(1):idx1(2))==min(hfit(2).YData(idx1(1):idx1(2))),1);
%         else
%             uppbound1 = idx1(2)+find(hfit(2).YData(idx1(2):min(length(hfit(2).YData),idx1(2)+30))==min(hfit(2).YData(idx1(2):min(length(hfit(2).YData),idx1(2)+30))),1);
%             lowbound1 = idx1(1)+find(hfit(2).YData(idx1(1):idx1(2))==min(hfit(2).YData(idx1(1):idx1(2))),1);
%         end
%     
% else if length(idx1)==2 & max(idx1~=idx12)==1 & max(idx1~=idx11)==1
%         
%         if dummyindx(1)==idx1(1)
%             lowbound1 = 2+find(hfit(2).YData(2:idx1(1))==min(hfit(2).YData(2:idx1(1))),1);
%             uppbound1 = idx1(1)+find(hfit(2).YData(idx1(1):idx1(2))==min(hfit(2).YData(idx1(1):idx1(2))),1);
%         else
%             uppbound1 = idx1(2)+find(hfit(2).YData(idx1(2):min(length(hfit(2).YData),idx1(2)+30))==...
%                 min(hfit(2).YData(idx1(2):min(length(hfit(2).YData),idx1(2)+30))),1);
%             lowbound1 = idx1(1)+find(hfit(2).YData(idx1(1):idx1(2))==min(hfit(2).YData(idx1(1):idx1(2))),1);
%         end
%     end
%     end
%     
%         
%     end
% end


% monocyteidx = find(remainingMono(:,2)>h3.BinEdges(uppbound1));
% stromaidx = find(remainingMono(:,2)<h3.BinEdges(lowbound1));
% monocytes = remainingMono(monocyteidx,:);
% stroma = remainingMono(stromaidx,:);
% remainingMononew = remainingMono(find(~ismember(1:length(remainingMono),[stromaidx;monocyteidx])),:);

%% 
figure
plot(remainingMononew(:,1),remainingMononew(:,2),'.')
LineFit=lsline;
slopeOfB = (LineFit.YData(2)-LineFit.YData(1))/(LineFit.XData(2)-LineFit.XData(1));
slopeOfA = -1/slopeOfB;
addpath('/Users/anweshachaudhury/Desktop/Anwesha research/MATLABfuns/geom2d-2014.10.27/geom2d/geom2d')
indx=find(remainingMononew(:,1)==max(remainingMononew(:,1)));
lineB = createLine(LineFit.XData(2),LineFit.YData(2),1,slopeOfB);
lineA=createLine(remainingMononew(indx,1),remainingMononew(indx,2),1,slopeOfA);
for i=1:length(remainingMononew(:,1))
    projectionsA(i,:) =  projPointOnLine([remainingMononew(i,1),remainingMononew(i,2)], lineA);
end

for i=1:length(remainingMononew(:,1))
    projectionsB(i,:) =  projPointOnLine([remainingMononew(i,1),remainingMononew(i,2)], lineB);
end
figure
plot(neutro(:,1),neutro(:,2),'r.',monocytes(:,1),monocytes(:,2),'b.',eosino(:,1),eosino(:,2),'g.',remainingMononew(:,1),remainingMononew(:,2),'y.')
ylim([0 3.5e4]);
legend('Neutrophil','Monocyte','Eosinophil','Lymphocytes')
totWBC = (length(neutro)+length(projectionsB)+...
    length(monocytes)+length(eosino));
% keyboard
% Fitting Gaussian models to populations
for i = 1:5
[muN(i),sigmaN(i)]=normfit(neutro(:,i));
[muL(i),sigmaL(i)]=normfit(remainingMononew(:,i));
[muM(i),sigmaM(i)]=normfit(monocytes(:,i));
end
[muNoverall,sigmaNoverall] = normfit(neutro(:,1:5));
covN = cov(neutro(:,1:5));
[muLoverall,sigmaLoverall] = normfit(remainingMononew(:,1:5));
covL = cov(remainingMononew(:,1:5));
[muMoverall,sigmaMoverall] = normfit(monocytes(:,1:5));
covM = cov(monocytes(:,1:5));
% disp(['The patient filename is ',filename]);
% disp(['Total cells ',num2str(length(fcsdat))]);
disp(['Neutrophils% ',num2str(length(neutro)*100/(length(neutro)+length(projectionsB)+...
    length(monocytes)+length(eosino)))]);
disp(['Lymphocytes% ',num2str(length(projectionsB)*100/(length(neutro)+length(projectionsB)+...
    length(monocytes)+length(eosino)))]);
disp(['Monocytes% ',num2str(length(monocytes)*100/(length(neutro)+length(projectionsB)+...
    length(monocytes)+length(eosino)))]);
disp(['Eosinophils% ',num2str(length(eosino)*100/(length(neutro)+length(projectionsB)+...
    length(monocytes)+length(eosino)))]);
% disp(['NRBC% ',num2str(length(NRBC)*100/(length(neutro)+length(projectionsB)+...
%     length(monocytes)+length(eosino)))]);
% keyboard
popDyn(end+1,:) = [muN,sigmaN,muL,sigmaL,muM,sigmaM];
popDynoverall(end+1,:) = [muNoverall,sigmaNoverall,muLoverall,...
    sigmaLoverall,muMoverall,sigmaMoverall,covN(:)',covL(:)',covM(:)'];
data(end+1,:)=[length(fcsdat),length(neutro)*100/(totWBC),length(remainingMononew)*100/totWBC,...
    length(monocytes)*100/totWBC,length(eosino)*100/totWBC,length(NRBC)*100/totWBC,str2num(filename(end-3:end))];
cd ('/Users/anweshachaudhury/Desktop/Anwesha research/FCSfiles/WBCresults')

save(filename(end-3:end))