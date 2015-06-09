% clear
% close all
% clc
function [data] = clusteringOfFCSfiles (filename,data)
close all
% Indices are  1-IAS, 2-ALL, 3-PSS, 4-DSS, 5-FL3, 6-time
addpath('/Users/anweshachaudhury/Desktop/Anwesha research/MATLABfuns/fcs_read')
cd ('/Users/anweshachaudhury/Desktop/Anwesha research/FCSfiles/FCS_6200-6249_42973AZ.20121107.121214')
% filename='42318AZ96.20120622.190108.3945';
[fcsdat, fcshdr, fcsdatscaled, fcsdatcomp]=fca_readfcs(filename);
% [FCSData, FCSHdrStd, FCSHdrAbt] = ReadFCS_abbott(filename);
fcsdat(:,3)=[]; %remove ALL_W
if length(fcsdat)<500
    disp('Too few cells')
    return
end

%% Following the protocol provided by Don for WBC classification

% Classifying mono-nuclear and poly-nuclear cells
[T,C] = kmeans(fcsdat(:,[1,3]),2); % finding mono and poly clusters in PSS vs IAS plot
plot(fcsdat(T==1,1),fcsdat(T==1,3),'r*')
hold on
plot(fcsdat(T==2,1),fcsdat(T==2,3),'b*')
xlabel('IAS')
ylabel('PSS')
fidx(1,1)=max(fcsdat(T==1,3));
fidx(2,2)=min(fcsdat(T==2,3));
fidx(1,2)=max(fcsdat(T==2,3));
fidx(2,1)=min(fcsdat(T==1,3));
polyCluster = find(fidx(1,:)==max(max(fidx))); %identify the cluster number for the "polynucleated cells"
poly=fcsdat(ismember(T,polyCluster),:); %polynucleated cells
mono=fcsdat(~ismember(T,polyCluster),:); %mononucleated cells
figure
h = histogram((mono(:,5)),128); %128 instead of 256
xlabel('FL3')
cutOffFL3 = h.BinEdges(55); %defining the cutoff using FL3 histogram at 110 channel (used 55 instead of 110)
highFL3=mono(mono(:,5)>=cutOffFL3,:);
lowFL3=mono(mono(:,5)<cutOffFL3,:);

% manipulation using ALL (for NRBC)
figure
h1 = histogram((highFL3(:,2)),180); %128 channels instead of 256
xlabel('ALL')
% addpath('/Users/anweshachaudhury/Desktop/Anwesha research/MATLABfuns/peakfinder')
% idx = peakfinder(h1.Values,(max(h1.Values)-min(h1.Values))/4,max(h1.Values)-1); %Finding the peaks of the histogram
[pks,idx,w,p] = findpeaks(h1.Values,'MinPeakProminence',max(h1.Values)-max(h1.Values)/4,'MinPeakDistance',20);
if isempty(pks)==1 || length(pks)<2
    disp('Two peaks not identified in ALL plot for highFL3 population, please modify the findpeaks parameter')
% keyboard
    [pks,idx,w,p] = findpeaks(h1.Values,'MinPeakDistance',length(h1.Values)/4);
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
uppBounddum=ceil(length(find(h1.Values(middle:upp)==...
    min(h1.Values(middle:upp))))/2);
uppBound=middle+max(find(h1.Values(middle:upp)==min(h1.Values(middle:upp)),uppBounddum));
NRBC=highFL3(highFL3(:,2)>=h1.BinEdges(lowBound) & highFL3(:,2)<=h1.BinEdges(uppBound),:);
dead = highFL3(highFL3(:,2)<h1.BinEdges(lowBound),:);
mononuclear=highFL3(highFL3(:,2)>h1.BinEdges(uppBound),:);
remainingMono = [lowFL3;mononuclear];

% further analysis using the poly nuclear cells/unknown cells
arctanVector = atand((poly(:,4)-5)./(poly(:,3)-25)); % tan inv of DSS/PSS
figure
h2 = histogram(arctanVector,50); %histogram of the arctan plot into 50 bins
xlabel('arctan(DSS/PSS)')
lowPop = h2.BinEdges<22; % 22 is the default angle
highPop = h2.BinEdges>22;
findPeak1 = max(h2.Values(find(lowPop)));
findPeak2 = max(h2.Values(find(highPop(1:end-1))));
if max(findPeak1,findPeak2)/min(findPeak1,findPeak2)<4 %ratio between 2 peaks should not be more than 1/4
    resp = 2;
    disp('Two peaks')
else
    resp = 1;
    disp('One peak')
end
if resp == 1
    peak=h2.BinEdges(find(h2.Values==max(findPeak1,findPeak2)));
    LineAngle = peak;
else
    peak = [h2.BinEdges(find(h2.Values==findPeak1)),h2.BinEdges(find(h2.Values==findPeak2,1,'last'))];
    LineAngle = 22;
end

DSS0 = 5;
PSS0 = 25;
slope = tand(22);
DSS1 = 12; %one segment
PSS1= (DSS1-DSS0)/slope+PSS0;

DSS2=36;
PSS2= (DSS2-DSS0)/slope+PSS0;
slopeLine2 = 0.8*slope;
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

% Classifying monocytes vs lymphocytes
figure
h3 = histogram((remainingMono(:,2)),128);
xlabel('ALL')
% idx1 = peakfinder(h3.Values,(max(h3.Values)-min(h3.Values))/4); %Finding the peaks of the histogram
[pks1,idx1,w,p] = findpeaks(h3.Values,'MinPeakProminence',8,'MinPeakDistance',10);
if mod(length(idx1),2)>0
mid1 = median(idx1);
else
 mid1=idx1(length(idx1)/2+1);   
end
low1 = idx1(1);
upp1 = idx1(end);
lowbound1 = low1+find(h3.Values(low1:mid1)==min(h3.Values(low1:mid1)),1);
uppbound1 = mid1+find(h3.Values(mid1:upp1)==min(h3.Values(mid1:upp1)),1);
monocyteidx = find(remainingMono(:,2)>h3.BinEdges(uppbound1));
stromaidx = find(remainingMono(:,2)<h3.BinEdges(lowbound1));
monocytes = remainingMono(monocyteidx,:);
stroma = remainingMono(stromaidx,:);
remainingMononew = remainingMono(find(~ismember(1:length(remainingMono),[stromaidx;monocyteidx])),:);

figure
plot(remainingMononew(:,1),remainingMononew(:,2),'.')
LineFit=lsline;
slopeOfB = (LineFit.YData(2)-LineFit.YData(1))/(LineFit.XData(2)-LineFit.XData(1));
slopeOfA = -1/slopeOfB;
addpath('/Users/anweshachaudhury/Desktop/Anwesha research/MATLABfuns/geom2d-2014.10.27/geom2d/geom2d')
indx=find(remainingMononew(:,1)==max(remainingMononew(:,1)));
lineB = createLine(LineFit.XData(2),LineFit.YData(2),1,slopeOfB);
lineA=createLine(remainingMononew(indx,1),remainingMononew(indx,2),1,slopeOfA);
for i=1:length(remainingMononew)
    projectionsA(i,:) =  projPointOnLine([remainingMononew(i,1),remainingMononew(i,2)], lineA);
end

for i=1:length(remainingMononew)
    projectionsB(i,:) =  projPointOnLine([remainingMononew(i,1),remainingMononew(i,2)], lineB);
end

disp(['The patient filename is ',filename]);
disp(['Total WBC ',num2str(length(fcsdat))]);
disp(['Neutrophils% ',num2str(length(neutro)*100/length(fcsdat))]);
disp(['Lymphocytes% ',num2str(length(projectionsB)*100/length(fcsdat))]);
disp(['Monocytes% ',num2str(length(monocytes)*100/length(fcsdat))]);
disp(['Eosinophils% ',num2str(length(eosino)*100/length(fcsdat))]);
disp(['NRBC% ',num2str(length(NRBC)*100/length(fcsdat))]);
data(end+1,:)=[length(fcsdat),length(neutro)*100/length(fcsdat),length(projectionsB)*100/length(fcsdat),...
    length(monocytes)*100/length(fcsdat),length(eosino)*100/length(fcsdat),length(NRBC)*100/length(fcsdat)];