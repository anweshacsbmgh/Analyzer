% clear
% close all
% clc
% 
% addpath('/Users/anweshachaudhury/Desktop/Anwesha research/MATLABfuns/fcs_read')
% fnames = dir('/Users/anweshachaudhury/Desktop/Anwesha research/FCSfiles/allFCS');
% 
function data = PltoIdentifier (filename,data)
close all
% filename = fnames(12).name;
addpath('/Users/anweshachaudhury/Desktop/Anwesha research/MATLABfuns/fcs_read')
cd ('/Users/anweshachaudhury/Desktop/Anwesha research/FCSfiles/allFCS')

[FCSData, FCSHdrStd, FCSHdrAbt] = ReadFCS_abbott(filename);
AbtHdrTbl = struct2table(FCSHdrAbt);
fcsdat = FCSData(4).RawData;
figure
plot(fcsdat(:,1),fcsdat(:,2),'.')
LineReg = lsline;
slope = ((LineReg.YData(2)-LineReg.YData(1))/(LineReg.XData(2)-LineReg.XData(1)));
[b,bint,res] = regress(fcsdat(:,2),[ones(size(fcsdat(:,1))),fcsdat(:,1)]);
Bounds = quantile(res,[0.02,0.5,0.98]);
Bounds1 = [0.95*max(res),0,0.95*min(res)];

cIntercept = LineReg.YData(2)-slope*LineReg.XData(2);
lowBound = Bounds(1);
UppBound = Bounds(3);
% noPLTIdx = find(res>=UppBound | res<=lowBound);
% PLTIdx = find(~ismember(1:length(fcsdat),noPLTIdx));
bbound = min(abs(Bounds(3)),abs(Bounds(1)));
c_upp = cIntercept+bbound/cos(atan(slope))*sign(Bounds(3));
c_low = cIntercept+bbound/cos(atan(slope))*sign(Bounds(1));
% hold on
% plot(fcsdat(PLTIdx,1),fcsdat(PLTIdx,2),'b.',fcsdat(noPLTIdx,1),fcsdat(noPLTIdx,2),'r.')
% addpath('/Users/anweshachaudhury/Desktop/Anwesha research/MATLABfuns/geom2d-2014.10.27/geom2d/geom2d')
% lineLow=createLine(0,c_low,1,slope);
% lineUp=createLine(0,c_upp,1,slope);
hrefLow = refline(slope,c_low);
hrefHigh = refline(slope,c_upp);
yCalcUpp = slope.*fcsdat(:,1)+c_upp;
yCalcLow = slope.*fcsdat(:,1)+c_low;
PLTIdx = find(fcsdat(:,2)>=yCalcLow & fcsdat(:,2)<=yCalcUpp);
noPLTIdx = find(~ismember(1:length(fcsdat),PLTIdx));
hold on
plot(fcsdat(PLTIdx,1),fcsdat(PLTIdx,2),'b.',fcsdat(noPLTIdx,1),fcsdat(noPLTIdx,2),'r.')
% [h,stats] = cdfplot(fcsdat(PLTIdx,1:2));
muO = mean(fcsdat(PLTIdx,1:2));
sigmaO = std(fcsdat(PLTIdx,1:2));
% [muO,sigmaO]=normfit(fcsdat(PLTIdx,1:2));
coRelate = corrcoef(fcsdat(PLTIdx,1:2));
kurt = kurtosis(fcsdat(PLTIdx,1:2));
skew = skewness(fcsdat(PLTIdx,1:2));
modal = mode((fcsdat(PLTIdx,1:2)));
medn = median((fcsdat(PLTIdx,1:2)));
interQuart = iqr((fcsdat(PLTIdx,1:2)));
%% Impedence Platelet
fcsdat1 = FCSData(3).RawData;
LowRange = min(fcsdat1(:,1));
UppRange = max(fcsdat1(:,1));
hImp = histfit(fcsdat1(:,1),30,'lognormal');
gridTrans = (hImp(2).XData-hImp(2).XData(1))./(hImp(2).XData(end)-hImp(2).XData(1))*(70-2)+2;
% h = histfit(gridTrans,30,'lognormal');

% notPltsIdx = find(fcsdat1(:,1)<=grid(2) | fcsdat1(:,1)>=grid(20));
% fcsdat1(notPltsIdx,:) = [];

hLog = fitdist(fcsdat1(:,1),'lognormal');
muTransform = exp(hLog.mu);
MPV = (muTransform-hImp(2).XData(1))./(hImp(2).XData(end)-hImp(2).XData(1))*(70-2)+2;

if MPV <6 | MPV >15
    disp('MPV out of range')
end
% keyboard
data(end+1,:) = [str2num(filename(end-3:end)),hLog.mu,hLog.sigma,MPV,length(fcsdat1),length(fcsdat),muO(1),muO(2),sigmaO(1),sigmaO(2),kurt,skew,modal,medn,interQuart,coRelate(:)'];