clear
close all
clc
load ('TNTpatients.mat')
load('tropdata.mat')
k=1;
for i=2:length(vals)
    if vals(i,end)>0
%         k=k+1;
        idx(i)=find(Troponinexcel.VarName14==vals(i,end),1);

compositions(i,:) = [vals(i,end),Troponinexcel.pNEU(idx(i)), Troponinexcel.pLYM(idx(i)),...
    Troponinexcel.pMON(idx(i)),Troponinexcel.pEOS(idx(i)),vals(i,2),...
    vals(i,3),vals(i,4),vals(i,5)];
        errorneu(i,:) = (vals(i,2)-Troponinexcel.pNEU(idx(i)))*100/Troponinexcel.pNEU(idx(i));
        errorlym(i,:) = (vals(i,3)-Troponinexcel.pLYM(idx(i)))*100/Troponinexcel.pLYM(idx(i));
        errormon(i,:) = (vals(i,4)-Troponinexcel.pMON(idx(i)))*100/Troponinexcel.pMON(idx(i));
        erroreos(i,:) = (vals(i,5)-Troponinexcel.pEOS(idx(i)))*100/Troponinexcel.pEOS(idx(i));
%     else
%         idx(i)=0;
%         errorneu(i) = 0;
%         errorlym(i) = 0;
%         errormon(i) = 0;
%         erroreos(i) = 0;
    end
%     %   errorneu = (vals(1+i,2)-Troponinexcel.pNEU(idx(i)))/Troponinexcel.pNEU(idx(i));
%     %   errorlym = (vals(1+i,3)-Troponinexcel.pLYM(idx(i)))/Troponinexcel.pLYM(idx(i));
%     %   errormon = (vals(1+i,4)-Troponinexcel.pMON(idx(i)))/Troponinexcel.pMON(idx(i));
%     %   erroreos = (vals(1+i,5)-Troponinexcel.pEOS(idx(i)))/Troponinexcel.pEOS(idx(i));
end
errors = [errorneu, errorlym, errormon, erroreos];

%% Compare the error with Changi Excel
ChangiExcel = xlsread('/Users/anweshachaudhury/Desktop/Anwesha research/Troponindata/Changi_DataLog_6301-7400_4.xlsx');
% ChangiExcel = xlsread('/Users/anweshachaudhury/Desktop/Anwesha research/Troponindata/Changi_DataLog_5101-6300_3');
    cd('/Users/anweshachaudhury/Desktop/Anwesha research/Analyzers');
load('allresultsInVal.mat')
for i=1:length(SIndex)
    if ismember(SIndex(i),ChangiExcel(:,1))==1
        idx(i)=find(ChangiExcel(:,1)==SIndex(i));
        errorneu(i,1) = (vals(i+1,2)-ChangiExcel(idx(i),37))*100/ChangiExcel(idx(i),37);
        errorlym(i,1) = (vals(i+1,3)-ChangiExcel(idx(i),41))*100/ChangiExcel(idx(i),41);
        errormon(i,1) = (vals(i+1,4)-ChangiExcel(idx(i),44))*100/ChangiExcel(idx(i),44);
        erroreos(i,1) = (vals(i+1,5)-ChangiExcel(idx(i),47))*100/ChangiExcel(idx(i),47);
    else
        idx(i)=0;
        errorneu(i,1) = 0;
        errorlym(i,1) = 0;
        errormon(i,1) = 0;
        erroreos(i,1) = 0;
    end
end