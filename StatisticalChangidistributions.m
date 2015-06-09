clear 
close all
clc
load('tropdata.mat')
distinfo = xlsread('ChangiCompareWBC.xlsx','CompareChangiDist');
distinfo(1,:)=[];
Safe_pat_xls = [12 13 18 22 23 29 31 33 37 40 41 54 55 56 59 60 61 68 71 75 77 79 82 88 89 93 95 100 101 104 108]'; %green markers on excel sheet

% a = [-2; find(~isnan(table2array(patientnum)))]';

for i=1:length(distinfo)
%     idx(i)=find(distinfo(i,1)==Troponinexcel.VarName14,1);
if ismember(distinfo(i,2),Safe_pat_xls)==1
    output(i)=1;
else
    output(i)=0;
end
end
% output=fliplr(output);
[bstep,se,pval,inmodel,statsstep,nextstep,history] = stepwisefit(distinfo(:,3:end),...
    output');

%% Statistical platelets
clear
load('tropdata.mat')
pltinfo=xlsread('ChangiCompareWBC.xlsx','PLTchangi');
Safe_pat_xls = [12 13 18 22 23 29 31 33 37 40 41 54 55 56 59 60 61 68 71 75 77 79 82 88 89 93 95 100 101 104 108]'; %green markers on excel sheet
pltinfo(1,:)=[];
for i=1:length(pltinfo)
%     idx(i)=find(distinfo(i,1)==Troponinexcel.VarName14,1);
if ismember(pltinfo(i,2),Safe_pat_xls)==1
    output(i)=1;
else
    output(i)=0;
end
end

[bstep,se,pval,inmodel,statsstep,nextstep,history] = stepwisefit(pltinfo(:,3:end),...
    output');