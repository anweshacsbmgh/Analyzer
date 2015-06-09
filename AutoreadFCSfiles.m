clear 
close all
clc
% cd ('/Users/anweshachaudhury/Desktop/Anwesha research/FCSFiles/FCS_6200-6249_42973AZ.20121107.121214');
% fnames = dir('/Users/anweshachaudhury/Desktop/Anwesha research/FCSFiles/FCS_6992-7107_42973AZ.20121211.075949');
load ('SIndex6992_7107.mat')
numfids = length(SIndex);
% % vals = cell(1,numfids);
vals = zeros(1,6);
for K = 3:numfids
    if finding(K-2)>0
    cd('/Users/anweshachaudhury/Desktop/Anwesha research/Analyzers');
vals = clusteringOfFCSfiles(fnames(K).name,vals);
    end
% SIndex(K-2) = str2num(fnames(K).name(end-3:end));
end
save('results6992_71071.mat')