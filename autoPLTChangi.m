clear

close all
clc
load('tropdata.mat')
% cd ('/Users/anweshachaudhury/Desktop/Anwesha research/FCSfiles/allFCS');
fnames = dir('/Users/anweshachaudhury/Desktop/Anwesha research/FCSfiles/allFCS');
numfids = length(fnames);
vals = zeros(1,24);
MPV=0;
for K = 4:numfids
    cd('/Users/anweshachaudhury/Desktop/Anwesha research/Analyzers');
    if ismember(str2num(fnames(K).name(end-3:end)),Troponinexcel.VarName14)==1
        vals = PltoIdentifier(fnames(K).name,vals);
% id = find(Troponinexcel.VarName14==str2num(fnames(K).name(end-3:end)),1)
%         MPV(end+1) = Troponinexcel.MPV(id);
%         file(K,1)=str2num(fnames(K).name(end-3:end));
    else
        continue
    end
    close all
    K
    SIndex(K-3) = str2num(fnames(K).name(end-3:end));
end
