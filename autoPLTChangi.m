clear

close all
clc
load('/Users/anweshachaudhury/Desktop/somepaths.mat')
% cd ('/Users/anweshachaudhury/Desktop/Anwesha research/FCSfiles/allFCS');
% fnames = dir('/Users/anweshachaudhury/Desktop/Anwesha research/FCSfiles/allFCS');
fnames = cellstr(assignFile);
numfids = length(fnames);
vals = zeros(1,24);
MPV=0;
for K = 1:numfids
    if strcmp(assignFile(K,1),'0') || cellfun(@isempty,cellstr(assignFile(K,:)))
        vals(end+1,:) = zeros(1,24);
    else
    cd('/Users/anweshachaudhury/Desktop/Anwesha research/Analyzers');
    kname = strcat(pathFile{K},'/',fnames{K});
%     if ismember(str2num(fnames(K).name(end-3:end)),Troponinexcel.VarName14)==1
        vals = PltoIdentifier(kname,vals);
% id = find(Troponinexcel.VarName14==str2num(fnames(K).name(end-3:end)),1)
%         MPV(end+1) = Troponinexcel.MPV(id);
%         file(K,1)=str2num(fnames(K).name(end-3:end));
%     else
%         continue
%     end
    end
    close all
    K
%     SIndex(K) = str2num(fnames{K});
end
