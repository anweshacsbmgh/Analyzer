clear

close all
clc
load('/Users/anweshachaudhury/Desktop/Anwesha research/Data Files/Troponindata/tropdata.mat')
% cd ('/Users/anweshachaudhury/Desktop/Anwesha research/FCSfiles/allFCS');
fnames = dir('/Users/anweshachaudhury/Desktop/Anwesha research/FCSfiles/allFCS');
numfids = length(fnames);
vals = zeros(1,7);
popDyn = zeros(1,30);
popDynoverall = zeros(1,105);
for K = 3+1
    cd('/Users/anweshachaudhury/Desktop/Anwesha research/Analyzers');
    if ismember(str2num(fnames(K).name(end-3:end)),Troponinexcel.VarName14)==1
[vals,popDyn,popDynoverall] = clusteringOfFCSfiles(fnames(K).name,vals,popDyn,popDynoverall);
%             idx=find(Troponinexcel.VarName14==vals(K,1),1);
% patientNum(K-3) = Troponinexcel.VarName12(idx);
    else
        continue
    end
close all
K
SIndex(K-3) = str2num(fnames(K).name(end-3:end));
end

% save('results.mat')

%% Notes
% Could not open "42973AZ.20121112.174959.6659, 42973AZ.20121116.102401.6824...
...42973AZ.20121119.155914.6980, 42973AZ.20121120.100138.7018" using fcsread
% Check initial clustering for 42973AZ.20121115.114853.6744, 42973AZ.20121118.103628.6947 (mono vs poly)