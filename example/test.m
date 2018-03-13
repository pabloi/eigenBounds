%Test:

%% Load some data
load ./data/alignedEMG.mat
%%
aux=allEMG{1,1}.Data;
aux=aux-nanmean(aux);
C=aux*aux';
C(isnan(C))=0;

%% 
[Ct,f]=toeplitize(C);
pr=PRbound(f)
pr1=PRdef(eig(Ct))

%% Empirical participation ratio:
[p,c,a]=pca(C);
empPR=sum(a)^2/sum(a.^2);