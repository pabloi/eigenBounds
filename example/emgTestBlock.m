%% Toy test:
%load ./data/EMG1.mat
x=EMG1{1}.Data';
x(isnan(x))=0;
blockSize=114; %Empirical
x=x(:,1:blockSize*31);  %Discarding data after walking stopped
%%
aux=x-mean(x,2);
C=aux'*aux;

%%
[Ct2,f2]=toeplitize(C);
%[Ct,f,a]=blockToeplitize(C,114);

%%
actualPR=PRdef(C) %Definition of PR for the original matrix
disp('Should be smaller than:')

%Block approximation:
%looseBlockBound=PRbound(f)*PRbound(a) %Gao bound
%exactBlockBound=PRboundExact(f)*PRboundExact(a) %My bound
%blockToeplitzPR=PRdef(Ct) %Definition of PR for the block-toeplitz matrix


%Non-block approximation:
looseToepBound=PRbound(f2)
exactToepBound=PRboundExact(f2)
toepPR=PRdef(Ct2) %Definition of PR for the toeplitz matrx

%% 
figure
subplot(1,3,1)
title('Original covar')
imagesc(C)
subplot(1,3,2)
title('Block toeplitz covar')
imagesc(Ct)
subplot(1,3,3)
title('Toeplitz covar')
imagesc(Ct2)