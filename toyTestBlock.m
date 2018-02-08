%% Toy test:
M=25;
tau=30;
N=tau*21;
x=nan(M,N);
x=(sin(2*pi*[0:N-1]/tau+.5*pi*[0:M-1]'/M + .03*randn(size(x)))) + .5*randn(size(x));

%%
aux=x-mean(x,2);
C=aux'*aux;

%%
[Ct2,f2]=toeplitize(C);
[Ct,f,a]=blockToeplitize(C,tau);

%%
actualPR=PRdef(C) %Definition of PR for the original matrix
disp('Should be smaller than:')

%Block approximation:
looseBlockBound=PRbound(f)*PRbound(a) %Gao bound
exactBlockBound=PRboundExact(f)*PRboundExact(a) %My bound
blockToeplitzPR=PRdef(Ct) %Definition of PR for the block-toeplitz matrix


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