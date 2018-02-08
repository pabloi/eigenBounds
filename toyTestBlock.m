%% Toy test:
M=14;
tau=30;
N=tau*11;
x=nan(M,N);
x=(sin(2*pi*[0:N-1]/tau+.2*pi*[0:M-1]'/M + .03*randn(size(x)))) + .5*randn(size(x));

%%
aux=x-mean(x,2);
C=aux'*aux;

%%
[Ct,f,a]=blockToeplitize(C,tau);

%%
pr=PRbound(f)*PRbound(a)
pr1=PR(eig(Ct))
pr2=PR(eig(C))

%% 
figure
subplot(1,2,1)
imagesc(C)
subplot(1,2,2)
imagesc(Ct)