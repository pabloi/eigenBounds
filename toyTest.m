%% Toy test:
N=2000;
M=14;
x=nan(M,N);
x(:,1)=randn(M,1);
tau=2; a=exp(-1/tau);
s=2;
for k=2:N
    x(:,k)=a*x(:,k-1)+s*randn(M,1);
end

%%
aux=x-mean(x,2);
C=aux'*aux;

%%
[Ct,f]=toeplitize(C);
pr=PRbound(f)
pr1=PR(eig(Ct))
pr2=PR(eig(C))

%% 
figure
subplot(1,2,1)
imagesc(C)
subplot(1,2,2)
imagesc(Ct)