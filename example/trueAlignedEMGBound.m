%% 
%load ../data/EMG2.mat %Need to load aligned EMG, or post-hoc align
%%
D=EMG2{5}.Data; %40s (4k samples) of data for 30 muscles
D=D./sqrt(sum(D.^2,1));
C=D'*D;
%% PCA results
l=eig(C);
%% Empirical power spectral density
F=fft(D,[],1);
s=sum(abs(F).^2,2);
G=dct(D,[],1);
t=sum(abs(G).^2,2);

%%
fh=figure;
M=length(l);
plot(cumsum(sort(l,'descend'))/sum(l),'DisplayName','Eigenvalues')
hold on
b=cumsum(sort(s,'descend'))/sum(s);
plot(b,'DisplayName','PSD bound')
plot(b(1:M)+(1-b(1:M)).*[1:M]'/M,'DisplayName','Experimental DFT bound')
c=cumsum(sort(t,'descend'))/sum(t);
plot(c,'DisplayName','DCT bound');
plot(c(1:M)+(1-c(1:M)).*[1:M]'/M,'DisplayName','Experimental DCT bound')
axis([1 30 0 1])