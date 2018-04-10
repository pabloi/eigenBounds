%% Example of Fourier-based bound for autocorrelated data

%% Construct: a data matrix with independent but autocorrelated columns
N=1e4;
M=30;
dc=.1;
X=(randn(N,M)+dc);
for j=1:M
sigma=10*(1+5*abs(randn));
Ga=exp(-[0:N-1]'/sigma); %Frequency-space shaping of signal
Ga(2:end)=Ga(2:end)+conj(flipud(Ga(2:end))); %Ensuring a real-valued filter
G(:,j)=Ga;

X(:,j)=ifft(fft(X(:,j),[],1).*G(:,j)); %Filtering
g=ifft(G(:,j));
g=g(1:200); %200-th order FIR approximation;
end
X=X./sqrt(sum(X.^2,1)); %Normalizing energy per column: this is not needed, but if not done the bounds are trivial for the filters used (the dominant difference between columns is the energy, not their autocorrelations)

%% Plots
figure('Units','Normalized','OuterPosition',[0 0 1 1])
subplot(2,2,1) %Data
plot(X)
title('Data')

subplot(2,2,2) %Data autocorrelation
hold on
r=0;
for j=1:M
rj=autocorr(X(:,j),'NumLags',249);
plot(rj)
r=r+rj*sum(X(:,j).^2);
end
%plot(g/g(1),'k','LineWidth',2)
plot(r/r(1),'r','LineWidth',2)
title('Autocov')
xlabel('Lag (samples)')

subplot(2,2,3) %Bounds
l=VAF(X'*X); %Computing X'*X because it is cheaper
plot(l,'LineWidth',2,'DisplayName','Empirical eigenvalues')
hold on
G1=G;
G1(1,:)=G1(1,:)*(dc*sqrt(N)); %G is the filter spectrum. But the original signal spectrum was flat except for the DC component. This is the theoretical DC component after filtering.
b=VAF(sum(abs(G1).^2,2)); %Fourier-based bound
plot(b,'LineWidth',2,'DisplayName','*True* column spectrum based bound')
F=fft(X,[],1);
c=VAF(sum(abs(F).^2,2));
plot(c,'LineWidth',2,'DisplayName','Empirical column spectrum based bound')
d=VAF(sum(X.^2));
plot(d,'LineWidth',2,'DisplayName','Empirical column energy based bound')
C=dct(X,[],1);
e=VAF(sum(C.^2,2));
plot(e,'LineWidth',2,'DisplayName','Empirical column DCT based bound')
plot(d+(1-d).*c(1:length(d))','LineWidth',2,'DisplayName','Experimental')
axis([1 M 0 1])
title('Eigenvalues')
legend('Location','SouthEast')
ylabel('Variance Accounted For (%)')
xlabel('Dimensions')