%% Example of Fourier-based bound for autocorrelated data
close all
%TODO: generate data NOT by multiplying  white noise with the filter in
%frequency space, but by generating the frequency representation directly:
%this way we can have the same expected value, but have fixed variance
%(current version has signal-dependent noise, proportional)
%% Construct: a data matrix with independent but autocorrelated columns
N=1e4;
M=30;
X=(randn(N,M)); %M realizations of N-sample white noise
for j=1:M
sigma=20*(1+5*abs(randn));
H=exp(-[0:N-1].^2'/(2*sigma.^2))+.1; %Frequency-space shaping of signal
Ga=ones(N,1);% White noise
Ga(20:20:N)=Ga(20:20:N)+10./sqrt([1:20:N]'); %1/f signal
Ga=Ga.*H;
Ga((N/2)+2:N)=conj(flipud(Ga(2:(N/2)))); %Ensuring a real-valued filter
%Ga(N/2+1)=0;
Ga(1)=abs(10*randn); %Adding more DC
G(:,j)=Ga;

X(:,j)=ifft(fft(X(:,j),[],1).*Ga); %Filtering
%g=ifft(G(:,j));
%g=g(1:200); %200-th order FIR approximation;
end
X=X./sqrt(sum(X.^2,1)); %Normalizing energy per column: this is not needed, but if not done the bounds are trivial for the filters used (the dominant difference between columns is the energy, not their autocorrelations)
G=G./sqrt(sum(G.^2,1));
sqrt(sum(X.^2,1))
sqrt(sum(abs(G).^2,1))
%% Plots
figure('Units','Normalized','OuterPosition',[0 0 1 1])
subplot(2,2,1) %Data
plot(X)
title('Data')

subplot(2,2,2) %Data autocorrelation
hold on
r=0;
for j=1:M
rj=autocorr(X(:,j),'NumLags',999);
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
G1=abs(G).^2;
b=VAF(sum(G1,2)); %Fourier-based bound
plot(b,'LineWidth',2,'DisplayName','*True* column spectrum based bound')
F=fft(X,[],1);
c=VAF(sum(abs(F).^2,2));
plot(c,'LineWidth',2,'DisplayName','Empirical column spectrum based bound')
d=VAF(sum(X.^2));
plot(d,'LineWidth',2,'DisplayName','Empirical column energy based bound')
C=dct(X,[],1);
e=VAF(sum(C.^2,2));
plot(e,'LineWidth',2,'DisplayName','Empirical column DCT based bound')
%plot(d+(1-d).*c(1:length(d))','LineWidth',2,'DisplayName','Experimental') %This works as bound IF all columns have same energy (ie., d is a straight line)
plot([1:M]/M+(1-[1:M]/M).*c(1:M)','LineWidth',2,'DisplayName','Experimental 2') 
axis([1 M 0 1])
title('Eigenvalues')
legend('Location','SouthEast')
ylabel('Variance Accounted For (%)')
xlabel('Dimensions')

subplot(2,2,4)
plot(sum(G1,2),'LineWidth',2,'DisplayName','True spectral density')
hold on
plot(sum(abs(F).^2,2)/N,'LineWidth',2,'DisplayName','Empirical spectral density')
set(gca,'YScale','log')
legend