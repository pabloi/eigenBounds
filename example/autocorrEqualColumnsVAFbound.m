%% Example of Fourier-based bound for autocorrelated data

%% Construct: a data matrix with independent but autocorrelated columns
N=1e3;
M=30;
dc=.11;
v=3+randn(1,M);
X=(randn(N,M)+dc).*v;

sigma=20;
G=exp(-[0:N-1]'/sigma); %Frequency-space shaping of signal

G(2:end)=G(2:end)+conj(flipud(G(2:end))); %Ensuring a real-valued filter

X=ifft(fft(X,[],1).*G); %Filtering
g=ifft(G);
g=g(1:200); %200-th order FIR approximation;
% for j=1:M
% X(:,j)=conv(X(:,j),g,'same'); %Alternative, causal, filtering
% end

%% Plots
figure('Units','Normalized','OuterPosition',[0 0 1 1])
subplot(2,2,1) %Data
plot(X)
title('Data')

subplot(2,2,2) %Data autocorrelation
hold on
r=0;
for j=1:M
rj=autocorr(X(:,j),'NumLags',200);
plot(rj)
r=r+rj*sum(X(:,j).^2);
end
plot(g/g(1),'k','LineWidth',2)
plot(r/r(1),'r','LineWidth',2)
title('Autocov')
xlabel('Lag (samples)')

subplot(2,2,3) %Bounds
l=VAF(X'*X); %Computing X'*X because it is cheaper
plot(l,'LineWidth',2,'DisplayName','Empirical eigenvalues')
hold on
G1=G;
G1(1)=G1(1)*(dc*sqrt(N)); %G is the filter spectrum. But the original signal spectrum was flat except for the DC component. This is the theoretical DC component after filtering.
b=VAF(abs(G1).^2); %Fourier-based bound
plot(b,'LineWidth',2,'DisplayName','*True* column spectrum based bound')
F=fft(X,[],1);
c=VAF(sum(abs(F).^2,2));
plot(c,'LineWidth',2,'DisplayName','Empirical column spectrum based bound')
d=VAF(sum(X.^2));
plot(d,'LineWidth',2,'DisplayName','Empirical column energy based bound')
C=dct(X,[],1);
e=VAF(sum(C.^2,2));
plot(e,'LineWidth',2,'DisplayName','Empirical column DCT based bound')
axis([1 M 0 1])
title('Eigenvalues')
legend('Location','SouthEast')
ylabel('Variance Accounted For (%)')
xlabel('Dimensions')