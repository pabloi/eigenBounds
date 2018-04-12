%% Example of Fourier-based bound for autocorrelated data
close all
%% First, create some expected spectral density for each column
N=1e4;
M=30;
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
end

%% Construct: 1e5 repetitions of matrix with the above *true* spectral density, compute eigenvalues
Nsamp=1e3;
meanF2=0;
meanS=0;
meanT=0;
allL=nan(Nsamp,M);
allExcessSVAF=nan(Nsamp,M);
allSvaf=nan(Nsamp,N);
for k=1:Nsamp
    X=(randn(N,M)); %M realizations of N-sample white noise
    for j=1:M
    X(:,j)=ifft(fft(X(:,j),[],1).*G(:,j)); %Filtering
    end
    %Compute empirical density:
    F=fft(X,[],1);
    F2=abs(F).^2; %Sample spectral density
    meanF2=meanF2+F2/Nsamp;
    S=F2/sum(F2(:)); %Normalizing each sample spectral density individually
    meanS=meanS+S/Nsamp;
    s=VAF(sum(S,2)); %This is the same as VAF(sum(F2,2)), normalization doesnt matter
    allSvaf(k,:)=s;
    %Compute eigenvalues:
    l=VAF(X'*X);
    allL(k,:)=l;
    allExcessSVAF(k,:)=l-s(1:M); %This can be negative for any one realization, but its average should be strictly positive
    %Compute trace (for normalization purposes)
    t=norm(X,'fro')^2;
    meanT=meanT+t/Nsamp;
end

meanL=mean(allL);
stdL=std(allL);
meanSvaf=mean(allSvaf);
stdSvaf=std(allSvaf);
meanE=mean(allExcessSVAF);
stdE=std(allExcessSVAF);
%% Plots
figure('Units','Normalized','OuterPosition',[0 0 1 1])
subplot(2,2,1) %Data
p1=plot(sum(meanF2,2)/N,'LineWidth',3,'DisplayName','Empirical mean spectral density');
hold on
p2=plot(sum(abs(G).^2,2),'LineWidth',1,'DisplayName','Expected spectral density');
set(gca,'Yscale','log')
title('Spectral density (DFT)')
legend
axis([1 6e2 .1 1e4])

subplot(2,2,2)
plot(sum(meanS,2),'LineWidth',3,'DisplayName','Empirical normalized mean spectral density','Color',p1.Color)
hold on
plot(sum(abs(G).^2,2)/norm(abs(G),'fro')^2,'LineWidth',1,'DisplayName','Expected normalized spectral density','Color',p2.Color)
set(gca,'Yscale','log')
title('Normalized spectral density (DFT)')
legend
axis([1 6e2 1e-6 1])
%These plots appear to match exactly, somthing I didn't anticipate. May be
%there is point convergence (thus the appearance) but no distribution
%convergence (thus the slight differences between the expectation bound and
%the 

subplot(2,2,3) %Bounds
hold on
patch([1:M M:-1:1],[meanL+stdL fliplr(meanL-stdL)],'k','FaceAlpha',.3,'EdgeColor','none')
plot(meanL,'LineWidth',3,'DisplayName','Empirical mean VAF','Color','k')

patch([1:M M:-1:1],[meanSvaf(1:30)+stdSvaf(1:30) fliplr(meanSvaf(1:30)-stdSvaf(1:30))],p1.Color,'FaceAlpha',.3,'EdgeColor','none')
plot(meanSvaf,'LineWidth',3,'DisplayName','Empirical mean VAF bound','Color',p1.Color)


%c=VAF(sum(meanF2,2)); %This should match b exactly because F2 should converge on the average to G.
%plot(c,'LineWidth',2,'DisplayName','Empirical column spectrum based bound')
G1=abs(G).^2;
b=VAF(sum(G1,2)); %Fourier-based bound
plot(b,'LineWidth',1,'DisplayName','*True* column spectrum based bound','Color',p2.Color)
c=[1:M]/M+(1-[1:M]/M).*b(1:M)'; %Experimental bound: take the spectral bound, and assume that independently of it, 1/M of the excess energy is explained by every dimension added
p3=plot(c,'LineWidth',1,'DisplayName','Experimental');
axis([1 M 0 1])
legend('Location','SouthEast')
title('VAF')
ylabel('Variance Accounted For (%)')
xlabel('Dimensions')


subplot(2,2,4)
hold on
patch([1:M M:-1:1],[meanE(1:30)+stdE(1:30) fliplr(meanE(1:30)-stdE(1:30))],p1.Color,'FaceAlpha',.3,'EdgeColor','none')
plot(meanE,'LineWidth',3,'DisplayName','Empirical mean excess VAF bound','Color',p1.Color)
patch([1:M M:-1:1],[meanL+stdL-b(1:M)' fliplr(meanL-b(1:M)'-stdL)],p2.Color,'FaceAlpha',.3,'EdgeColor','none')
plot(meanL-b(1:30)','LineWidth',3,'DisplayName','Excess expected VAF bound','Color',p2.Color)
patch([1:M M:-1:1],[meanL+stdL-c(1:M) fliplr(meanL-c(1:M)-stdL)],p3.Color,'FaceAlpha',.3,'EdgeColor','none')
plot(meanL-c(1:30),'LineWidth',3,'DisplayName','Excess experimental bound','Color',p3.Color)

title('VAF excess')
axis([1 M 0 .5])
ylabel('Empirical VAF above bound')
xlabel('Dimensions')
legend('Location','NorthWest')