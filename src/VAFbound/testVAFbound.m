%% Test VAF decomp
M=30;
realDC=2+2*[1:M]/M; %Different levels of DC
T=100;
phi=2*pi*[0:T-1]'/T;
alpha=1.1;
fourierSeriesABS=1./[1:5]'.^alpha;
for i=1:M
    realPeriodic(:,i)=sin([1:5].*phi+100*pi*rand(1,5)) * (5*fourierSeriesABS+.5*randn(size(fourierSeriesABS))); %Periodic component: superposition of 4
end
Nperiods=20;
realNoise=randn(T,Nperiods,M);
X=reshape(realDC,1,1,M)+reshape(realPeriodic,T,1,M)+realNoise;
X=reshape(X,T*Nperiods,M);

%expectedDCenergy= Nperiods*T*mean(realDC.^2);
%expectedPeriodicEnergy=
%expectedNoiseEnergy=
%%
[energyDecomp,periodicSpectrumAlpha,Ncomponents]=estimateTemporalCharacteristics(X,T);
VAF_emp=fourierVAFbound(energyDecomp,periodicSpectrumAlpha,Ncomponents);
%VAF_theo=fourierVAFbound(energyDecomp,2.2,5);
 
 %% 
 figure
 plot(VAF_emp,'DisplayName','Spectrum-based bound')
 hold on
 l=eig(X'*X);
 l=sort(l,'descend');
 plot(cumsum(l)/sum(l),'DisplayName','Actual VAF')
 
 %plot(0,mean(realDC),'o')
 plot(0,(energyDecomp(1)),'kx','DisplayName','Estimated DC')
 legend