function [energyDecomp,periodicSpectrumAlpha,Nharmonics]=estimateTemporalCharacteristics(X,T)
%Estimates energy decomposition and periodic component spectrum, for use in
%fourierVAFbound
%INPUTS:
%X NxM matrix, assumed to be M independent timeseries of length N
%T periodicity of the signal

if nargin<2 || isempty(T) %Period not given, estimating
    error('Unimplemented')
end

[N,M]=size(X);
Nperiods=floor(N/T);
X=X(1:Nperiods*T,:);
energyTotal=sum(X.^2);
DC=mean(X);
energyDecomp(1,:)=(Nperiods*T*DC.^2);
X=X-DC;
Xp=reshape(X,T,Nperiods,M);
periodic=(mean(Xp,2));
T2T=Xp-periodic;
T2T=reshape(T2T,Nperiods*T,M);
periodic=squeeze(periodic);

energyDecomp(3,:)=(sum(T2T.^2));
energyDecomp(2,:)=(Nperiods*sum(periodic.^2));

energyDecomp=energyDecomp./energyTotal;
%Check: energy adds up
if any(abs(sum(energyDecomp)-1)>1e-5)
    error('Energy doesn''t add up')
end
energyDecomp=mean(energyDecomp,2);

F=abs(fft(periodic)).^2;
F=F(2:ceil(T/2),:); %Only first half of spectrum
%Estimate Nharmonics:
aux=cumsum(F)/sum(F);
Nharmonics=find(aux>.98,1,'first');
y=log(mean(F(1:Nharmonics,:),2));
x=log([1:length(y)]');
%plot(x,y)
p=polyfit(x,y,1);
periodicSpectrumAlpha=-p(1);

