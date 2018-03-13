function VAF=fourierVAFbound(energyDecomp,periodicSpectrumAlpha,Nharmonics)
%Computes a bound on the expected(?) VAF of a matrix with independent
%columns with some expected spectral characteristics given by the inputs.

%INPUTS:
%energyDecomp: Expected eenrgy decomp. A 3x1 vector adding to 1. First element represents DC
%energy, second element represents periodic (centered) signal energy, third
%element represents non-periodic (centered) energy.
%periodicSpectrumAlpha: assumes that the power-spectrum of the periodic component
%follows a 1/f^alpha law and uses it.
%Ncomponents: number of relevant harmonics for the periodic signal
%(number of non-zero elements of the Fourier series)

if nargin<3 || isempty(Nharmonics)
    Nharmonics=10;
end

%Check: if energyDecomp does not sum to 1, normalize
if abs(sum(energyDecomp)-1)>1e-5
    warning('Energy is not normalized. Normalizing.')
    energyDecomp=energyDecomp/sum(energyDecomp);
end

periodicDecomposition=1./[1:Nharmonics].^periodicSpectrumAlpha;
periodicDecomposition=periodicDecomposition./sum(periodicDecomposition);

f(1)=energyDecomp(1);
f(2:2:2*Nharmonics+1)=.5*energyDecomp(2)*periodicDecomposition;
f(3:2:2*Nharmonics+1)=.5*energyDecomp(2)*periodicDecomposition;

VAF=cumsum(sort(f,'descend'));




end