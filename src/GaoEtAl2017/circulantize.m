function [C] = circulantize(A,antiDiagFlag,noAliasFlag)
%CIRCULANTIZE Summary of this function goes here
%   Detailed explanation goes here
if nargin<2 || isempty(antiDiagFlag)
    antiDiagFlag=false;
end
if nargin<3
    noAliasFlag=false;
end
if noAliasFlag
    M=2*max(size(A))-1;
else
    M=size(A,1); %Assume square
end
F=fft(eye(size(A)),M,1);
F=sqrt(sqrt(size(A,1)/M))*F./sqrt(size(A,1)); %Normalization for unitarity
if antiDiagFlag %This generates a Hankel matrix, rather than Toeplitz
    %C=ifft(ifft(diag(diag(fft(fft(A,M,1),M,2))),[],2),[],1);
    C=F'*diag(diag(F*A*F.'))*conj(F);
else
    %C=ifft(fft(diag(diag(ifft(fft(A,M,1),M,2))),[],2),[],1);
    C=F'*diag(diag(F*A*F'))*F;
end
%C=sqrt(numel(C)/numel(A))*C(1:size(A,1),1:size(A,2));
%Note: inner fft(fft()) and outer ifft(ifft()) result in anti-diagonal circularization
%ifft(fft()) and fft(ifft()) result in diagonal circularization.

end