function pr=PRbound(autocorr)
%Implements Eq. 16 from Gao et al. 2017
%autocorr is an length N vector, that represents the auto-correlation starting at t=0, and until
%t=N. It is assumed to be symmetric around 0.

%pr= 2*numel(autocorr)*autocorr(1)^2 / (2*sum(autocorr.^2));

%My exact computation:
N=numel(autocorr);
pr= N / (1 + 2*sum([N-1:-1:1]'.*autocorr(2:end).^2)/(N*autocorr(1)^2));

end