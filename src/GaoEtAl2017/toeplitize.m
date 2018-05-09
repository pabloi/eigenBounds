function [Ct,f]=toeplitize(C)
%Computes the LS best-fit symmetric toeplitz matrix Ct that approximates C
%Useful to find the 'best' stationary approximation to a covariance matrix
N=size(C,1);
t=nan(N,1);
for k=0:N-1
    t(k+1)=mean([diag(C,k); diag(C,-k)]);
end
Ct=toeplitz(t);
f=t;
end

%Computation through efficient circulantize(): toeplitizev2()