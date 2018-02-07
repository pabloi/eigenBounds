function [Ct,f]=toeplitize(C)
%Computes the LS best-fit toeplitz matrix Ct that approximates C
%Useful to find the 'best' stationary approximation to a covariance matrix
N=size(C,1);
t=nan(N,1);
for k=1:N
    t(k)=mean(diag(C,mod(k-1,N)));
end
Ct=toeplitz(t);
f=t;
end