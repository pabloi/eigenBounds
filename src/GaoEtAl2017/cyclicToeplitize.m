function [Ct,f]=cyclicToeplitize(C)
%Not symmetric necessarily

%N=size(C,1);
%t=ifft(sum(abs(fft(fft(C,[],1),[],2)).^2));
%t=real(t); %Only necessary because of numeric issues
%f=t*trace(C)/(N*t(1));
%Ct=toeplitz(f);

%Alt:
N=size(C,1);
t=nan(N,1);
r=nan(1,N);
for k=0:N-1
    t(k+1)=mean([diag(C,k); diag(C,N-k)]);
    r(k+1)=mean([diag(C,-k); diag(C,k-N)]);
end
f=t;

Ct=toeplitz(f,r);

%Alternative computation: call circulantize()
