function [T] = toeplitizeV2(A)
%TOEPLITIZE Summary of this function goes here
%   Detailed explanation goes here

%First: 
M=max(size(A));
m=min(size(A));
D=circulantize(A,[],true);

%Third: re-normalize
if M==size(A,2)
    aux=toeplitz(m:-1:1,[m*ones(1,M-m), m:-1:1]);
else
    aux=toeplitz([m*ones(1,M-m), m:-1:1],m:-1:1);
end
T=m*D./aux;

%Q: can this whole process be done faster 
%by using a weighted fft with some given sample size (2M-1)? It seems it should 
%(i.e. toeplitizing is akin to a circulantization of 
%zero-padded matrix with amplification of high-freq components)
end