function CS=fcst(X,N,dim)
%Computes a sine/cosine decomposition of the Fourier transform of REAL data
%On real data, the fft is hermitic (conjugate symmetric), which implies
%redundancy. Further, this representation loses phase information if we
%only look at the energy/power spectral density (a function of abs(fft)).
%Here we save this information by representing the fft differently: the
%first half will be real(fft(X)), for f>=0, and the second half will be
%imag(fft(X)), also for f>0 (strict, as DC component has no imag part). 

%There is some relation to the hilbert transform, I think

dim=1; %Default, ignoring dim argument for now. TODO: change to first non-singleton
%if nargin<3
%    dim=1; %DEfault dimension to act. 
%end
if nargin<2
    N=size(X,dim);
end
F=fft(X,N,dim);
firstHalfLength=floor(N/2)+1;
secondHalfLength=N-firstHalfLength;
C=real(F(1:firstHalfLength,:)); 
S=-imag(F(firstHalfLength+1:N,:)); %Taking imag of second half, which already doesn't contain DC or the f=fsamp/2 component
C(1,:)=C(1,:)/sqrt(2); %This is so all the coefficients are equally scaled (all represent sqrt(energy) in module)
CS=sqrt(2)*cat(dim,C,S)/sqrt(N); %Normalized to be unitary 


%Check:
%C(1,:)=sqrt(2)*C(1,:);
%realF=cat(dim,C,flipud(C(2:secondHalfLength+1,:)));
%imagF=cat(dim,[zeros(1,size(S,2));flipud(S);zeros(firstHalfLength-secondHalfLength-1,size(S,2))],-S);
%newF=realF+1i*imagF;

%figure
%subplot(1,2,1)
%hold on
%plot(real(F(:,1)))
%plot(realF(:,1))
%subplot(1,2,2)
%hold on
%plot(imag(F(:,1)))
%plot(imagF(:,1))