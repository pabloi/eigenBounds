function X=ifcst(CS,N,dim)
dim=1;
%if nargin<3
%    dim=1;
%end
if nargin<2
    N=size(CS,dim);
end
firstHalfLength=floor(N/2)+1;
secondHalfLength=N-firstHalfLength;
C=CS(1:firstHalfLength,:);
C(1,:)=sqrt(2)*C(1,:);
S=CS(firstHalfLength+1:end,:);
realF=cat(dim,C,flipud(C(2:secondHalfLength+1,:)));
imagF=cat(dim,[zeros(1,size(S,2));flipud(S);zeros(firstHalfLength-secondHalfLength-1,size(S,2))],-S);
newF=realF+1i*imagF;
newF=sqrt(N)*newF/sqrt(2);
X=ifft(newF,N,dim,'symmetric'); %Symmetric flag is here to avoid complex results because of numerical issues
