function [Ct,f,a]=blockToeplitize(C,blockSize)
%Computes the LS best-fit block-toeplitz matrix Ct that approximates C
%Ct is a matrix made of blocks of size blockSize
%each block is the product of a weight with a symmetric toeplitz (the same
%for all blocks)
%Thus, the structure of the matrix Ct is the Kroneker product of a two
%smaller symmetric toeplitz matrices.

%First: check that blockSize is a divisor of the size of C, if not, crop C
N=size(C,1);
k=floor(N/blockSize);
if k*blockSize<N
    warning('Block size is not a divisor of matrix size, cropping matrix')
    N=k*blockSize;
    C=C(1:N,1:N);
end

%Second: do the factorization
%2a: average across diagonal blocks
subblock=zeros(blockSize,blockSize,k);
for j=1:k
    offsetK=(j-1)*blockSize;
    for i=1:(k-j+1)
         offsetRow=(blockSize*(i-1)); %Averaging only along upper triangle, since it is symmetric
         subblock(:,:,j)= subblock(:,:,j)+C(offsetRow+[1:blockSize],offsetRow+[1:blockSize]+offsetK);
    end
    subblock(:,:,j)=subblock(:,:,j)/(k-j+1);
end

%2b: factorize the blocks
%figure; imagesc(reshape(subblock,blockSize,blockSize*k));axis equal;
subblock=reshape(subblock,blockSize^2,k);
[a,F]=pca(subblock,'Centered','off');
F=F(:,1);
a=a(:,1);

%2c: toeplitize the main subblock
[Ft,f]=toeplitize(reshape(F,blockSize,blockSize));

%Third: reconstruct the block-toeplitz matrix
Ct=nan(size(C));
for j=1:k
    offsetK=(j-1)*blockSize;
    for i=1:(k-j+1)
         offsetRow=(blockSize*(i-1)); %Averaging only along upper triangle, since it is symmetric
         Ct(offsetRow+[1:blockSize],offsetRow+[1:blockSize]+offsetK)=Ft*a(j);
    end
end
CtT=Ct';
Ct(isnan(Ct))=CtT(isnan(Ct));

end