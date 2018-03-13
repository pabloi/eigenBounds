function [Ct,g,f,a]=periodicToeplitize(C,period,separableFlag)
%Computes the LS best-fit (residual frobenius norm) pseudo-periodic toeplitz matrix Ct that approximates C

%The pseudo-periodic toeplitz matrix is defined as a matrix that satisfies:
%Ct(i,i+k*period+j)=Ct(i,i+k*period-j)=Ct(i,i-k*period+j)=Ct(i,i-k*period-j) for all i
%Ct(i,i+k*period+j)=g(k,j); (j<= period/2)
%Notice that this is a symmetric toeplitz matrix.
%If the separableFlag is set, then we also enforce g(k,j)=f(j)*a(k).
%INPUT:
%period needs to be an integer, and there is NO need that M*period=size(C,1) for some integer M
%separableFlag is [0,1,2]

%Notice that if period is even then BOTH g(k,period/2) and g(k+1,period/2)
%define the same diagonals (overlap of bands). 
%Two solutions: 
%1) arbitrarily assign Ct(i,i+(k-.5)*period)=g(k,period/2) AND Ct(i,i+(k+.5)*period)=g(k+1,period/2)
%I think this won't work, because then I lose the property PR(Ct)>=PR(C)
%2) average all those diagonals together, and enforce g(k,period/2)=g(k+1,period/2) for all k
%This second solution means that the anti-period diagonals are all the same
%and they are separable iff g(k,period/2)=0 OR a(k)=a for all k;

if nargin<3
    separableFlag=0;
end
N=size(C,1);
M=ceil(N/period -.5)+1; %Number of periods present in data, equal to bands that need estimating
T=floor(period/2)+1;
g=nan(T,M);
f=[];
a=[];

for k=1:M
    for j=1:T
        %Get difference in subscripts:
        gap=(k-1)*period+(j-1);
        mirrorGap=(k-1)*period-(j-1);
        %Make list of all subscript combinations possible:
        aux=[1:N]';
        relevantSubscriptsI=[aux;aux];
        relevantSubscriptsJ=[aux+gap; aux+mirrorGap];
        %Remove out-of-range subscripts:
        bad=relevantSubscriptsJ<1 | relevantSubscriptsJ>N;
        relevantSubscriptsI=relevantSubscriptsI(~bad);
        relevantSubscriptsJ=relevantSubscriptsJ(~bad);
        %Get elements:
        data=C(sub2ind([N,N],relevantSubscriptsI,relevantSubscriptsJ));
        %Average:
        g(j,k)=mean(data(:));
    end
end
if mod(period,2)==0 && ~separableFlag %Even period, need to fix anti-period diagonals
    %Get all anti-period diagonal elements:
    aux=[1:N]';
    relevantSubscriptsI=repmat(aux,M,1);
    antiPeriodGaps=[1:M]*period-period/2;
    relevantSubscriptsJ=reshape(aux+antiPeriodGaps,N*M,1);
    %Remove out-of-range subscripts:
     bad=relevantSubscriptsJ<1 | relevantSubscriptsJ>N;
     relevantSubscriptsI=relevantSubscriptsI(~bad);
     relevantSubscriptsJ=relevantSubscriptsJ(~bad);
     %Get elements:
     data=C(sub2ind([N,N],relevantSubscriptsI,relevantSubscriptsJ));
     %Average:
     g(T,:)=mean(data(:));
end

%Re-construct Ct:
Taux=T-1+mod(period,2); %Equal to T if mod(period,2)==1, T-1 otherwise;
if separableFlag>0
    warning('This is not right: it computes a factorization, but not one that upperbounds PR')
    %(I think). We should do a weighted average, where the weight is
    %the number of elements of the matrix that were originally averaged to
    %compute each element of g initially.
    g(T,:)=0;
    f=nanmean(g,2); %Averaging over periods
    a=nanmean(g,1)/nanmean(g(:)); %The normalization is to preserve mean(g) = mean(f*a)
    if separableFlag==2 %Assuming that a(2:end)=a, and only a(1) can be different
        a(2:end)=mean(a(2:end)); %I think this works
    end
    g=f.*a;
    f=[flipud(f(2:Taux));f];
    aux=reshape(f.*a,1,numel(f)*numel(a));
    firstRow=aux(Taux:end);
else
    middlePart=[flipud(g(2:Taux,2:end));g(:,2:end)];
    firstRow=[g(:,1)',reshape(middlePart,1,numel(middlePart))];
end
Ct=toeplitz(firstRow(1:N)); %Discarding part of final band to maintain matrix size