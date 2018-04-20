function newBound=improveBound(bound,columnEnergy)
%Given a VAF bound computed, and a known distribution of column energies,
%it reduces slack of the bound (i.e. increases the bound or leaves it the
%same) while provably staying a lower bound.

%If columnEnergy is a scalar and interger, it is taken to be the rank of
%the matrix, and column energies are taken to be equal. Otherwise the rank 
%is taken to be length(columnEnergy)
if numel(columnEnergy)==1
    columnEnergy=ones(columnEnergy,1)/columnEnergy; %Equal column energies
else
    columnEnergy=columnEnergy./sum(columnEnergy); %Making sure they add up to 1
end

cE=sort(columnEnergy,'descend');
dB=diff(bound);
newBound(1)=max(bound(1),cE(1));
rank=length(cE);
for i=1:rank-1
    minimumIncrease=(1./(rank-i))*(1-newBound(i));
    if dB(i)<=minimumIncrease %Increase in bound was less than 1/M of leftover energy, impossible!
        newBound(i+1)=newBound(i)+minimumIncrease;
    else
        newBound(i+1)=bound(i+1);
    end
end

end