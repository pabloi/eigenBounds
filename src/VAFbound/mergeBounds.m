function mixedBound=mergeBounds(bound1,bound2)
%This function takes two LOWER bounds for a VAF/eigenvalue curve and merges
%them into a single,better bound. It does this by first taking the max
%across the two bounds, and then imposing a concavity of bound constraint

M1=length(bound1);
M2=length(bound2);
%First, make bounds same size
if M1<M2
    bound2=bound2(1:M1);
else
    bound1=bound1(1:M2);
end
%Second, take max:
mixedBound=max(bound1,bound2);
%Third, impose concavity:
mixedBound=[0 mixedBound];
for i=1:length(mixedBound)-2
    dB=mixedBound(i+1)-mixedBound(i); %This eigenvalue
    project=mixedBound(i)+dB*[1:length(mixedBound)-i]; %Projection if all following eigenvalues are the same as this
    part=mixedBound(i+1:end); %Current estimation of eigenvalue cumulative sum
    if any(project<part) %If the projection falls below the current estimate, need to correct
        [m,idx]=max(part-project);
        m
        idx
        newEig=(part(idx)-mixedBound(i))/idx
        mixedBound(i+1:i+idx)=mixedBound(i)+newEig*[1:idx];
    end
end
mixedBound=mixedBound(2:end);
end