function [PR,looseBound]=PRboundPeriodicSeparable(f,a,N)
G=sum(f.^2);
period=numel(f); %If mod(T,2)==0, then f(end) should be 0.
offset=ceil(period/2);
T=period-offset+1;
fplus=f(offset:end);
Gplus=.5*fplus(1)^2+sum(fplus(2:end).^2); %This should be literally G/2
M=numel(a);
Gi=sum([0:T-1]'.*fplus.^2);

if nargin<4
    N=period*(M-1); %Assuming we are measuring an exact number of periods
end

%Check that T>=period
%Check N<= period*(M-.5);

A=sum(a(2:end-1).^2 .* (N-[1:M-2]*period));
numerator=.5*N^2*fplus(1)^2*a(1)^2 ;
firstPeriodContribution=a(1)^2*(N*G/2-Gi);
middlePartContribution=G*A;
lastPeriodContribution=a(end)^2*Gi;
PR=numerator/(firstPeriodContribution + middlePartContribution + lastPeriodContribution);

%Sanity check:
%autocorr=reshape(f.*a,1,numel(f)*numel(a));
%autocorr=autocorr(offset:offset+N-1);
%summands=[N:-1:1].*autocorr.^2;
%firstCheck=.5*summands(1)+sum(summands(2:T)); %==firstPeriod
%middleCheck=sum(summands(T+1:end-T+1)); %==middlePart
%lastCheck=sum(summands(end-T+2:end)); %==lastPeriod
am=mean(a(2:end));
looseBound=N*fplus(1)^2*a(1)^2 / (fplus(1)^2*(a(1)^2-am^2)+G*am^2*(M));
end