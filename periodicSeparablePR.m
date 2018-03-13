function PR=PRboundPeriodicSeparable(f,a,period,N)
G=f(1)^2 +2*sum(f(2:end).^2);
Gplus=sum(f.^2);
T=numel(f);
M=numel(a);
Gi=sum([0:T-1]'.*f(:).^2);

%Check that 2*T>=period
%Check N<= period*(M-.5);

A=sum(a(2:end-1).^2 .* (N-[1:M-2]*period));

PR=.5*N^2*f(1)^2*a(1)^2 / (G*A + N*a(1)^2*Gplus + (a(end)^2-a(1)^2)*Gi);

end