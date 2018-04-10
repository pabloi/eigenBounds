function [upper,lower]=PR2VAFbound(PR,M)
%Check: PR<M

if PR==M
    upper=[1:M]/M;
    lower=upper;
else
    upper=ones(1,M);

    %First: find the two-value eigenvalue distribution with given PR
    Mmax=min(floor(PR),M-1);%For i larger than floor(PR), the upper bound is 1.
    for i=1:Mmax 

       %Need to solve:   %1/PR=(i*a^2 + (M-i)*b^2);
       %Subject to:   %b=(1-i*a)/(M-i);
       %Results in: %1/PR=(i*a^2 + (1-i*a)^2/(M-i));
       a=roots([i+i^2/(M-i) -2*i/(M-i) 1/(M-i)-1/PR]);
       b=(1-i*a)/(M-i);
       %Keep the positive solution only:
       aux=a>0 & b>=-1e-9/M & a>=b;
       a=a(aux);
       b=b(aux);

       %Check results:   
       Prinv=(i*a^2 + (M-i)*b^2);
       if abs(Prinv*PR-1)>1e-6
           error('')
       end

       %Second, find upper and lower bounds for each dimension:
       upper(i)=i*a;

       if i==1
           a1=a;
           b1=b;
       end
    end

    lower=min([1-b1*[M-1:-1:0];a*[1:M]]);
    %Lower should be i*a for the first few elements

end