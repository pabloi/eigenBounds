function p=PR(x)
if numel(x)==length(x) %Vector
    p=sum(x).^2/sum(x.^2);
else %Matrix
    %This works but is inefficient:
    %x=eig(x); %Needs to be square
    %v=PR(x);
    %This is efficient:
    p=tr(x)^2/norm(x,'fro')^2;
end
end