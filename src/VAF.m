function v=VAF(x)

if numel(x)==length(x) %Vector
    x=sort(x,'descend');
    v=cumsum(x)/sum(x);
else %Matrix
    x=eig(x); %Needs to be square
    v=VAF(x);
end

end