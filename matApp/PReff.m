function pr=PReff(C)
if numel(C)==length(C) %Assuming it was vector of eigenvalues
pr=sum(eigenvalues)^2/sum(eigenvalues.^2);
else
pr=trace(C)^2/norm(C,'fro')^2;
%Equivalent to:
%    pr=PRdef(eig(eigenvalues));
end

end
