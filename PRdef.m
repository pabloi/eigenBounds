function pr=PRdef(eigenvalues)
if numel(eigenvalues)==length(eigenvalues) %Vector, assuming it is eigenvalues
    pr=sum(eigenvalues)^2/sum(eigenvalues.^2);
else %MAtrix, assuming it is the full matrix
    pr=PRdef(eig(eigenvalues));
end

end