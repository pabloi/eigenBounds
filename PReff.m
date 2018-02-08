function pr=PReff(C)
pr=sum(diag(C))^2/norm(C,'fro')^2;

end