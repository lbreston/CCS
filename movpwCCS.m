function Y = movpwCCS(X,l,tau,dim,k,o)
Y=movfun(@(X1,X2)pwCCS({X1,X2},l,tau,dim),X,k,o);
end

