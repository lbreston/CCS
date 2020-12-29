function Y = movCCS(X,l,tau,dim,k,o)
Y=movfun(@(w,z)CCS(w,z,l,tau,dim),X,k,o);
end