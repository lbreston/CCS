function Y=pwCCS(X,l,tau,dim)
triSymOut={@(c)num2cell(c,2),@(d,e)sub2ind(e,[d(1) d(2)],[d(2) d(1)])};
Y=tuplefun(@(w,z)CCS(w,z,l,tau,dim),X,1,false,triSymOut);
end

