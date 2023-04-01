function A = CSR_to_sparse( A)

nnz = A.ptr(end)-1;

Ai = zeros( nnz, 1);
Aj = zeros( nnz, 1);
Av = zeros( nnz, 1);

k = 0;
for i = 1: A.m
  for ii = A.ptr( i): A.ptr( i+1)-1
    j = A.ind( ii);
    v = A.val(   ii);
    k = k+1;
    Ai( k) = i;
    Aj( k) = j;
    Av( k) = v;
  end
end

A = sparse(Ai,Aj,Av,A.m,A.n);

end