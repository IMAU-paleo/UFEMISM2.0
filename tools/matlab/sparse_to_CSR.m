function [m, n, A_ptr, A_index, A_val] = sparse_to_CSR( A)

m = size(A,1);
n = size(A,2);

A_ptr   = zeros( m+1,1);
A_index = zeros( nnz(A),1);
A_val   = zeros( nnz(A),1);

k = 1;
A_ptr( 1) = 1;
for i = 1: m
  for j = 1: n
    if (nnz(A(i,j)==1))
      A_index( k) = j;
      A_val(   k) = A(i,j);
      k = k+1;
    end
  end
  A_ptr( i+1) = k;
end

end