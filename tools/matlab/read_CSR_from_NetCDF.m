function A = read_CSR_from_NetCDF( filename)

% Get matrix size
f = ncinfo( filename);
A.m   = [];
A.n   = [];
A.nnz = [];
for di = 1: length( f.Dimensions)
  if     strcmpi(f.Dimensions(di).Name,'m')
    A.m   = f.Dimensions(di).Length;
  elseif strcmpi(f.Dimensions(di).Name,'n')
    A.n   = f.Dimensions(di).Length;
  elseif strcmpi(f.Dimensions(di).Name,'nnz')
    A.nnz = f.Dimensions(di).Length;
  end
end

% Safety
if isempty(A.m) || isempty(A.n) || isempty(A.nnz)
  error('Couldnt find matrix size in file!')
end

% Read matrix data
A.ptr = ncread( filename, 'ptr');
A.ind = ncread( filename, 'ind');
A.val = ncread( filename, 'val');

% Convert to Matlab sparse matrix
A = CSR_to_sparse( A);

end