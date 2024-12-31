function mesh = read_matrix_operators_from_file( mesh, filename)

f = ncinfo( filename);

for gi = 1: length( f.Groups)

  name = f.Groups( gi).Name;

  if ~startsWith( name,{'M_','M2_'})
    continue
  end

  disp(['Reading matrix operator ' name '...'])

  % Read dimension sizes
  A.m = [];
  A.n = [];
  A.nnz = [];
  for di = 1: length( f.Groups( gi).Dimensions)
    switch f.Groups( gi).Dimensions( di).Name
      case 'm'
        A.m   = f.Groups( gi).Dimensions( di).Length;
      case 'n'
        A.n   = f.Groups( gi).Dimensions( di).Length;
      case 'nnz'
        A.nnz = f.Groups( gi).Dimensions( di).Length;
    end
  end

  % Safety
  if isempty(A.m) || isempty(A.n) || isempty(A.nnz)
    error('Couldnt find matrix size in file!')
  end
  
  % Read matrix data
  A.ptr = ncread( filename, [name '/ptr']);
  A.ind = ncread( filename, [name '/ind']);
  A.val = ncread( filename, [name '/val']);
  
  % Convert to Matlab sparse matrix
  mesh.(name) = CSR_to_sparse( A);

end

end