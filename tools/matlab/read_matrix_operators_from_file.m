function mesh = read_matrix_operators_from_file( mesh, filename)

f = ncinfo( filename);

%% Read mesh translation tables and convert them to sparse matrices

% Read tables
mesh.n2vi     = ncread( filename, 'mesh_translation_tables/n2vi');
mesh.n2viuv   = ncread( filename, 'mesh_translation_tables/n2viuv');
mesh.n2vik    = ncread( filename, 'mesh_translation_tables/n2vi');
mesh.n2vikuv  = ncread( filename, 'mesh_translation_tables/n2viuv');
mesh.n2viks   = ncread( filename, 'mesh_translation_tables/n2vi');
mesh.n2viksuv = ncread( filename, 'mesh_translation_tables/n2viuv');

mesh.vi2n     = ncread( filename, 'mesh_translation_tables/vi2n');
mesh.viuv2n   = ncread( filename, 'mesh_translation_tables/viuv2n');
mesh.vik2n    = ncread( filename, 'mesh_translation_tables/vik2n');
mesh.vikuv2n  = ncread( filename, 'mesh_translation_tables/vikuv2n');
mesh.viks2n   = ncread( filename, 'mesh_translation_tables/viks2n');
mesh.viksuv2n = ncread( filename, 'mesh_translation_tables/viksuv2n');

mesh.n2ti     = ncread( filename, 'mesh_translation_tables/n2ti');
mesh.n2tiuv   = ncread( filename, 'mesh_translation_tables/n2tiuv');
mesh.n2tik    = ncread( filename, 'mesh_translation_tables/n2ti');
mesh.n2tikuv  = ncread( filename, 'mesh_translation_tables/n2tiuv');
mesh.n2tiks   = ncread( filename, 'mesh_translation_tables/n2ti');
mesh.n2tiksuv = ncread( filename, 'mesh_translation_tables/n2tiuv');

mesh.ti2n     = ncread( filename, 'mesh_translation_tables/ti2n');
mesh.tiuv2n   = ncread( filename, 'mesh_translation_tables/tiuv2n');
mesh.tik2n    = ncread( filename, 'mesh_translation_tables/tik2n');
mesh.tikuv2n  = ncread( filename, 'mesh_translation_tables/tikuv2n');
mesh.tiks2n   = ncread( filename, 'mesh_translation_tables/tiks2n');
mesh.tiksuv2n = ncread( filename, 'mesh_translation_tables/tiksuv2n');

mesh.n2ei     = ncread( filename, 'mesh_translation_tables/n2ei');
mesh.n2eiuv   = ncread( filename, 'mesh_translation_tables/n2eiuv');
mesh.n2eik    = ncread( filename, 'mesh_translation_tables/n2ei');
mesh.n2eikuv  = ncread( filename, 'mesh_translation_tables/n2eiuv');
mesh.n2eiks   = ncread( filename, 'mesh_translation_tables/n2ei');
mesh.n2eiksuv = ncread( filename, 'mesh_translation_tables/n2eiuv');

mesh.ei2n     = ncread( filename, 'mesh_translation_tables/ei2n');
mesh.eiuv2n   = ncread( filename, 'mesh_translation_tables/eiuv2n');
mesh.eik2n    = ncread( filename, 'mesh_translation_tables/eik2n');
mesh.eikuv2n  = ncread( filename, 'mesh_translation_tables/eikuv2n');
mesh.eiks2n   = ncread( filename, 'mesh_translation_tables/eiks2n');
mesh.eiksuv2n = ncread( filename, 'mesh_translation_tables/eiksuv2n');

% Create sparse matrices

% b-bu
jj = int32((1:mesh.nTri)');
ii = mesh.tiuv2n( jj,1);
vv = double(jj*0) + 1;
mesh.M_b_bu = sparse( ii,jj,vv,mesh.nTri*2,mesh.nTri);
mesh.M_bu_b = mesh.M_b_bu';

% b-bv
jj = int32((1:mesh.nTri)');
ii = mesh.tiuv2n( jj,2);
vv = double(jj*0) + 1;
mesh.M_b_bv = sparse( ii,jj,vv,mesh.nTri*2,mesh.nTri);
mesh.M_bv_b = mesh.M_b_bv';

%% Read matrix operators
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