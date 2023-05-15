clc
clear all
close all

foldername = '../../results_20230424_001';

%% Read stuff

mesh = read_mesh_from_file( [foldername '/testfile.nc']);

% M_map_a_a = read_CSR_from_NetCDF( [foldername '/M_map_a_a.nc']);
M_ddx_a_a = read_CSR_from_NetCDF( [foldername '/M_ddx_a_a.nc']);
M_ddy_a_a = read_CSR_from_NetCDF( [foldername '/M_ddy_a_a.nc']);

M_map_a_b = read_CSR_from_NetCDF( [foldername '/M_map_a_b.nc']);
M_ddx_a_b = read_CSR_from_NetCDF( [foldername '/M_ddx_a_b.nc']);
M_ddy_a_b = read_CSR_from_NetCDF( [foldername '/M_ddy_a_b.nc']);

M_map_a_c = read_CSR_from_NetCDF( [foldername '/M_map_a_c.nc']);
M_ddx_a_c = read_CSR_from_NetCDF( [foldername '/M_ddx_a_c.nc']);
M_ddy_a_c = read_CSR_from_NetCDF( [foldername '/M_ddy_a_c.nc']);


M_map_b_a = read_CSR_from_NetCDF( [foldername '/M_map_b_a.nc']);
M_ddx_b_a = read_CSR_from_NetCDF( [foldername '/M_ddx_b_a.nc']);
M_ddy_b_a = read_CSR_from_NetCDF( [foldername '/M_ddy_b_a.nc']);

% M_map_b_b = read_CSR_from_NetCDF( [foldername '/M_map_b_b.nc']);
M_ddx_b_b = read_CSR_from_NetCDF( [foldername '/M_ddx_b_b.nc']);
M_ddy_b_b = read_CSR_from_NetCDF( [foldername '/M_ddy_b_b.nc']);

M_map_b_c = read_CSR_from_NetCDF( [foldername '/M_map_b_c.nc']);
M_ddx_b_c = read_CSR_from_NetCDF( [foldername '/M_ddx_b_c.nc']);
M_ddy_b_c = read_CSR_from_NetCDF( [foldername '/M_ddy_b_c.nc']);


M_map_c_a = read_CSR_from_NetCDF( [foldername '/M_map_c_a.nc']);
M_ddx_c_a = read_CSR_from_NetCDF( [foldername '/M_ddx_c_a.nc']);
M_ddy_c_a = read_CSR_from_NetCDF( [foldername '/M_ddy_c_a.nc']);

M_map_c_b = read_CSR_from_NetCDF( [foldername '/M_map_c_b.nc']);
M_ddx_c_b = read_CSR_from_NetCDF( [foldername '/M_ddx_c_b.nc']);
M_ddy_c_b = read_CSR_from_NetCDF( [foldername '/M_ddy_c_b.nc']);

% M_map_c_c = read_CSR_from_NetCDF( [foldername '/M_map_c_c.nc']);
M_ddx_c_c = read_CSR_from_NetCDF( [foldername '/M_ddx_c_c.nc']);
M_ddy_c_c = read_CSR_from_NetCDF( [foldername '/M_ddy_c_c.nc']);


M2_ddx_b_b    = read_CSR_from_NetCDF( [foldername '/M2_ddx_b_b.nc']);
M2_ddy_b_b    = read_CSR_from_NetCDF( [foldername '/M2_ddy_b_b.nc']);
M2_d2dx2_b_b  = read_CSR_from_NetCDF( [foldername '/M2_d2dx2_b_b.nc']);
M2_d2dxdy_b_b = read_CSR_from_NetCDF( [foldername '/M2_d2dxdy_b_b.nc']);
M2_d2dy2_b_b  = read_CSR_from_NetCDF( [foldername '/M2_d2dy2_b_b.nc']);

%% Solve Poisson equation

% Calculate N
N_a = zeros( mesh.nV, 1);
amp = 7.0;
for vi = 1: mesh.nV
  x = mesh.V( vi,1);
  y = mesh.V( vi,2);
  xp = (x - mesh.xmin) / (mesh.xmax - mesh.xmin);
  yp = (y - mesh.ymin) / (mesh.ymax - mesh.ymin);
  N_a( vi) = 1.0 + exp(amp * sin( 2.0 * pi * xp) * sin( 2.0 * pi * yp));
end

N_b     = M_map_a_b * N_a;
dN_dx_b = M_ddx_a_b * N_a;
dN_dy_b = M_ddy_a_b * N_a;

% Calculate stiffness matrix A
A_b = M_ddx_a_b * diag_sparse( N_a) * M_ddx_b_a + M_ddy_a_b * diag_sparse( N_a) * M_ddy_b_a;
A_c = M_ddx_a_c * diag_sparse( N_a) * M_ddx_c_a + M_ddy_a_c * diag_sparse( N_a) * M_ddy_c_a;

A2_b = diag_sparse( N_b) * M2_d2dx2_b_b + diag_sparse( dN_dx_b) * M2_ddx_b_b + ...
       diag_sparse( N_b) * M2_d2dy2_b_b + diag_sparse( dN_dy_b) * M2_ddy_b_b;

% Add boundary conditions
for ti = 1: mesh.nTri
  if (mesh.TriBI( ti) > 0)
    A_b( ti,:) = 0;
    A_b( ti,ti) = 1;
    A2_b( ti,:) = 0;
    A2_b( ti,ti) = 1;
  end
end
for ei = 1: mesh.nE
  if (mesh.EBI( ei) > 0)
    A_c( ei,:) = 0;
    A_c( ei,ei) = 1;
  end
end

% Calculate load vector b
b_b = zeros( mesh.nTri,1);
for ti = 1: mesh.nTri
  if (mesh.TriBI( ti) == 0)
    % Free triangle

    b_b( ti) = 0.0;

  else
    % Border triangle

    if (mesh.TriBI( ti) == 6 || mesh.TriBI( ti) == 7 || mesh.TriBI( ti) == 8)
      % Western border: bump
      y = mesh.TriGC( ti,2);
      y = (y - mesh.ymin) / (mesh.ymax - mesh.ymin);
      b_b( ti) = 0.5 * (1.0 - cos( 2.0 * pi * y));
    else
      % Other borders: zero
      b_b( ti) = 0.0;
    end

  end

end

b_c = zeros( mesh.nE,1);
for ei = 1: mesh.nE
  if (mesh.EBI( ei) == 0)
    % Free edge

    b_c( ei) = 0.0;

  else
    % Border edge

    if (mesh.EBI( ei) == 6 || mesh.EBI( ei) == 7 || mesh.EBI( ei) == 8)
      % Western border: bump
      y = mesh.E( ei,2);
      y = (y - mesh.ymin) / (mesh.ymax - mesh.ymin);
      b_c( ei) = 0.5 * (1.0 - cos( 2.0 * pi * y));
    else
      % Other borders: zero
      b_c( ei) = 0.0;
    end

  end

end

% Solve
x_b  = A_b  \ b_b;
x_c  = A_c  \ b_c;
x2_b = A2_b \ b_b;

% Plot
plot_mesh_data_b( mesh, x_b);
plot_mesh_data_c( mesh, x_c);
plot_mesh_data_b( mesh, x2_b);

function A = diag_sparse( A_vec)
A = sparse( 1:length(A_vec), 1:length(A_vec), A_vec, length(A_vec), length(A_vec));
end