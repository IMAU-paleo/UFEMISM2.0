clc
clear all
close all

filename1 = '../../results_20230605_001/test_remapping_mesh2mesh_output1.nc';
filename2 = '../../results_20230605_001/test_remapping_mesh2mesh_output2.nc';

mesh1 = read_mesh_from_file( filename1);
mesh2 = read_mesh_from_file( filename2);

wa = 200;
ha = 200;
margins_hor = [25, 25, 25, 25, 25];
margins_ver = [25, 25, 25];
H = setup_multipanel_figure( wa, ha, margins_hor, margins_ver);

for i = 1: 1
  for j = 1: 2
    set( H.Ax{ i,j},'xtick',[],'ytick',[],'clim',[-1,1]);
  end
end

d1_ex      = ncread( filename1,'d1_ex');
d21_nn     = ncread( filename1,'d21_nn');
d21_trilin = ncread( filename1,'d21_trilin');
d21_cons   = ncread( filename1,'d21_cons');

plot_mesh_data_a( H.Ax{ 1,1}, mesh1, d1_ex);
plot_mesh_data_a( H.Ax{ 1,2}, mesh1, d21_nn);
plot_mesh_data_a( H.Ax{ 1,3}, mesh1, d21_trilin);
plot_mesh_data_a( H.Ax{ 1,4}, mesh1, d21_cons);

d2_ex      = ncread( filename2,'d2_ex');
d12_nn     = ncread( filename2,'d12_nn');
d12_trilin = ncread( filename2,'d12_trilin');
d12_cons   = ncread( filename2,'d12_cons');

plot_mesh_data_a( H.Ax{ 2,1}, mesh2, d2_ex);
plot_mesh_data_a( H.Ax{ 2,2}, mesh2, d12_nn);
plot_mesh_data_a( H.Ax{ 2,3}, mesh2, d12_trilin);
plot_mesh_data_a( H.Ax{ 2,4}, mesh2, d12_cons);

function plot_mesh_data_a( ax, mesh, d)
  patch('parent',ax,'vertices',mesh.V(1:mesh.nV,:),'faces',mesh.Tri(1:mesh.nTri,:),...
    'facecolor','interp','facevertexcdata',d,'edgecolor','none');
  set(ax,'xlim',[mesh.xmin,mesh.xmax],'ylim',[mesh.ymin,mesh.ymax]);
end