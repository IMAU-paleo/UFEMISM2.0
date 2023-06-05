clc
clear all
close all

filename = '../../results_20230605_001/test_remapping_lonlat2mesh_output.nc';

mesh = read_mesh_from_file( filename);

wa = 200;
ha = 200;
margins_hor = [25, 25, 25];
margins_ver = [25, 25];
H = setup_multipanel_figure( wa, ha, margins_hor, margins_ver);

for i = 1: 1
  for j = 1: 2
    set( H.Ax{ i,j},'xtick',[],'ytick',[],'clim',[-1,1]);
  end
end

d_mesh_ex = ncread( filename,'d_mesh_ex');
d_mesh    = ncread( filename,'d_mesh');

plot_mesh_data_a( H.Ax{ 1,1}, mesh, d_mesh_ex);
plot_mesh_data_a( H.Ax{ 1,2}, mesh, d_mesh   );

function plot_mesh_data_a( ax, mesh, d)
  patch('parent',ax,'vertices',mesh.V(1:mesh.nV,:),'faces',mesh.Tri(1:mesh.nTri,:),...
    'facecolor','interp','facevertexcdata',d,'edgecolor','none');
  set(ax,'xlim',[mesh.xmin,mesh.xmax],'ylim',[mesh.ymin,mesh.ymax]);
end