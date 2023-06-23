clc
clear all
close all

filename = '../../results_unit_tests/test_thickness_evolution_Halfar_dome_explicit_output.nc';

mesh = read_mesh_from_file( filename);

%% Plot

wa = 300;
ha = 300;
margins_hor = [25, 25, 25, 100];
margins_ver = [50, 25];
H = setup_multipanel_figure( wa, ha, margins_hor, margins_ver);

for i = 1: size( H.Ax,1)
  for j = 1: size( H.Ax,2)
    set( H.Ax{ 1,j},'xtick',[],'ytick',[],'clim',[-50,50],'fontsize',24);
  end
end

title( H.Ax{ 1,1},'Explicit');
title( H.Ax{ 1,2},'Implicit');
title( H.Ax{ 1,3},'Semi-implicit');

pos = get( H.Ax{ 1,3},'position');
H.Cbar = colorbar( H.Ax{ 1,3},'location','eastoutside');
ylabel( H.Cbar,'Ice thickness error (m)')
set( H.Ax{ 1,3},'position',pos);

% Explicit
filename = '../../results_unit_tests/test_thickness_evolution_Halfar_dome_explicit_output.nc';
Hi_modelled   = ncread( filename,'Hi_modelled');
Hi_analytical = ncread( filename,'Hi_analytical');
dHi = Hi_modelled - Hi_analytical;
plot_mesh_data_b( H.Ax{ 1,1}, mesh, dHi);

% Implicit
filename = '../../results_unit_tests/test_thickness_evolution_Halfar_dome_implicit_output.nc';
Hi_modelled   = ncread( filename,'Hi_modelled');
Hi_analytical = ncread( filename,'Hi_analytical');
dHi = Hi_modelled - Hi_analytical;
plot_mesh_data_b( H.Ax{ 1,2}, mesh, dHi);

% Semi-implicit
filename = '../../results_unit_tests/test_thickness_evolution_Halfar_dome_semi-implicit_output.nc';
Hi_modelled   = ncread( filename,'Hi_modelled');
Hi_analytical = ncread( filename,'Hi_analytical');
dHi = Hi_modelled - Hi_analytical;
plot_mesh_data_b( H.Ax{ 1,3}, mesh, dHi);

function plot_mesh_data_b( ax, mesh, d)
  patch('parent',ax,'vertices',mesh.V(1:mesh.nV,:),'faces',mesh.Tri(1:mesh.nTri,:),...
    'facecolor','flat','facevertexcdata',d,'edgecolor','none');
  set(ax,'xlim',[mesh.xmin,mesh.xmax],'ylim',[mesh.ymin,mesh.ymax]);
end