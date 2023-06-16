clc
clear all
close all

filename = '../../results_20230616_001/test_ice_velocities_Halfar_dome_output.nc';

mesh = read_mesh_from_file( filename);

%% Plot

wa = 300;
ha = 300;
margins_hor = [90, 25, 25, 100];
margins_ver = [50, 25, 80];
H = setup_multipanel_figure( wa, ha, margins_hor, margins_ver);

for i = 1: size( H.Ax,1)
  for j = 1: size( H.Ax,2)
    set( H.Ax{ 1,j},'xtick',[],'ytick',[],'clim',[0.6667,150],'colorscale','log','fontsize',24);
    set( H.Ax{ 2,j},'xlim',[0,550],'ylim',[0,250],'fontsize',24);
  end
end

set( H.Ax{ 2,2},'yticklabels',[]);
set( H.Ax{ 2,3},'yticklabels',[]);

xlabel( H.Ax{ 2,1},'Distance from divide (km)')
xlabel( H.Ax{ 2,2},'Distance from divide (km)')
xlabel( H.Ax{ 2,3},'Distance from divide (km)')
ylabel( H.Ax{ 2,1},'Surface velocity (m/yr)')

title( H.Ax{ 1,1},'SIA');
title( H.Ax{ 1,2},'DIVA');
title( H.Ax{ 1,3},'BPA');

pos = get( H.Ax{ 1,3},'position');
H.Cbar = colorbar( H.Ax{ 1,3},'location','eastoutside');
ylabel( H.Cbar,'Surface velocity (m/yr)')
set( H.Ax{ 1,3},'position',pos);

R = sqrt( mesh.TriGC( :,1).^2 + mesh.TriGC( :,2).^2) / 1e3; % in km

% SIA
u_3D_b_SIA = ncread( filename,'u_3D_b_SIA');
v_3D_b_SIA = ncread( filename,'v_3D_b_SIA');
uabs_surf_SIA = sqrt( u_3D_b_SIA( :,1).^2 + v_3D_b_SIA( :,1).^2);
plot_mesh_data_b( H.Ax{ 1,1}, mesh, uabs_surf_SIA);
line('parent',H.Ax{ 2,1},'xdata',R,'ydata',uabs_surf_SIA,'linestyle','none','marker','.')

% DIVA
u_3D_b_DIVA = ncread( filename,'u_3D_b_DIVA');
v_3D_b_DIVA = ncread( filename,'v_3D_b_DIVA');
uabs_surf_DIVA = sqrt( u_3D_b_DIVA( :,1).^2 + v_3D_b_DIVA( :,1).^2);
plot_mesh_data_b( H.Ax{ 1,2}, mesh, uabs_surf_DIVA);
line('parent',H.Ax{ 2,2},'xdata',R,'ydata',uabs_surf_DIVA,'linestyle','none','marker','.')

% BPA
u_3D_b_BPA = ncread( filename,'u_3D_b_BPA');
v_3D_b_BPA = ncread( filename,'v_3D_b_BPA');
uabs_surf_BPA = sqrt( u_3D_b_BPA( :,1).^2 + v_3D_b_BPA( :,1).^2);
plot_mesh_data_b( H.Ax{ 1,3}, mesh, uabs_surf_BPA);
line('parent',H.Ax{ 2,3},'xdata',R,'ydata',uabs_surf_BPA,'linestyle','none','marker','.')

function plot_mesh_data_b( ax, mesh, d)
  patch('parent',ax,'vertices',mesh.V(1:mesh.nV,:),'faces',mesh.Tri(1:mesh.nTri,:),...
    'facecolor','flat','facevertexcdata',d,'edgecolor','none');
  set(ax,'xlim',[mesh.xmin,mesh.xmax],'ylim',[mesh.ymin,mesh.ymax]);
end