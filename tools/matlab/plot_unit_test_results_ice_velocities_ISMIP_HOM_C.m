clc
clear all
close all

foldername = '../../results_20230619_001';

Ls    = [160,80,40,20,10,5];
umins = [0,0,10,13,13.5,6];
umaxs = [150,70,35,20,17,18];

for li = 1: length( Ls)
  
  L    = Ls(    li);
  umin = umins( li);
  umax = umaxs( li);

  filename = [foldername '/test_ISMIP_HOM_C_' num2str( L) 'km_output.nc'];

  mesh = read_mesh_from_file( filename);

  %% Plot

  wa = 300;
  ha = 300;
  margins_hor = [25, 25, 25, 100];
  margins_ver = [50, 25];
  H = setup_multipanel_figure( wa, ha, margins_hor, margins_ver);

  for i = 1: size( H.Ax,1)
    for j = 1: size( H.Ax,2)
      set( H.Ax{ 1,j},'xtick',[],'ytick',[],'clim',[umin,umax],'fontsize',24);
    end
  end

  title( H.Ax{ 1,1},'SIA/SSA');
  title( H.Ax{ 1,2},'DIVA');
  title( H.Ax{ 1,3},'BPA');
  pos = get( H.Ax{ 1,3},'position');
  H.Cbar = colorbar( H.Ax{ 1,3},'location','eastoutside');
  ylabel( H.Cbar,'Surface velocity (m/yr)')
  set( H.Ax{ 1,3},'position',pos);

  % SIA/SSA
  u_3D_b_SIASSA = ncread( filename,'u_3D_b_SIASSA');
  plot_mesh_data_b( H.Ax{ 1,1}, mesh, u_3D_b_SIASSA( :,1));

  % DIVA
  u_3D_b_DIVA = ncread( filename,'u_3D_b_DIVA');
  plot_mesh_data_b( H.Ax{ 1,2}, mesh, u_3D_b_DIVA( :,1));

  % BPA
  u_3D_b_BPA = ncread( filename,'u_3D_b_BPA');
  plot_mesh_data_b( H.Ax{ 1,3}, mesh, u_3D_b_BPA( :,1));

end

function plot_mesh_data_b( ax, mesh, d)
  patch('parent',ax,'vertices',mesh.V(1:mesh.nV,:),'faces',mesh.Tri(1:mesh.nTri,:),...
    'facecolor','flat','facevertexcdata',d,'edgecolor','none');
  set(ax,'xlim',[mesh.xmin,mesh.xmax],'ylim',[mesh.ymin,mesh.ymax]);
end