clc
clear all
close all

filename = '../../results_20230612_001/test_mesh_operators_3D_output.nc';

mesh = read_mesh_from_file( filename);
zeta = ncread( filename,'zeta');
nz = length( zeta);

%% Actual values

wa = 130;
ha = 130;
margins_hor = [5, 5, 5, 5, 5, 5, 5, 5];
margins_ver = [5, 5, 45, 5, 45, 5, 5];
H = setup_multipanel_figure( wa, ha, margins_hor, margins_ver);

for i = 1: size(H.Ax,1)
  for j = 1: size(H.Ax,2)
    set( H.Ax{ i,j},'xtick',[],'ytick',[]);
  end
end

fieldnames = {
  'df_dx_bk'
  'df_dy_bk'
  'df_dz_bk'
  'd2f_dx2_bk'
  'd2f_dxdy_bk'
  'd2f_dy2_bk'
  'd2f_dz2_bk'
  };

for fi = 1: length( fieldnames)
  
  fieldname = fieldnames{ fi};
  fieldname_ex = [fieldname '_ex'];
  
  d    = ncread( filename, fieldname);
  d_ex = ncread( filename, fieldname_ex);
  
  % Surface
  plot_mesh_data_b( H.Ax{ 1,fi}, mesh, d_ex( :,1))
  plot_mesh_data_b( H.Ax{ 2,fi}, mesh, d(    :,1))
  clim = [min( d_ex( :,1)), max( d_ex( :,1))];
  if clim(2) == clim(1)
    if clim(2)==0
      clim = [-1,1];
    else
      clim(2) = clim(2) + 0.15*clim(1);
      clim(1) = clim(1) - 0.15*clim(1);
    end
  end
  set( H.Ax{ 1,fi},'clim',clim);
  set( H.Ax{ 2,fi},'clim',clim);
  
  % Column
  plot_mesh_data_b( H.Ax{ 3,fi}, mesh, d_ex( :,round(nz/2)))
  plot_mesh_data_b( H.Ax{ 4,fi}, mesh, d(    :,round(nz/2)))
  clim = [min( d_ex( :,round(nz/2))), max( d_ex( :,round(nz/2)))];
  if clim(2) == clim(1)
    if clim(2)==0
      clim = [-1,1];
    else
      clim(2) = clim(2) + 0.15*clim(1);
      clim(1) = clim(1) - 0.15*clim(1);
    end
  end
  set( H.Ax{ 3,fi},'clim',clim);
  set( H.Ax{ 4,fi},'clim',clim);
  
  % Base
  plot_mesh_data_b( H.Ax{ 5,fi}, mesh, d_ex( :,nz))
  plot_mesh_data_b( H.Ax{ 6,fi}, mesh, d(    :,nz))
  clim = [min( d_ex( :,nz)), max( d_ex( :,nz))];
  if clim(2) == clim(1)
    if clim(2)==0
      clim = [-1,1];
    else
      clim(2) = clim(2) + 0.15*clim(1);
      clim(1) = clim(1) - 0.15*clim(1);
    end
  end
  set( H.Ax{ 5,fi},'clim',clim);
  set( H.Ax{ 6,fi},'clim',clim);
  
end

%% Errors

wa = 200;
ha = 200;
margins_hor = [5, 5, 5, 5, 5, 5, 5, 5];
margins_ver = [5, 5, 5, 5];
H = setup_multipanel_figure( wa, ha, margins_hor, margins_ver);

for i = 1: size(H.Ax,1)
  for j = 1: size(H.Ax,2)
    set( H.Ax{ i,j},'xtick',[],'ytick',[],'clim',[0,0.1]);
  end
end

fieldnames = {
  'df_dx_bk'
  'df_dy_bk'
  'df_dz_bk'
  'd2f_dx2_bk'
  'd2f_dxdy_bk'
  'd2f_dy2_bk'
  'd2f_dz2_bk'
  };

for fi = 1: length( fieldnames)
  
  fieldname = fieldnames{ fi};
  fieldname_ex = [fieldname '_ex'];
  
  d    = ncread( filename, fieldname);
  d_ex = ncread( filename, fieldname_ex);
  err  = abs( d - d_ex) ./ max( abs(d_ex(:)));
  
  % Surface
  plot_mesh_data_b( H.Ax{ 1,fi}, mesh, err( :,2))
  
  % Column
  plot_mesh_data_b( H.Ax{ 2,fi}, mesh, err( :,round(nz/2)))
  
  % Base
  plot_mesh_data_b( H.Ax{ 3,fi}, mesh, err( :,nz-1))
  
end

function plot_mesh_data_b( ax, mesh, d)
  patch('parent',ax,'vertices',mesh.V(1:mesh.nV,:),'faces',mesh.Tri(1:mesh.nTri,:),...
    'facecolor','flat','facevertexcdata',d,'edgecolor','none');
  set(ax,'xlim',[mesh.xmin,mesh.xmax],'ylim',[mesh.ymin,mesh.ymax]);
end