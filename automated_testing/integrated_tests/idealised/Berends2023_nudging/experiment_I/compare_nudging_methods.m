clc
clear all
close all

cmap_absolute = flipud( crameri( 'lajolla'));
cmap_error    = crameri( 'vik');

clim_phi_fric     = [0.8,2.0] .* [0.9,1/0.9];
clim_phi_fric_err = [0.5,2.0];

clim_Hs           = [0,3000];
clim_Hs_err       = [-50,50];

clim_u            = [0.1,100];
clim_u_err        = [0.5,2.0];

%% Read mesh and target quantities
filename_target = 'results_04_exp_I_spinup_5km/main_output_ANT_00001.nc';
mesh = read_mesh_from_file( filename_target);
time = ncread( filename_target,'time'); ti = length( time);
target.Hi        = ncread( filename_target,'Hi'                 ,[1,ti],[Inf,1]);
target.Hs        = ncread( filename_target,'Hs'                 ,[1,ti],[Inf,1]);
target.phi_fric  = ncread( filename_target,'till_friction_angle',[1,ti],[Inf,1]);
target.uabs_surf = ncread( filename_target,'uabs_surf'          ,[1,ti],[Inf,1]);

mask_a = target.Hi > 0.1;
mask_b = false( mesh.nTri,1);
for ti = 1: mesh.nTri
  all_are_true = true;
  for n = 1: 3
    vi = mesh.Tri( ti,n);
    if ~mask_a( vi); all_are_true = false; end
  end
  if all_are_true; mask_b( ti) = true; end
end

%% Read model results

model_names = {
  'results_05_exp_I_inversion_5km_H_dHdt_flowline'
  'results_06_exp_I_inversion_5km_H_dHdt_local'
  'results_07_exp_I_inversion_5km_H_u_flowline'
  };

for mi = 1: length( model_names)
  filename_model = [model_names{ mi} '/main_output_ANT_00001.nc'];
  time = ncread( filename_model,'time'); ti = length( time);
  model( mi).Hs        = ncread( filename_model,'Hs'                 ,[1,ti],[Inf,1]);
  model( mi).phi_fric  = ncread( filename_model,'till_friction_angle',[1,ti],[Inf,1]);
  model( mi).uabs_surf = ncread( filename_model,'uabs_surf'          ,[1,ti],[Inf,1]);
  model( mi).Hs_err        = model( mi).Hs - target.Hs;
  model( mi).phi_fric_err  = model( mi).phi_fric ./ target.phi_fric;
  model( mi).uabs_surf_err = model( mi).uabs_surf ./ target.uabs_surf;
end

%% Setup figure
ha = 120;
wa = ha * (mesh.xmax - mesh.xmin) / (mesh.ymax - mesh.ymin);
H = setup_multipanel_figure( wa, ha, [100,35,0,35,0,35,0,100],[25,25,25,25]);

for i = 1: size( H.Ax,1)
  for j = 1: size( H.Ax,2)

    ax = H.Ax{ i,j};
    set( ax, 'xlim',[mesh.xmin,mesh.xmax]*0.8,'ylim',[mesh.ymin,mesh.ymax]*0.8,...
      'xtick',[],'ytick',[],'xgrid','off','ygrid','off');
    H.Ax{ i,j}.XAxis.Visible = 'off';
    H.Ax{ i,j}.YAxis.Visible = 'off';

    if i < 3
      H.Patch{ i,j} = patch( 'parent',ax,'vertices',mesh.Vor,'faces',...
        changem( double( mesh.VVor),NaN),'facecolor','flat',...
        'edgecolor','none','facevertexalphadata',double( mask_a),...
        'alphadatamapping','none','facealpha','flat');
    else
      H.Patch{ i,j} = patch( 'parent',ax,'vertices',mesh.V,'faces',...
        mesh.Tri,'facecolor','flat',...
        'edgecolor','none','facevertexalphadata',double( mask_b),...
        'alphadatamapping','none','facealpha','flat');
    end

    drawnow('update')

  end
end

%% Till friction angle

% Colorbar for absolute values
ax = H.Ax{ 1,1};
pos = get( ax,'position');
H.CBar_phi_fric = colorbar( ax,'location','westoutside');
set( ax,'position',pos);
ylabel( H.CBar_phi_fric,'\phi [deg]')

% Colorbar for errors
ax = H.Ax{ 1,end};
pos = get( ax,'position');
H.CBar_phi_fric_err = colorbar( ax,'location','eastoutside');
set( ax,'position',pos);
ylabel( H.CBar_phi_fric_err,'r_{\phi} [frac]')

% Plot target
colormap( H.Ax{ 1,1}, cmap_absolute)
set( H.Ax{ 1,1},'clim',clim_phi_fric,'colorscale','log')
set( H.Patch{ 1,1},'facevertexcdata',target.phi_fric)
drawnow('update')

% Plot models
for mi = 1: length( model)

  ax_absolute    = H.Ax{    1,2*mi};
  patch_absolute = H.Patch{ 1,2*mi};
  ax_error       = H.Ax{    1,2*mi+1};
  patch_error    = H.Patch{ 1,2*mi+1};

  colormap( ax_absolute, cmap_absolute)
  set( ax_absolute,'clim',clim_phi_fric,'colorscale','log')
  set( patch_absolute,'facevertexcdata',model( mi).phi_fric)

  colormap( ax_error, cmap_error)
  set( ax_error,'clim',clim_phi_fric_err,'colorscale','log')
  set( patch_error,'facevertexcdata',model( mi).phi_fric_err)

  drawnow('update')
end

%% Surface elevation

% Colorbar for absolute values
ax = H.Ax{ 2,1};
pos = get( ax,'position');
H.CBar_Hs = colorbar( ax,'location','westoutside');
set( ax,'position',pos);
ylabel( H.CBar_Hs,'s [m]')

% Colorbar for errors
ax = H.Ax{ 2,end};
pos = get( ax,'position');
H.CBar_Hs_err = colorbar( ax,'location','eastoutside');
set( ax,'position',pos);
ylabel( H.CBar_Hs_err,'\Deltas [m]')

% Plot target
colormap( H.Ax{ 2,1}, cmap_absolute)
set( H.Ax{ 2,1},'clim',clim_Hs)
set( H.Patch{ 2,1},'facevertexcdata',target.Hs)
drawnow('update')

% Plot models
for mi = 1: length( model)

  ax_absolute    = H.Ax{    2,2*mi};
  patch_absolute = H.Patch{ 2,2*mi};
  ax_error       = H.Ax{    2,2*mi+1};
  patch_error    = H.Patch{ 2,2*mi+1};

  colormap( ax_absolute, cmap_absolute)
  set( ax_absolute,'clim',clim_Hs)
  set( patch_absolute,'facevertexcdata',model( mi).Hs)

  colormap( ax_error, cmap_error)
  set( ax_error,'clim',clim_Hs_err)
  set( patch_error,'facevertexcdata',model( mi).Hs_err)
  drawnow('update')
end

%% Velocity

% Colorbar for absolute values
ax = H.Ax{ 3,1};
pos = get( ax,'position');
H.CBar_u = colorbar( ax,'location','westoutside');
set( ax,'position',pos);
ylabel( H.CBar_u,'u [m/yr]')

% Colorbar for errors
ax = H.Ax{ 3,end};
pos = get( ax,'position');
H.CBar_u_err = colorbar( ax,'location','eastoutside');
set( ax,'position',pos);
ylabel( H.CBar_u_err,'r_u [frac]')

% Plot target
colormap( H.Ax{ 3,1}, cmap_absolute)
set( H.Ax{ 3,1},'clim',clim_u,'colorscale','log')
set( H.Patch{ 3,1},'facevertexcdata',target.uabs_surf)
drawnow('update')

% Plot models
for mi = 1: length( model)

  ax_absolute    = H.Ax{    3,2*mi};
  patch_absolute = H.Patch{ 3,2*mi};
  ax_error       = H.Ax{    3,2*mi+1};
  patch_error    = H.Patch{ 3,2*mi+1};

  colormap( ax_absolute, cmap_absolute)
  set( ax_absolute,'clim',clim_u,'colorscale','log')
  set( patch_absolute,'facevertexcdata',model( mi).uabs_surf)

  colormap( ax_error, cmap_error)
  set( ax_error,'clim',clim_u_err,'colorscale','log')
  set( patch_error,'facevertexcdata',model( mi).uabs_surf_err)
  drawnow('update')
end