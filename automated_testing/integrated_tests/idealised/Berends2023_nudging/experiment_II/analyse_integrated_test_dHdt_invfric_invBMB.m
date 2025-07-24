function analyse_integrated_test_dHdt_invfric_invBMB( varargin)
% Analyse the results of the MISMIP_mod integrated test

disp('Analysing integrated test idealised/Berends2023_nudging/experiment_II...')
disp('')
  
test_name = 'experiment_II';
test_path = 'integrated_tests/idealised/Berends2023_nudging';

%%

% In the GitHub Workflow, provide the component test results folders as
% input; but retain the option of running without input (i.e. as a script)
% locally, with user-defined folders

input_args = varargin;
if isempty( input_args)
  % Assume this is a local run

  foldername_automated_testing = '/Users/Beren017/Documents/GitHub/UFEMISM2.0/automated_testing';
  addpath('/Users/Beren017/Documents/GitHub/UFEMISM2.0/tools/matlab/')
  do_figures = true;

elseif isscalar( input_args)
  % Assume this is a GitHub Workflow run

  foldername_automated_testing = varargin{1};
  addpath('tools/matlab/')
  do_figures = false;

else
  error('need either foldername_automated_testing, or nothing as input!')
end

foldername_test = [foldername_automated_testing '/' test_path '/' test_name];
addpath([foldername_automated_testing '/scoreboard/scripts'])

%%

% Target values
filename = [foldername_test '/results_05_exp_II_warm_ocean_retreat_10yr_5km/main_output_ANT_00001.nc'];
mesh = read_mesh_from_file( filename);
time = ncread( filename,'time'); ti = length( time);
phi_fric_target  = ncread( filename,'till_friction_angle',[1,ti],[Inf,1]);
Hi_target        = ncread( filename,'Hi'                 ,[1,ti],[Inf,1]);
Hb_target        = ncread( filename,'Hb'                 ,[1,ti],[Inf,1]);
Hs_target        = ncread( filename,'Hs'                 ,[1,ti],[Inf,1]);
uabs_surf_target = ncread( filename,'uabs_surf'          ,[1,ti],[Inf,1]);
R_shear_target   = ncread( filename,'R_shear'            ,[1,ti],[Inf,1]);
BMB_target       = ncread( filename,'BMB'                ,[1,ti],[Inf,1]);

% Inverted values
filename = [foldername_test '/results_06_exp_II_dHdt_invfric_invBMB_5km/main_output_ANT_00001.nc'];
time = ncread( filename,'time'); ti = length( time); %ti = find( time==2000);
phi_fric_inverted  = ncread( filename,'till_friction_angle',[1,ti],[Inf,1]);
Hs_inverted        = ncread( filename,'Hs'                 ,[1,ti],[Inf,1]);
uabs_surf_inverted = ncread( filename,'uabs_surf'          ,[1,ti],[Inf,1]);
BMB_inverted       = ncread( filename,'BMB'                ,[1,ti],[Inf,1]);

% Ice masks
seawater_density = 1028;
ice_density = 910;
TAF_target = Hi_target - max( 0, -Hb_target * (seawater_density / ice_density));

% Sliding grounded ice masks (R_shear > 0.05)
mask_a = TAF_target > 0 & R_shear_target > 0.05;
mask_b = false( mesh.nTri,1);
mask_shelf = TAF_target < 0 & Hi_target > 0.1;

for ti = 1: mesh.nTri
  has_icefree = false;
  for n = 1: 3
    vi = mesh.Tri( ti,n);
    if ~mask_a( vi)
      has_icefree = true;
      break
    end
  end
  if ~has_icefree
    mask_b( ti) = true;
  end
end

% Transparency maps for figures
a = zeros( mesh.nV,1) + 1;
a( ~mask_a) = 0.1;
b = zeros( mesh.nTri,1) + 1;
b( ~mask_b) = 0.1;
a_shelf = zeros( mesh.nV,1) + 1;
a_shelf( ~mask_shelf) = 0.1;

if (do_figures)

  % Till friction angle
  d = phi_fric_inverted ./ phi_fric_target;
  H = plot_mesh_data( mesh, d);
  set( H.Patch,'facevertexalphadata',a,'alphadatamapping','none','facealpha','flat');
  set( H.Ax,'clim',[0.8,1.25]);
  title( H.Ax,'Till friction angle')
  ylabel( H.Cbar,'inverted / target')
  colormap( H.Ax, flipud( crameri('roma')))
  set( H.Ax,'colorscale','log')

  % Ice geometry
  d = Hs_inverted - Hs_target;
  H = plot_mesh_data( mesh, d);
  % set( H.Patch,'facevertexalphadata',a,'alphadatamapping','none','facealpha','flat');
  set( H.Ax,'clim',[-50, 50]);
  title( H.Ax,'Ice thickness')
  ylabel( H.Cbar,'inverted - target')
  colormap( H.Ax, crameri('vik'))

  % Ice velocity
  d = (5 + uabs_surf_inverted) ./ (5 + uabs_surf_target);
  H = plot_mesh_data( mesh, d);
  set( H.Patch,'facevertexalphadata',b,'alphadatamapping','none','facealpha','flat');
  set( H.Ax,'clim',[0.8, 1.25]);
  title( H.Ax,'Ice velocity')
  ylabel( H.Cbar,'inverted / target')
  colormap( H.Ax, flipud( crameri('bam')))
  set( H.Ax,'colorscale','log')

  % Basal melt rate
  d = BMB_inverted - BMB_target;
  H = plot_mesh_data( mesh, d);
  set( H.Patch,'facevertexalphadata',a_shelf,'alphadatamapping','none','facealpha','flat');
  set( H.Ax,'clim',[-1, 1]);
  title( H.Ax,'Basal melt rate')
  ylabel( H.Cbar,'inverted - target')
  colormap( H.Ax, flipud( crameri('bam')))

end

% Calculate 95th percentiles
r95_till_friction_angle = calc_95th_percentile_ratio( phi_fric_target(  mask_a),     phi_fric_inverted(  mask_a));
p95_ice_thickness       = calc_95th_percentile(       Hs_target(        mask_a),     Hs_inverted(        mask_a));
r95_ice_velocity        = calc_95th_percentile_ratio( uabs_surf_target( mask_b)+5,   uabs_surf_inverted( mask_b)+5);
p95_melt_rate           = calc_95th_percentile(       BMB_target(       mask_shelf), BMB_inverted(       mask_shelf));

% Read stability info for the last spin-up simulation
filename = [foldername_test '/results_06_exp_II_dHdt_invfric_invBMB_5km/scalar_output_ANT_00001.nc'];
nskip = 5; % Skip the first few values, as the model is still relaxing there
stab = read_stability_info( filename, nskip);

% Write to scoreboard file
write_to_scoreboard_file( r95_till_friction_angle, ...
  p95_ice_thickness, r95_ice_velocity, p95_melt_rate, stab);

  function write_to_scoreboard_file( r95_till_friction_angle, ...
    p95_ice_thickness, r95_ice_velocity, p95_melt_rate, stab)

    % Set up a scoreboard results structure
    single_run = initialise_single_test_run( [test_name '_dHdt_invfric_invBMB'], test_path);
  
    % Add cost functions to results structure
    single_run = add_cost_function_to_single_run( single_run, ...
      'r95_till_friction_angle', '95% of till friction is within this fraction of its target', r95_till_friction_angle);
    single_run = add_cost_function_to_single_run( single_run, ...
      'p95_ice_thickness', '95% of ice thickness is within this range of its target [m]', p95_ice_thickness);
    single_run = add_cost_function_to_single_run( single_run, ...
      'r95_ice_velocity', '95% of ice velocity is within this fraction of its target', r95_ice_velocity);
    single_run = add_cost_function_to_single_run( single_run, ...
      'p95_melt_rate', '95% of melt rates are within this range of its target [m/yr]', p95_melt_rate);

    single_run = add_stability_info_cost_functions( single_run, stab);
    
    % Write to scoreboard file
    write_scoreboard_file( foldername_automated_testing, single_run);
  
  end
  
end