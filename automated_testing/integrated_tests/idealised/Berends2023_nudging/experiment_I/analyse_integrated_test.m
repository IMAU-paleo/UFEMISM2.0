function analyse_integrated_test( varargin)
% Analyse the results of the MISMIP_mod integrated test

disp('Analysing integrated test idealised/Berends2023_nudging/experiment_I...')
disp('')
  
test_name = 'experiment_I';
test_path = 'integrated_tests/idealised/Berends2023_nudging/experiment_I';

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

foldername_test = [foldername_automated_testing '/' test_path];
addpath([foldername_automated_testing '/scoreboard/scripts'])

%%

% Target values
filename = [foldername_test '/results_04_exp_I_spinup_5km/main_output_ANT_00001.nc'];
mesh = read_mesh_from_file( filename);
time = ncread( filename,'time'); ti = length( time);
phi_fric_target  = ncread( filename,'till_friction_angle',[1,ti],[Inf,1]);
Hs_target        = ncread( filename,'Hs'                 ,[1,ti],[Inf,1]);
uabs_surf_target = ncread( filename,'uabs_surf'          ,[1,ti],[Inf,1]);

% Inverted values
filename = [foldername_test '/results_05_exp_I_inversion_5km/main_output_ANT_00001.nc'];
time = ncread( filename,'time'); ti = length( time);
phi_fric_inverted  = ncread( filename,'till_friction_angle',[1,ti],[Inf,1]);
Hs_inverted        = ncread( filename,'Hs'                 ,[1,ti],[Inf,1]);
uabs_surf_inverted = ncread( filename,'uabs_surf'          ,[1,ti],[Inf,1]);

% Ice masks
mask_a = Hs_target > 2;
mask_b = false( mesh.nTri,1);
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

if (do_figures)

  % Till friction angle
  d = phi_fric_inverted ./ phi_fric_target;
  d( ~mask_a) = NaN;
  H = plot_mesh_data( mesh, d);
  set( H.Ax,'clim',[0.8,1.25]);
  title( H.Ax,'Till friction angle')
  ylabel( H.Cbar,'inverted / target')
  colormap( H.Ax, crameri('roma'))
  set( H.Ax,'colorscale','log')

  % Ice geometry
  d = Hs_inverted - Hs_target;
  d( ~mask_a) = NaN;
  H = plot_mesh_data( mesh, d);
  set( H.Ax,'clim',[-100, 100]);
  title( H.Ax,'Ice thickness')
  ylabel( H.Cbar,'inverted - target')
  colormap( H.Ax, crameri('vik'))

  % Ice velocity
  d = (0.5 + uabs_surf_inverted) ./ (0.5 + uabs_surf_target);
  d( ~mask_b) = NaN;
  H = plot_mesh_data( mesh, d);
  set( H.Ax,'clim',[0.8, 1.25]);
  title( H.Ax,'Ice velocity')
  ylabel( H.Cbar,'inverted / target')
  colormap( H.Ax, crameri('bam'))
  set( H.Ax,'colorscale','log')

end

% Calculate 95th percentiles
r95_till_friction_angle = calc_95th_percentile_ratio( phi_fric_target(  mask_a), phi_fric_inverted(  mask_a));
p95_ice_thickness       = calc_95th_percentile(       Hs_target(        mask_a), Hs_inverted(        mask_a));
r95_ice_velocity        = calc_95th_percentile_ratio( uabs_surf_target( mask_b), uabs_surf_inverted( mask_b));

% Read stability info for the last spin-up simulation
filename = [foldername_test '/results_05_exp_I_inversion_5km/scalar_output_ANT_00001.nc'];
nskip = 5; % Skip the first few values, as the model is still relaxing there
stab = read_stability_info( filename, nskip);

% Write to scoreboard file
write_to_scoreboard_file( r95_till_friction_angle, ...
  p95_ice_thickness, r95_ice_velocity, stab);

  function write_to_scoreboard_file( r95_till_friction_angle, ...
    p95_ice_thickness, r95_ice_velocity, stab)

    % Set up a scoreboard results structure
    single_run = initialise_single_test_run( test_name, test_path);
  
    % Add cost functions to results structure
    single_run = add_cost_function_to_single_run( single_run, ...
      'r95_till_friction_angle', '95% of till friction is within this fraction of its target', r95_till_friction_angle);
    single_run = add_cost_function_to_single_run( single_run, ...
      'p95_ice_thickness', '95% of ice thickness is within this range of its target', p95_ice_thickness);
    single_run = add_cost_function_to_single_run( single_run, ...
      'r95_ice_velocity', '95% of ice velocity is within this fraction of its target', r95_ice_velocity);

    single_run = add_stability_info_cost_functions( single_run, stab);
    
    % Write to scoreboard file
    write_scoreboard_file( foldername_automated_testing, single_run);
  
  end
  
end