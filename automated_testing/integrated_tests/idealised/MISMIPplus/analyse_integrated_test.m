function analyse_integrated_test( varargin)
% Analyse the results of the MISMIP_mod integrated test

disp('Analysing integrated test idealised/MISMIPplus...')
disp('')
  
test_name = 'MISMIPplus';
test_path = 'integrated_tests/idealised/MISMIPplus';

%%

% In the GitHub Workflow, provide the component test results folders as
% input; but retain the option of running without input (i.e. as a script)
% locally, with user-defined folders

input_args = varargin;
if isempty( input_args)
  % Assume this is a local run

  foldername_automated_testing = '/Users/Beren017/Documents/GitHub/UFEMISM2.0/automated_testing';
  addpath('/Users/Beren017/Documents/GitHub/UFEMISM2.0/tools/matlab/')

elseif isscalar( input_args)
  % Assume this is a GitHub Workflow run

  foldername_automated_testing = varargin{1};
  addpath('tools/matlab/')

else
  error('need either foldername_automated_testing, or nothing as input!')
end

foldername_test = [foldername_automated_testing '/' test_path];
addpath([foldername_automated_testing '/scoreboard/scripts'])

%%

% Read output of ice1r_4km
foldername = [foldername_test '/results_4km_ice1r'];
timeframes = get_UFEMISM_filelist( foldername, 'ANT');

time = zeros( length( timeframes),1);
x_GL = zeros( length( timeframes),1);

filename = '';
for tfi = 1: length( timeframes)

  % Read and interpolate UFEMISM output

  filename_prev = filename;
  filename      = timeframes(tfi).filename;
  ti            = timeframes(tfi).ti;

  if ~strcmp( filename_prev, filename)
    mesh = read_mesh_from_file( filename);
    xt = (mesh.xmin: 1000: mesh.xmax)';
    yt = xt * 0;
    A = calc_transect_matrix( mesh, xt, yt);
  end

  Hi = A * ncread( filename, 'Hi', [1,ti], [Inf,1]);
  Hb = A * ncread( filename, 'Hb', [1,ti], [Inf,1]);

  TAF = Hi - max( 0, (-Hb) * (1028 / 910));

  % Calculate GL position

  ir = 1;
  while TAF(ir) > 0
    ir = ir+1;
  end
  il = ir-1;

  x1 = xt(il);
  x2 = xt(ir);

  f0 = 0;
  f1 = TAF(il);
  f2 = TAF(ir);

  lambda = (f2 - f1) / (x2 - x1);
  x0 = x1 + (f0 - f1) / lambda;

  time( tfi) = timeframes( tfi).time;
  x_GL( tfi) = x0;
 
end

% Smooth grounding-line position over time
x_GL_smooth = x_GL;
for it = 1:7
  x_GL_smooth( 2:end-1) = ...
    0.25 * x_GL_smooth( 1:end-2) + ...
    0.5  * x_GL_smooth( 2:end-1) + ...
    0.25 * x_GL_smooth( 3:end  );
end

err_x_GL_init     = abs( x_GL(1) - 450e3);
err_x_GL_final_lo = abs( min( 0, x_GL( end) - 350e3));
err_x_GL_final_hi = abs( max( 0, x_GL( end) - 420e3));
var_x_GL          = max( abs( x_GL_smooth - x_GL));

% Read stability info for the last spin-up simulation
% (the retreat simulations have a very fast-changing ice sheet,
%  so the time step is pretty much always at the minimum)
filename = [foldername_test '/results_4km_spinup/scalar_output_ANT_00001.nc'];
nskip = 5; % Skip the first few values, as the model is still relaxing there
stab = read_stability_info( filename, nskip);

% Write to scoreboard file
write_to_scoreboard_file( ...
  err_x_GL_init, err_x_GL_final_lo, err_x_GL_final_hi, var_x_GL, stab);

  function write_to_scoreboard_file( ...
      err_x_GL_init, err_x_GL_final_lo, err_x_GL_final_hi, var_x_GL, stab)

    % Set up a scoreboard results structure
    single_run = initialise_single_test_run( test_name, test_path);
  
    % Add cost functions to results structure
    single_run = add_cost_function_to_single_run( single_run, ...
      'err_x_GL_init', 'abs( x_GL(1) - 450e3)', err_x_GL_init);
    single_run = add_cost_function_to_single_run( single_run, ...
      'err_x_GL_final_lo', 'abs( min( 0, x_GL( end) - 350e3))', err_x_GL_final_lo);
    single_run = add_cost_function_to_single_run( single_run, ...
      'err_x_GL_final_hi', 'abs( max( 0, x_GL( end) - 420e3))', err_x_GL_final_hi);
    single_run = add_cost_function_to_single_run( single_run, ...
      'var_x_GL', 'max( abs( x_GL_smooth - x_GL))', var_x_GL);

    single_run = add_stability_info_cost_functions( single_run, stab);
    
    % Write to scoreboard file
    write_scoreboard_file( foldername_automated_testing, single_run);
  
  end
  
end