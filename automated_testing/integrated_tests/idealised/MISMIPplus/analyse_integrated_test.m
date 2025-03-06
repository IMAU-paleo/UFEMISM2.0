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
filename = [foldername_test '/results_4km_ice1r/transect_westeast.nc'];
x_GL = ncread( filename,'grounding_line_distance_from_start');

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