function analyse_integrated_test( varargin)
% Analyse the results of the MISMIP_mod integrated test

disp('Analysing integrated test idealised/MISMIP_mod...')
disp('')
  
test_name = 'MISMIP_mod';
test_path = 'integrated_tests/idealised/MISMIP_mod';

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

ors = {'east','northeast','north','northwest','west','southwest','south','southeast'};

for ori = 1: length( ors)

  or = ors{ ori};

  % Read model output
  filename_spinup  = [foldername_test '/results_spinup_10km/transect_'  or '.nc'];
  filename_advance = [foldername_test '/results_advance_10km/transect_' or '.nc'];
  filename_retreat = [foldername_test '/results_retreat_10km/transect_' or '.nc'];

  time_spinup = ncread( filename_spinup,'time');
  rGL_spinup  = ncread( filename_spinup,'grounding_line_distance_from_start');

  time_advance = ncread( filename_advance,'time');
  rGL_advance  = ncread( filename_advance,'grounding_line_distance_from_start');

  time_retreat = ncread( filename_retreat,'time');
  rGL_retreat  = ncread( filename_retreat,'grounding_line_distance_from_start');
  
  abs_GL_hyst.(or) = abs(rGL_retreat(end) - rGL_spinup(end));

end

% Read stability info for the last spin-up simulation
filename = [foldername_test '/results_spinup_10km/scalar_output_ANT_00001.nc'];
nskip = 5; % Skip the first few values, as the model is still relaxing there
stab = read_stability_info( filename, nskip);

% Write to scoreboard file
write_to_scoreboard_file( abs_GL_hyst, stab);

  function write_to_scoreboard_file( abs_GL_hyst, stab)

    % Set up a scoreboard results structure
    single_run = initialise_single_test_run( test_name, test_path);
  
    % Add cost functions to results structure
    single_run = add_cost_function_to_single_run( single_run, ...
      'GL_hyst_east', 'abs( rGL_retreat(end) - rGL_spinup(end) )', abs_GL_hyst.east);
    single_run = add_cost_function_to_single_run( single_run, ...
      'GL_hyst_northeast', 'abs( rGL_retreat(end) - rGL_spinup(end) )', abs_GL_hyst.northeast);
    single_run = add_cost_function_to_single_run( single_run, ...
      'GL_hyst_north', 'abs( rGL_retreat(end) - rGL_spinup(end) )', abs_GL_hyst.north);
    single_run = add_cost_function_to_single_run( single_run, ...
      'GL_hyst_northwest', 'abs( rGL_retreat(end) - rGL_spinup(end) )', abs_GL_hyst.northwest);
    single_run = add_cost_function_to_single_run( single_run, ...
      'GL_hyst_west', 'abs( rGL_retreat(end) - rGL_spinup(end) )', abs_GL_hyst.west);
    single_run = add_cost_function_to_single_run( single_run, ...
      'GL_hyst_southwest', 'abs( rGL_retreat(end) - rGL_spinup(end) )', abs_GL_hyst.southwest);
    single_run = add_cost_function_to_single_run( single_run, ...
      'GL_hyst_south', 'abs( rGL_retreat(end) - rGL_spinup(end) )', abs_GL_hyst.south);
    single_run = add_cost_function_to_single_run( single_run, ...
      'GL_hyst_southeast', 'abs( rGL_retreat(end) - rGL_spinup(end) )', abs_GL_hyst.southeast);

    single_run = add_stability_info_cost_functions( single_run, stab);
    
    % Write to scoreboard file
    write_scoreboard_file( foldername_automated_testing, single_run);
  
  end
  
end