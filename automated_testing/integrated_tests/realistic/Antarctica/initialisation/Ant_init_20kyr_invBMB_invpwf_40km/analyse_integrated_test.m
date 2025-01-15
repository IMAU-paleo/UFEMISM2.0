function analyse_integrated_test( varargin)
% Analyse the results of the Antarctic initialisation run

disp('Analysing integrated test realistic/Antarctica/initialisation/Ant_init_20kyr_invBMB_invpwf_40km...')
disp('')
  
test_name = 'Ant_init_20kyr_invBMB_invpwf_40km';
test_path = 'integrated_tests/realistic/Antarctica/initialisation';

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

foldername_test = [foldername_automated_testing '/' test_path '/' test_name];

filename_results = [foldername_test '/results/main_output_ANT_00001.nc'];
filename_scalars = [foldername_test '/results/scalar_output_ANT_00001.nc'];

addpath([foldername_automated_testing '/scoreboard/scripts'])

%%

% Read model output
time       = ncread( filename_scalars,'time');
Hi         = ncread( filename_results,'Hi');
ice_volume = ncread( filename_scalars, 'ice_volume');

% Calculate cost functions
m_last5000yr = time > time(end) - 5000;
RMSE_Hi = sqrt( mean( (Hi( :,end) - Hi( :,1)).^2 ));
peak_volume_difference = max( abs( ice_volume - ice_volume(1)));
final_volume_difference = ice_volume( end) - ice_volume(1);
ice_volume_var = max( ice_volume( m_last5000yr)) - min( ice_volume( m_last5000yr));

% Read stability info
stab = read_stability_info( filename_scalars, 0);

% Write to file
write_to_scoreboard_file( RMSE_Hi, peak_volume_difference, ...
  final_volume_difference, ice_volume_var, stab);

  function write_to_scoreboard_file( RMSE_Hi, peak_volume_difference, ...
  final_volume_difference, ice_volume_var)

    % Set up a scoreboard results structure
    single_run = initialise_single_test_run( test_name, test_path);
  
    % Add cost functions to results structure
    single_run = add_cost_function_to_single_run( single_run, ...
      'rmse', 'sqrt( mean( (Hi( :,end) - Hi( :,1)).^2 ))', RMSE_Hi);
    single_run = add_cost_function_to_single_run( single_run, ...
      'peak_volume_difference', 'max( abs( ice_volume - ice_volume(1)))', peak_volume_difference);
    single_run = add_cost_function_to_single_run( single_run, ...
      'final_volume_difference', 'ice_volume( end) - ice_volume(1)', final_volume_difference);
    single_run = add_cost_function_to_single_run( single_run, ...
      'ice_volume_var', 'max( ice_volume( last 5000 yr)) - min( ice_volume( last 5000 yr))', ice_volume_var);

    single_run = add_stability_info_cost_functions( single_run, stab);
    
    % Write to scoreboard file
    write_scoreboard_file( foldername_automated_testing, single_run);
  
  end
  
end