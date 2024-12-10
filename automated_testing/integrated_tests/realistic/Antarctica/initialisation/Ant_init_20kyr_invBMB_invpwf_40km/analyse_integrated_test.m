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

addpath([foldername_automated_testing '/scoreboard/scripts'])

%%

% Read model output
Hi = ncread( filename_results,'Hi');

RMSE_Hi = sqrt( mean( (Hi( :,end) - Hi( :,1)).^2 ));

write_to_scoreboard_file( RMSE_Hi);

  function write_to_scoreboard_file( RMSE_Hi)

    % Set up a scoreboard results structure
    single_run = initialise_single_test_run( test_name, test_path);
  
    % Add cost functions to results structure
    single_run = add_cost_function_to_single_run( single_run, ...
      'rmse', 'sqrt( mean( (Hi( :,end) - Hi( :,1)).^2 ))', RMSE_Hi);
    
    % Write to scoreboard file
    write_scoreboard_file( foldername_automated_testing, single_run);
  
  end
  
end