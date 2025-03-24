function analyse_integrated_test( varargin)
% Analyse the results of the SSA_icestream integrated test

disp('Analysing integrated test idealised/SSA_icestream...')
disp('')
  
test_name = 'SSA_icestream';
test_path = 'integrated_tests/idealised/SSA_icestream';

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

% Read and process UFEMISM output
foldernames = {...
  [foldername_test '/results_32km']
  [foldername_test '/results_16km']
  [foldername_test '/results_8km']
  [foldername_test '/results_4km']};

for fi = 1: length( foldernames)
  
  foldername = foldernames{ fi};
  filename = [foldername '/transect_southnorth.nc'];
  V = ncread( filename,'V');
  yt = V(:,2);
  u_3D = ncread( filename,'u_ort_3D');
  u_surf = u_3D(:,1);

  % Compare to analytical solution
  u_an = Schoof2006_icestream( 1e-18, 2000, -0.0003, 150e3, 1, yt);

  RMSE(fi) = sqrt( mean( (u_surf - u_an).^2 ));

end

% Read stability info 
filename = [foldername_test '/results_4km/scalar_output_ANT_00001.nc'];
nskip = 0; % Skip the first few values, as the model is still relaxing there
stab = read_stability_info( filename, nskip);

% Write to scoreboard file
write_to_scoreboard_file( ...
  RMSE(1), RMSE(2), RMSE(3), RMSE(4), stab);

  function write_to_scoreboard_file( ...
      RMSE_32km, RMSE_16km, RMSE_8km, RMSE_4km, stab)

    % Set up a scoreboard results structure
    single_run = initialise_single_test_run( test_name, test_path);
  
    % Add cost functions to results structure
    single_run = add_cost_function_to_single_run( single_run, ...
      'RMSE_32km', 'sqrt( mean( (u_surf - u_an).^2 ))', RMSE_32km);
    single_run = add_cost_function_to_single_run( single_run, ...
      'RMSE_16km', 'sqrt( mean( (u_surf - u_an).^2 ))', RMSE_16km);
    single_run = add_cost_function_to_single_run( single_run, ...
      'RMSE_8km', 'sqrt( mean( (u_surf - u_an).^2 ))', RMSE_8km);
    single_run = add_cost_function_to_single_run( single_run, ...
      'RMSE_4km', 'sqrt( mean( (u_surf - u_an).^2 ))', RMSE_4km);

    single_run = add_stability_info_cost_functions( single_run, stab);
    
    % Write to scoreboard file
    write_scoreboard_file( foldername_automated_testing, single_run);
  
  end
  
end