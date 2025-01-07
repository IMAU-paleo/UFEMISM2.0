function analyse_integrated_test( varargin)
% Analyse the results of the ISMIP-HOM experiment B integrated test

varargin = varargin{1};

disp('Analysing integrated test idealised/ISMIP-HOM/experiment_B...')
disp('')

%%

% In the GitHub Workflow, provide the component test results folders as
% input; but retain the option of running without input (i.e. as a script)
% locally, with user-defined folders

input_args = varargin;
if isempty( input_args)
  % Assume this is a local run

  foldername_automated_testing = '/Users/Beren017/Documents/GitHub/UFEMISM2.0/automated_testing';
  addpath('/Users/Beren017/Documents/GitHub/UFEMISM2.0/tools/matlab/')
  foldername_ISMIP_HOM = '/Users/Beren017/Documents/GitHub/data/model_ensembles/ISMIP-HOM';

elseif isscalar( input_args)
  % Assume this is a GitHub Workflow run

  foldername_automated_testing = varargin{1};
  addpath('tools/matlab/')
  foldername_ISMIP_HOM = 'data/model_ensembles/ISMIP-HOM';

else
  error('need either foldername_automated_testing, or nothing as input!')
end

addpath([foldername_automated_testing '/scoreboard/scripts'])
addpath(foldername_ISMIP_HOM)

foldername_test = [foldername_automated_testing '/integrated_tests/idealised/' ...
  'ISMIP-HOM/exp_B'];

%%

% Read and process the Pattyn et al. (2008) ISMIP-HOM model ensemble
[HO,FS] = process_ISMIP_HOM_ensemble_experiment_B( [foldername_ISMIP_HOM '/ismip_all']);

% Read and process UFEMISM output
Ls = [5,10,20,40,80,160];
approxs = {'SIASSA','DIVA','BPA'};

for Li = 1:6
  L = Ls( Li);
  
  if L<10
    ex = ['L00' num2str(L)];
  elseif L<100
    ex = ['L0'  num2str(L)];
  else
    ex = ['L'   num2str(L)];
  end

  for approxi = 1:3
    approx = approxs{ approxi};

    foldername = [foldername_test '/results_ISMIP_HOM_B_' num2str(L) '_' approx];
    filename = [foldername '/main_output_ANT_00001.nc'];
    mesh = read_mesh_from_file( filename);

    % Calculate velocity transect
    xt = HO.(ex).x * L*1e3;
    yt = xt*0 + L*1e3 / 4;
    A = calc_transect_matrix( mesh, xt, yt);
    u_surf = ncread( filename,'u_surf');
    u_surf = A * u_surf;

    % Calculate error with respect to HO model ensemble
    err_pos = max( 0, u_surf - HO.(ex).u_max);
    err_neg = max( 0, HO.(ex).u_min - u_surf);
    err_abs = max( err_pos, err_neg);
    RMSE_u_surf = sqrt( mean( err_abs.^2));
  
    % Write to scoreboard file
    test_path = ['integrated_tests/idealised/ISMIP_HOM/experiment_B/' approx];
    test_name = ex;
    write_to_scoreboard_file( RMSE_u_surf)

  end
end

  function write_to_scoreboard_file( RMSE_u_surf)

    % Set up a scoreboard results structure
    single_run = initialise_single_test_run( test_name, test_path);
  
    % Add cost functions to results structure
    single_run = add_cost_function_to_single_run( single_run, ...
      'rmse', 'sqrt( mean( (u_surf - u_ensemble).^2 ))', RMSE_u_surf);
    
    % Write to scoreboard file
    write_scoreboard_file( foldername_automated_testing, single_run);
  
  end
end