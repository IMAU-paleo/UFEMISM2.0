function analyse_integrated_test( varargin)
% Analyse the results of the Halfar_5km integrated test

disp('Analysing integrated test idealised/Halfar_dome/Halfar_5km...')
disp('')
  
test_name = 'Halfar_5km';
test_path = 'integrated_tests/idealised/Halfar_dome';

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

filename_config  = [foldername_test '/config.cfg'];
filename_results = [foldername_test '/results/main_output_ANT_00001.nc'];

addpath([foldername_automated_testing '/scoreboard/scripts'])

%%

% Read model output
mesh      = read_mesh_from_file( filename_results);
time      = ncread( filename_results,'time'); nt = length( time);
Hi        = ncread( filename_results,'Hi');

% Read Halfar dome parameters from config file, calculate analytical solution
[A_flow, n_flow, H0, R0] = get_Halfar_dome_params_from_config( filename_config);

Hi_analytical = zeros( mesh.nV, nt);
for ti = 1: nt
  Hi_analytical( :,ti) = Halfar_solution( A_flow, n_flow, H0, R0, mesh.V(:,1), mesh.V(:,2), time(ti));
end

% Calculate thickness error at the end of the simulation
RMSE_Hi = sqrt( mean( (Hi( :,end) - Hi_analytical( :,end)).^2 ));

% Write to scoreboard file
write_to_scoreboard_file( RMSE_Hi);

  function [A_flow, n_flow, H0, R0] = get_Halfar_dome_params_from_config( filename_config)
    fid = fopen( filename_config,'r');
    temp = textscan(fid,'%s','delimiter','\n'); temp = temp{1};
    fclose( fid);
    
    found_A_flow = false;
    found_n_flow = false;
    found_H0     = false;
    found_R0     = false;
    
    for i = 1: length( temp)
      str = temp{i};
    
      if contains( str,'uniform_Glens_flow_factor_config')
        found_A_flow = true;
        ii = strfind( str,'='); ii = ii(1);
        str = str( ii+1:end);
        ii = strfind( str,'!'); ii = ii(1);
        str = str( 1:ii-1);
        str = strip( str,'both');
        A_flow = str2double( str);
      end
    
      if contains( str,'Glens_flow_law_exponent_config')
        found_n_flow = true;
        ii = strfind( str,'='); ii = ii(1);
        str = str( ii+1:end);
        ii = strfind( str,'!'); ii = ii(1);
        str = str( 1:ii-1);
        str = strip( str,'both');
        n_flow = str2double( str);
      end
    
      if contains( str,'refgeo_idealised_Halfar_H0_config')
        found_H0 = true;
        ii = strfind( str,'='); ii = ii(1);
        str = str( ii+1:end);
        ii = strfind( str,'!'); ii = ii(1);
        str = str( 1:ii-1);
        str = strip( str,'both');
        H0 = str2double( str);
      end
    
      if contains( str,'refgeo_idealised_Halfar_R0_config')
        found_R0 = true;
        ii = strfind( str,'='); ii = ii(1);
        str = str( ii+1:end);
        ii = strfind( str,'!'); ii = ii(1);
        str = str( 1:ii-1);
        str = strip( str,'both');
        R0 = str2double( str);
      end
    
    end

    if ~found_A_flow
      error('Couldnt find A_flow in config file!')
    end
    if ~found_n_flow
      error('Couldnt find n_flow in config file!')
    end
    if ~found_H0
      error('Couldnt find H0 in config file!')
    end
    if ~found_R0
      error('Couldnt find R0 in config file!')
    end

  end

  function H = Halfar_solution( A_flow, n_flow, H0, R0, x, y, time)
  
  ice_density  = 910;
  grav         = 9.81;
  sec_per_year = 31556926;
  
  Gamma = (2/5) * (A_flow / sec_per_year) * (ice_density * grav)^n_flow;
  t0 = 1 / ((5*n_flow + 3) * Gamma) * ((2*n_flow + 1)/(n_flow + 1))^n_flow * (R0^(n_flow + 1)) / (H0^(2*n_flow  + 1));

  tp = (time * sec_per_year) + t0;

  r = sqrt(x.^2 + y.^2);

  f1 = (t0 / tp)^(2 / (5*n_flow + 3));
  f2 = (t0 / tp)^(1 / (5*n_flow + 3));
  f3 = (r / R0);

  H = H0 * f1 * max( 0, (1 - (f2 * f3).^((n_flow + 1) / n_flow))).^(n_flow / (2*n_flow + 1));
end

  function write_to_scoreboard_file( RMSE_Hi)

    % Set up a scoreboard results structure
    single_run = initialise_single_test_run( test_name, test_path);
  
    % Add cost functions to results structure
    single_run = add_cost_function_to_single_run( single_run, ...
      'rmse', 'sqrt( mean( (Hi( :,end) - Hi_analytical( :,end)).^2 ))', RMSE_Hi);
    
    % Write to temporary scoreboard file
    filename_temporary_scoreboard_file = [foldername_test ...
      '/scoreboard_temp_' strrep( test_path, '/', '_') '_' test_name '.xml'];
    all_runs.single_run = single_run;
    write_scoreboard_file( all_runs, filename_temporary_scoreboard_file);
  
  end
  
end