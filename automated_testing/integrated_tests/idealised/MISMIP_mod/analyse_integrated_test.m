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

filename_results1 = [foldername_test '/results_spinup_10km/main_output_ANT_00001.nc'];
filename_results2 = [foldername_test '/results_retreat_10km/main_output_ANT_00012.nc'];

% Read model output
mesh1     = read_mesh_from_file( filename_results1);
mesh2     = read_mesh_from_file( filename_results2);
time1     = ncread( filename_results1,'time'); nt1 = length( time1);
time2     = ncread( filename_results2,'time'); nt2 = length( time2);
Hi1       = ncread( filename_results1,'Hi',[1,nt1],[Inf,1]);
Hi2       = ncread( filename_results2,'Hi',[1,nt2],[Inf,1]);
Hb1       = ncread( filename_results1,'Hb',[1,nt1],[Inf,1]);
Hb2       = ncread( filename_results2,'Hb',[1,nt2],[Inf,1]);

r_GL1 = calc_mean_grounding_line_radius( mesh1, Hi1, Hb1);
r_GL2 = calc_mean_grounding_line_radius( mesh2, Hi2, Hb2);

GL_hyst = r_GL2 - r_GL1;

% Read stability info for the last spin-up simulation
filename = 'results_spinup_10km/scalar_output_ANT_00001.nc';
nskip = 5; % Skip the first few values, as the model is still relaxing there
stab = read_stability_info( filename, nskip);

% Write to scoreboard file
write_to_scoreboard_file( GL_hyst, stab);

  function r_GL = calc_mean_grounding_line_radius( mesh, Hi, Hb)
    TAF = Hi - max( 0, -Hb * (1028 / 910));
    n_GL = 0;
    sum_r_GL = 0;
    for ei = 1: mesh.nE
      vi = mesh.EV( ei,1);
      vj = mesh.EV( ei,2);
      if TAF(vi)*TAF(vj) < 0
        f1 = TAF(vi);
        f2 = TAF(vj);
        d = norm(mesh.V(vi,:) - mesh.V(vj,:));
        lambda = (f2 - f1) / d;
        p_GL = lambda * mesh.V(vi,:) + (1-lambda) * mesh.V(vj,:);
        r_GL = norm(p_GL);
        n_GL = n_GL + 1;
        sum_r_GL = sum_r_GL + r_GL;
      end
    end
    r_GL = sum_r_GL / n_GL;
  end

  function write_to_scoreboard_file( GL_hyst, stab)

    % Set up a scoreboard results structure
    single_run = initialise_single_test_run( test_name, test_path);
  
    % Add cost functions to results structure
    single_run = add_cost_function_to_single_run( single_run, ...
      'GL_hyst', 'r_GL2 - r_GL1', GL_hyst);

    single_run = add_stability_info_cost_functions( single_run, stab);
    
    % Write to scoreboard file
    write_scoreboard_file( foldername_automated_testing, single_run);
  
  end
  
end