function analyse_integrated_test_ISMIP_HOM_SIASSA( foldername_automated_testing, foldername_ISMIP_HOM)
% Analyse the results of the Halfar_5km integrated test


%%

% Read and process the Pattyn et al. (2008) ISMIP-HOM model ensemble
[HO,FS] = process_ISMIP_HOM_ensemble_experiment_A( [foldername_ISMIP_HOM '/ismip_all']);

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

    foldername = ['results_ISMIP_HOM_A_' num2str(L) '_' approx];
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
    test_path = ['integrated_tests/idealised/ISMIP_HOM/experiment_A/' approx];
    test_name = ex;
    write_to_scoreboard_file( RMSE_u_surf)

  end
end

  function write_to_scoreboard_file( RMSE_Hi)

    % Set up a scoreboard results structure
    single_run = initialise_single_test_run( test_name, test_path);
  
    % Add cost functions to results structure
    single_run = add_cost_function_to_single_run( single_run, ...
      'rmse', 'sqrt( mean( (Hi( :,end) - Hi_analytical( :,end)).^2 ))', RMSE_Hi);
    
    % Write to scoreboard file
    write_scoreboard_file( foldername_automated_testing, single_run);
  
  end
  
end