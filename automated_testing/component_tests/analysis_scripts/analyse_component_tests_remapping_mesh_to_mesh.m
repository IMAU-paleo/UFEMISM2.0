function analyse_component_tests_remapping_mesh_to_mesh( foldername_automated_testing, do_print_figures)
% Analyse the results of all the mesh-to-mesh remapping component tests

disp('    Analysing mesh-to-mesh remapping component tests...')
disp('')

foldername_results = [foldername_automated_testing '/component_tests/results/remapping/mesh_to_mesh'];
foldername_figures = [foldername_automated_testing '/component_tests/figures'];

% List all the test results
filenames = dir( foldername_results);
i = 1;
while i <= length( filenames)
  if  contains( filenames(i).name, 'res_') && ...
      contains( filenames(i).name, '.nc')
    i = i+1;
  else
    filenames( i) = [];
  end
end

for fi = 1: length( filenames)
  analyse_remapping_test( filenames( fi).name)
end

disp('')

function analyse_remapping_test( filename_short)
  % Analyse the results of the complete set of mesh-to-mesh remapping component
  % tests for a single mesh-mesh combination

  disp(['      ' filename_short '...']);

  filename_full = [foldername_results '/' filename_short];

  % Read test results
  mesh1.A        = ncread( filename_full, 'mesh1_A');
  mesh2          = read_mesh_from_file( filename_full);
  d_mesh1_ex     = ncread( filename_full, 'd_mesh1_ex');
  d_mesh2_ex     = ncread( filename_full, 'd_mesh2_ex');
  d_mesh2_nn     = ncread( filename_full, 'd_mesh2_nn');
  d_mesh2_trilin = ncread( filename_full, 'd_mesh2_trilin');
  d_mesh2_cons   = ncread( filename_full, 'd_mesh2_cons');

  d_err_nn     = d_mesh2_nn     - d_mesh2_ex;
  d_err_trilin = d_mesh2_trilin - d_mesh2_ex;
  d_err_cons   = d_mesh2_cons   - d_mesh2_ex;

  clim = [-1,1] * 200;

  % Plot error
  if do_print_figures

    wa = 300;
    ha = wa;
    H = setup_multipanel_figure( wa, ha, [5,100], [5, 5]);

    set( H.Ax{1,1},'xtick',[],'ytick',[],'xgrid','off','ygrid','off')
    H.Ax{1,1}.XAxis.Visible = 'off';
    H.Ax{1,1}.YAxis.Visible = 'off';

    H.mesh = patch('parent',H.Ax{1,1},'vertices',mesh2.V,'faces',mesh2.Tri,...
      'facevertexcdata',d_err,'facecolor','interp','edgecolor','none');
  
    pos = get(H.Ax{1,1},'position');
    x_lo = pos(1) + pos(3) + 5;
    pos = get(H.Fig,'position');
    x_hi = pos(3) - 5;
    y_lo = 5;
    y_hi = pos(4) - 5;
    wa = x_hi - x_lo;
    ha = y_hi - y_lo;
    H.Cbaraxes = axes('parent',H.Fig,'units','pixels','position',[x_lo,y_lo,wa,ha],...
      'fontsize',get(H.Ax{1,1},'fontsize'));
    H.Cbaraxes.XAxis.Visible = 'off';
    H.Cbaraxes.YAxis.Visible = 'off';
    H.Cbar = colorbar(H.Cbaraxes,'location','west');
  
    cmap = bluewhiteredmap( 256);
    colormap( H.Ax{1,1},cmap);
    colormap( H.Cbaraxes,cmap);
  
    set( H.Ax{1,1},'clim',clim);
    set( H.Cbaraxes,'clim',clim);

    filename_png = strrep( filename_short, '.nc', '.png');
    print( H.Fig, [foldername_figures '/' filename_png], '-dpng');
    close( H.Fig)

  end

  write_to_scoreboard_file( filename_short, mesh1, d_mesh1_ex, mesh2, d_mesh2_ex, ...
    d_mesh2_nn, d_mesh2_trilin, d_mesh2_cons)
  
end

  function write_to_scoreboard_file( filename_short, mesh1, d_mesh1_ex, mesh2, d_mesh2_ex, ...
      d_mesh2_nn, d_mesh2_trilin, d_mesh2_cons)

  % Set up a scoreboard results structure
  test_name = filename_short(5:end-3);
  single_run = initialise_single_test_run( test_name, ...
    'component_tests/remapping/mesh_to_mesh');

  % Calculate cost functions
  rmse_nn     = sqrt( mean( (d_mesh2_nn     - d_mesh2_ex).^2));
  rmse_trilin = sqrt( mean( (d_mesh2_trilin - d_mesh2_ex).^2));
  rmse_cons   = sqrt( mean( (d_mesh2_cons   - d_mesh2_ex).^2));

  bounds_max_nn     = max( 0, max( d_mesh2_nn(    :)) - max( d_mesh1_ex(:)));
  bounds_max_trilin = max( 0, max( d_mesh2_trilin(:)) - max( d_mesh1_ex(:)));
  bounds_max_cons   = max( 0, max( d_mesh2_cons(  :)) - max( d_mesh1_ex(:)));

  bounds_min_nn     = max( 0, min( d_mesh2_ex(:)) - min( d_mesh2_nn(    :)));
  bounds_min_trilin = max( 0, min( d_mesh2_ex(:)) - min( d_mesh2_trilin(:)));
  bounds_min_cons   = max( 0, min( d_mesh2_ex(:)) - min( d_mesh2_cons(  :)));

  int_mesh1 = sum( d_mesh1_ex .* mesh1.A);

  int_mesh2_nn     = sum( d_mesh2_nn     .* mesh2.A);
  int_mesh2_trilin = sum( d_mesh2_trilin .* mesh2.A);
  int_mesh2_cons   = sum( d_mesh2_cons   .* mesh2.A);

  int_err_nn     = abs( 1 - int_mesh2_nn     / int_mesh1);
  int_err_trilin = abs( 1 - int_mesh2_trilin / int_mesh1);
  int_err_cons   = abs( 1 - int_mesh2_cons   / int_mesh1);

  % Add cost functions to results structure
  single_run = add_cost_function_to_single_run( single_run, 'rmse_nn'           , 'sqrt( mean( (d_mesh2_nn     - d_mesh2_ex).^2))'        , rmse_nn);
  single_run = add_cost_function_to_single_run( single_run, 'rmse_trilin'       , 'sqrt( mean( (d_mesh2_trilin - d_mesh2_ex).^2))'        , rmse_trilin);
  single_run = add_cost_function_to_single_run( single_run, 'rmse_cons'         , 'sqrt( mean( (d_mesh2_cons   - d_mesh2_ex).^2))'        , rmse_cons);

  single_run = add_cost_function_to_single_run( single_run, 'bounds_max_nn'     , 'max( 0, max( d_mesh2_nn(    :)) - max( d_mesh1_ex(:)))', bounds_max_nn);
  single_run = add_cost_function_to_single_run( single_run, 'bounds_max_trilin' , 'max( 0, max( d_mesh2_trilin(:)) - max( d_mesh1_ex(:)))', bounds_max_trilin);
  single_run = add_cost_function_to_single_run( single_run, 'bounds_max_cons'   , 'max( 0, max( d_mesh2_cons(  :)) - max( d_mesh1_ex(:)))', bounds_max_cons);

  single_run = add_cost_function_to_single_run( single_run, 'bounds_min_nn'     , 'max( 0, min( d_mesh2_ex(:)) - min( d_mesh2_nn(    :)))', bounds_min_nn);
  single_run = add_cost_function_to_single_run( single_run, 'bounds_min_trilin' , 'max( 0, min( d_mesh2_ex(:)) - min( d_mesh2_trilin(:)))', bounds_min_trilin);
  single_run = add_cost_function_to_single_run( single_run, 'bounds_min_cons'   , 'max( 0, min( d_mesh2_ex(:)) - min( d_mesh2_cons(  :)))', bounds_min_cons);

  single_run = add_cost_function_to_single_run( single_run, 'int_err_nn'        , 'abs( 1 - int_mesh2_nn     / int_mesh1)'                , int_err_nn);
  single_run = add_cost_function_to_single_run( single_run, 'int_err_trilin'    , 'abs( 1 - int_mesh2_trilin / int_mesh1)'                , int_err_trilin);
  single_run = add_cost_function_to_single_run( single_run, 'int_err_cons'      , 'abs( 1 - int_mesh2_cons   / int_mesh1)'                , int_err_cons);

  % Write to scoreboard file
  write_scoreboard_file( foldername_automated_testing, single_run);

end
  
end