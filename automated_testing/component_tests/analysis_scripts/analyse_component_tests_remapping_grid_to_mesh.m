function analyse_component_tests_remapping_grid_to_mesh( foldername_automated_testing, do_print_figures)
% Analyse the results of all the grid-to-mesh remapping component tests

disp('    Analysing mesh-to-grid remapping component tests...')
disp('')

foldername_results = [foldername_automated_testing '/component_tests/results/remapping/grid_to_mesh'];
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
  % Analyse the results of the complete set of grid-to-mesh remapping component
  % tests for a single grid-mesh combination

  disp(['      ' filename_short '...']);

  filename_full = [foldername_results '/' filename_short];

  % Read test results
  grid.x    = ncread( filename_full,'x');
  grid.y    = ncread( filename_full,'y');
  grid.dx = grid.x(2) - grid.x(1);
  mesh      = read_mesh_from_file( filename_full);
  d_grid_ex = ncread( filename_full, 'd_grid_ex');
  d_mesh    = ncread( filename_full, 'd_mesh');
  d_mesh_ex = ncread( filename_full, 'd_mesh_ex');
  d_err = d_mesh - d_mesh_ex;

  clim = [-1,1] * 200;

  % Plot error
  if do_print_figures

    wa = 300;
    ha = wa;
    H = setup_multipanel_figure( wa, ha, [5,100], [5, 5]);

    set( H.Ax{1,1},'xtick',[],'ytick',[],'xgrid','off','ygrid','off')
    H.Ax{1,1}.XAxis.Visible = 'off';
    H.Ax{1,1}.YAxis.Visible = 'off';

    H.mesh = patch('parent',H.Ax{1,1},'vertices',mesh.V,'faces',mesh.Tri,...
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

  write_to_scoreboard_file( filename_short, grid, d_grid_ex, mesh, d_mesh_ex, d_mesh)
  
end

function write_to_scoreboard_file( filename_short, grid, d_grid_ex, mesh, d_mesh_ex, d_mesh)

  % Set up a scoreboard results structure
  test_name = filename_short(5:end-3);
  single_run = initialise_single_test_run( test_name, ...
    'component_tests/remapping/grid_to_mesh');

  % Calculate cost functions
  rmse = sqrt( mean( (d_mesh - d_mesh_ex).^2));

  bounds_max = max( 0, max( d_mesh(:)) - max( d_grid_ex(:)));
  bounds_min = max( 0, min( d_grid_ex(:)) - min( d_mesh(:)));

  int_mesh = sum( d_mesh .* mesh.A);
  int_grid = sum( d_grid_ex(:)) * grid.dx^2;
  int_err = abs( 1 - int_mesh / int_grid);

  % Add cost functions to results structure
  single_run = add_cost_function_to_single_run( single_run, 'rmse'       , 'sqrt( mean( (d_mesh - d_mesh_ex).^2))'        , rmse);
  single_run = add_cost_function_to_single_run( single_run, 'bounds_max' , 'max( 0, max( d_mesh(:)) - max( d_grid_ex(:)))', bounds_max);
  single_run = add_cost_function_to_single_run( single_run, 'bounds_min' , 'max( 0, min( d_grid_ex(:)) - min( d_mesh(:)))', bounds_min);
  single_run = add_cost_function_to_single_run( single_run, 'int_err'    , 'abs( 1 - int_mesh / int_grid)'                , int_err);

  % Write to scoreboard file
  all_runs_new.single_run = single_run;
  filename_scoreboard_file = [foldername_automated_testing ...
    '/component_tests/temporary_scoreboard_files/' ...
    strrep( single_run.category, '/', '_') '_' test_name '.xml'];
  append_test_results_to_main_scoreboard_file( all_runs_new, filename_scoreboard_file);

end
  
end