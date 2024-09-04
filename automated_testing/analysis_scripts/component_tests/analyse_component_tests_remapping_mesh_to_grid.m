function analyse_component_tests_remapping_mesh_to_grid( foldername, foldername_automated_testing, do_print_figures)
% Analyse the results of all the mesh-to-grid remapping component tests

% List all the test results
filenames = dir( foldername);
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
  analyse_remapping_test( foldername, filenames( fi).name)
end

function analyse_remapping_test( foldername, filename)
  % Analyse the results of the complete set of mapping/derivative component
  % tests for a single mesh and a single test function

  disp(['Analysing ' filename '...']);

  filename_full = [foldername '/' filename];

  % Read test results
  grid.x    = ncread(filename_full,'x');
  grid.y    = ncread(filename_full,'y');
  grid.dx = grid.x(2) - grid.x(1);
  mesh      = read_mesh_from_file( filename_full);
  d_grid_ex = ncread( filename_full, 'd_grid_ex');
  d_grid    = ncread( filename_full, 'd_grid');
  d_mesh_ex = ncread( filename_full, 'd_mesh_ex');
  d_err = d_grid - d_grid_ex;

  clim = [-1,1] * 200;

  % Plot error
  if do_print_figures

    wa = 300;
    ha = wa;
    H = setup_multipanel_figure( wa, ha, [5,100], [5, 5]);

    set( H.Ax{1,1},'xtick',[],'ytick',[],'xgrid','off','ygrid','off')
    H.Ax{1,1}.XAxis.Visible = 'off';
    H.Ax{1,1}.YAxis.Visible = 'off';

    H.im = image('parent',H.Ax{1,1},'xdata',grid.x,'ydata',grid.y,...
      'cdata',d_err,'cdatamapping','scaled');
  
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

    filename_png = strrep( filename, '.nc', '.png');
    print( H.Fig, [foldername_automated_testing '/figures/' filename_png], '-dpng');
    close( H.Fig)

  end

  write_to_scoreboard( filename, mesh, d_mesh_ex, grid, d_grid_ex, d_grid);
  
end
  
function write_to_scoreboard( filename, mesh, d_mesh_ex, grid, d_grid_ex, d_grid)

  % Set up a scoreboard results structure
  test_name = ['remapping_mesh_to_grid_' filename(5:end-3)];
  res = initialise_test_results( test_name, 'remapping/mesh_to_grid');

  % Calculate cost functions
  rmse = sqrt( mean( (d_grid - d_grid_ex).^2));

  bounds_max = max( 0, max( d_grid(:)) - max( d_mesh_ex(:)));
  bounds_min = max( 0, min( d_mesh_ex(:)) - min( d_grid(:)));

  int_grid = sum( d_grid(:)) * grid.dx^2;
  int_mesh = sum( d_mesh_ex .* mesh.A);
  int_err = abs( 1 - int_grid / int_mesh);

  % Add cost functions to results structure
  res = add_result_to_test_results( res, 'rmse'       , 'sqrt( mean( (d_grid - d_grid_ex).^2))'        , rmse);
  res = add_result_to_test_results( res, 'bounds_max' , 'max( 0, max( d_grid(:)) - max( d_mesh_ex(:)))', bounds_max);
  res = add_result_to_test_results( res, 'bounds_min' , 'max( 0, min( d_mesh_ex(:)) - min( d_grid(:)))', bounds_min);
  res = add_result_to_test_results( res, 'int_err'    , 'abs( 1 - int_grid / int_mesh)'                , int_err);

  % Write to scoreboard file
  write_test_results_to_scoreboard_file( res, 'scoreboard');

end
  
end