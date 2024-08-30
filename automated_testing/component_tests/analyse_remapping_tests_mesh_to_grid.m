function analyse_remapping_tests_mesh_to_grid( foldername, do_print)

% Create text file with the analysis results
filename_analysis = [foldername '/remapping_analysis.txt'];

disp(['Creating analysis file "' filename_analysis '"...'])

fid = fopen( filename_analysis, 'w');
fprintf( fid, '%s\n', '==================================================================');
fprintf( fid, '%s\n', '===== Analysis of the mesh-to-grid remapping component tests =====');
fprintf( fid, '%s\n', '==================================================================');
fclose( fid);

% List all the test results
filenames = dir( foldername);
i = 1;
while i <= length( filenames)
  if  contains( filenames(i).name, 'results_remapping_') && ...
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
  if do_print

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
    print( H.Fig, [foldername '/' filename_png], '-dpng');
    close( H.Fig)

  end

  write_errors_to_text_file( filename, grid, d_grid_ex, mesh, d_grid, d_mesh_ex);
  
end
  
function write_errors_to_text_file( filename, grid, d_grid_ex, mesh, d_grid, d_mesh_ex)

  int_grid = sum( d_grid(:)) * grid.dx^2;
  int_mesh = sum( d_mesh_ex .* mesh.A);

  fid = fopen( filename_analysis,'a');
  fprintf( fid, '%s\n', '');
  fprintf( fid, '%s\n', strrep( filename, '.nc', ''));
  fprintf( fid, '%s %12.8e\n', '  min( d_mesh_ex)           = ', min( d_mesh_ex(:)));
  fprintf( fid, '%s %12.8e\n', '  max( d_mesh_ex)           = ', max( d_mesh_ex(:)));
  fprintf( fid, '%s %12.8e\n', '  int( d_mesh_ex)           = ', int_mesh);
  fprintf( fid, '%s %12.8e\n', '  min( d_grid)              = ', min( d_grid(:)));
  fprintf( fid, '%s %12.8e\n', '  max( d_grid)              = ', max( d_grid(:)));
  fprintf( fid, '%s %12.8e\n', '  int( d_grid)              = ', int_grid);
  fprintf( fid, '%s %12.8e\n', '  max( d_grid - d_grid_ex)  = ', max( d_grid(:) - d_grid_ex(:)));
  fprintf( fid, '%s %12.8e\n', '  RMS( d_grid - d_grid_ex)  = ', sqrt( mean( (d_grid(:) - d_grid_ex(:)).^2 )));
  fprintf( fid, '%s %12.8e\n', '  global conservation error = ', (int_grid - int_mesh) / int_mesh);
  fclose( fid);

end

end