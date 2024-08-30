function analyse_mapping_derivative_tests( foldername, do_print)

% Create text file with the analysis results
filename_analysis = [foldername '/map_deriv_analysis.txt'];

disp(['Creating analysis file "' filename_analysis '"...'])

fid = fopen( filename_analysis, 'w');
fprintf( fid, '%s\n', '===============================================================');
fprintf( fid, '%s\n', '===== Analysis of the mapping/derivatives component tests =====');
fprintf( fid, '%s\n', '===============================================================');
fclose( fid);

% List all the test results
filenames = dir( foldername);
i = 1;
while i <= length( filenames)
  if  contains( filenames(i).name, 'results_map_deriv_tests_') && ...
      contains( filenames(i).name, '.nc')
    i = i+1;
  else
    filenames( i) = [];
  end
end

for fi = 1: length( filenames)
  analyse_mapping_derivatives_tests_mesh_function( foldername, filenames( fi).name)
end

function analyse_mapping_derivatives_tests_mesh_function( foldername, filename)
  % Analyse the results of the complete set of mapping/derivative component
  % tests for a single mesh and a single test function

  disp(['Analysing ' filename '...'])
  
  [H_map, err_max_map] = analyse_mapping_derivatives_tests_mesh_function_map(  [foldername '/' filename]);
  
  [H_ddx, err_max_ddx] = analyse_mapping_derivatives_tests_mesh_function_ddxy( [foldername '/' filename],'ddx');
  
  [H_ddy, err_max_ddy] = analyse_mapping_derivatives_tests_mesh_function_ddxy( [foldername '/' filename],'ddy');
  
  [H_2nd, err_max_ddx_2nd, err_max_ddy_2nd, err_max_d2dx2_2nd, err_max_d2dxdy_2nd, err_max_d2dy2_2nd] = ...
    analyse_mapping_derivatives_tests_mesh_function_2nd( [foldername '/' filename]);
  
  write_errors_to_text_file( filename, err_max_map, err_max_ddx, err_max_ddy, ...
    err_max_ddx_2nd, err_max_ddy_2nd, err_max_d2dx2_2nd, err_max_d2dxdy_2nd, err_max_d2dy2_2nd);

  if do_print
    filename_map_png = strrep( filename, '.nc', '_map.png');
    filename_ddx_png = strrep( filename, '.nc', '_ddx.png');
    filename_ddy_png = strrep( filename, '.nc', '_ddy.png');
    filename_2nd_png = strrep( filename, '.nc', '_2nd.png');
  
    print( H_map.Fig, [foldername '/' filename_map_png], '-dpng');
    print( H_ddx.Fig, [foldername '/' filename_ddx_png], '-dpng');
    print( H_ddy.Fig, [foldername '/' filename_ddy_png], '-dpng');
    print( H_2nd.Fig, [foldername '/' filename_2nd_png], '-dpng');
  
    close( H_map.Fig)
    close( H_ddx.Fig)
    close( H_ddy.Fig)
    close( H_2nd.Fig)
  end
  
end
  
function write_errors_to_text_file( filename, err_max_map, err_max_ddx, err_max_ddy, ...
    err_max_ddx_2nd, err_max_ddy_2nd, err_max_d2dx2_2nd, err_max_d2dxdy_2nd, err_max_d2dy2_2nd)

  fid = fopen( filename_analysis,'a');
  fprintf( fid, '%s\n', '');
  fprintf( fid, '%s\n', strrep( filename, '.nc', ''));
  fprintf( fid, '%s %12.8e\n', '  err_max_map        = ', err_max_map);
  fprintf( fid, '%s %12.8e\n', '  err_max_ddx        = ', err_max_ddx);
  fprintf( fid, '%s %12.8e\n', '  err_max_ddy        = ', err_max_ddy);
  fprintf( fid, '%s %12.8e\n', '  err_max_ddx_2nd    = ', err_max_ddx_2nd);
  fprintf( fid, '%s %12.8e\n', '  err_max_ddy_2nd    = ', err_max_ddy_2nd);
  fprintf( fid, '%s %12.8e\n', '  err_max_d2dx2_2nd  = ', err_max_d2dx2_2nd);
  fprintf( fid, '%s %12.8e\n', '  err_max_d2dxdy_2nd = ', err_max_d2dxdy_2nd);
  fprintf( fid, '%s %12.8e\n', '  err_max_d2dy2_2nd  = ', err_max_d2dy2_2nd);
  fclose( fid);

end

function [H, err_max] = analyse_mapping_derivatives_tests_mesh_function_map( filename)

  H = [];

  % Read the mesh
  mesh = read_mesh_from_file( filename);
  
  %% Set up the plot
  if do_print

    wa = 150;
    ha = 150;
    margins_hor = [50, 50, 5, 5, 50, 5, 5, 150];
    margins_ver = [75, 5, 5, 5];
    H = setup_multipanel_figure( wa, ha, margins_hor, margins_ver);
    
    for i = 1: size( H.Ax,1)
      for j = 1: size( H.Ax,2)
        set( H.Ax{ i,j},'xtick',[],'ytick',[],'clim',[-1,1]);
      end
    end
  
    for i = 1: 1
      for j = 1: size( H.Ax,2)
        set( H.Ax{ i,j},'xaxislocation','top')
      end
    end
  
    title( H.Ax{ 1,1},'Exact solution')
    xlabel( H.Ax{ 1,1},'-')
    title( H.Ax{ 1,3},'Discretised approximation')
    title( H.Ax{ 1,6},'Discretisation error')
  
    ylabel( H.Ax{ 1,1},'A')
    ylabel( H.Ax{ 2,1},'B')
    ylabel( H.Ax{ 3,1},'C')
  
    xlabel( H.Ax{ 1,2},'from A')
    xlabel( H.Ax{ 1,3},'from B')
    xlabel( H.Ax{ 1,4},'from C')
  
    ylabel( H.Ax{ 1,2},'to A')
    ylabel( H.Ax{ 2,2},'to B')
    ylabel( H.Ax{ 3,2},'to C')
  
    xlabel( H.Ax{ 1,5},'from A')
    xlabel( H.Ax{ 1,6},'from B')
    xlabel( H.Ax{ 1,7},'from C')
  
    ylabel( H.Ax{ 1,5},'to A')
    ylabel( H.Ax{ 2,5},'to B')
    ylabel( H.Ax{ 3,5},'to C')
  
    for i = 1:3
      for j = 5:7
        colormap( H.Ax{ i,j},turbo(256));
      end
    end

  end
  
  %% Read data
  d_a_ex = ncread( filename,'d_a_ex');
  d_b_ex = ncread( filename,'d_b_ex');
  d_c_ex = ncread( filename,'d_c_ex');
  
  % d_a_a  = ncread( filename,'d_a_a');
  d_a_b  = ncread( filename,'d_a_b');
  d_a_c  = ncread( filename,'d_a_c');
  
  d_b_a  = ncread( filename,'d_b_a');
  % d_b_b  = ncread( filename,'d_b_b');
  d_b_c  = ncread( filename,'d_b_c');
  
  d_c_a  = ncread( filename,'d_c_a');
  d_c_b  = ncread( filename,'d_c_b');
  % d_c_c  = ncread( filename,'d_c_c');

  %% Calculate errors
  % err_a_a = calc_log_discretisation_error( d_a_ex, d_a_a);
  err_a_b = calc_log10_discretisation_error( d_b_ex, d_a_b);
  err_a_c = calc_log10_discretisation_error( d_c_ex, d_a_c);

  err_b_a = calc_log10_discretisation_error( d_a_ex, d_b_a);
  % err_b_b = calc_log_discretisation_error( d_b_ex, d_b_b);
  err_b_c = calc_log10_discretisation_error( d_c_ex, d_b_c);

  err_c_a = calc_log10_discretisation_error( d_a_ex, d_c_a);
  err_c_b = calc_log10_discretisation_error( d_b_ex, d_c_b);
  % err_c_c = calc_log_discretisation_error( d_c_ex, d_c_c);

  %% Calculate ranges

  % Calculate range for data
  clim = calc_color_limits( d_a_ex);

  % Calculate range for errors
  err_min = min([ ...
    % min( err_a_a( err_a_a > -Inf & mesh.VBI   == 0))
    min( err_a_b( err_a_b > -Inf & mesh.TriBI == 0))
    min( err_a_c( err_a_c > -Inf & mesh.EBI   == 0))
    min( err_b_a( err_b_a > -Inf & mesh.VBI   == 0))
    % min( err_b_b( err_b_b > -Inf & mesh.TriBI == 0))
    min( err_b_c( err_b_c > -Inf & mesh.EBI   == 0))
    min( err_c_a( err_c_a > -Inf & mesh.VBI   == 0))
    min( err_c_b( err_c_b > -Inf & mesh.TriBI == 0))
    % min( err_c_c( err_c_c > -Inf & mesh.EBI   == 0))
    ]);
  err_max = max([ ...
    % max( err_a_a( err_a_a > -Inf & mesh.VBI   == 0))
    max( err_a_b( err_a_b > -Inf & mesh.TriBI == 0))
    max( err_a_c( err_a_c > -Inf & mesh.EBI   == 0))
    max( err_b_a( err_b_a > -Inf & mesh.VBI   == 0))
    % max( err_b_b( err_b_b > -Inf & mesh.TriBI == 0))
    max( err_b_c( err_b_c > -Inf & mesh.EBI   == 0))
    max( err_c_a( err_c_a > -Inf & mesh.VBI   == 0))
    max( err_c_b( err_c_b > -Inf & mesh.TriBI == 0))
    % max( err_c_c( err_c_c > -Inf & mesh.EBI   == 0))
    ]);
  clim_err = [err_min, err_max];

  %% Plot data and errors
  
  if do_print

    plot_mesh_data_a( H.Ax{ 1,1}, mesh, d_a_ex, clim);
    % plot_mesh_data_a( H.Ax{ 1,2}, mesh, d_a_a, clim);
    plot_mesh_data_a( H.Ax{ 1,3}, mesh, d_b_a, clim);
    plot_mesh_data_a( H.Ax{ 1,4}, mesh, d_c_a, clim);
    % plot_mesh_data_a( H.Ax{ 1,5}, mesh, err_a_a, clim_err);
    plot_mesh_data_a( H.Ax{ 1,6}, mesh, err_b_a, clim_err);
    plot_mesh_data_a( H.Ax{ 1,7}, mesh, err_c_a, clim_err);
    
    plot_mesh_data_b( H.Ax{ 2,1}, mesh, d_b_ex, clim);
    plot_mesh_data_b( H.Ax{ 2,2}, mesh, d_a_b, clim);
    % plot_mesh_data_b( H.Ax{ 2,3}, mesh, d_b_b, clim);
    plot_mesh_data_b( H.Ax{ 2,4}, mesh, d_c_b, clim);
    plot_mesh_data_b( H.Ax{ 2,5}, mesh, err_a_b, clim_err);
    % plot_mesh_data_b( H.Ax{ 2,6}, mesh, err_b_b, clim_err);
    plot_mesh_data_b( H.Ax{ 2,7}, mesh, err_c_b, clim_err);
  
    plot_mesh_data_c( H.Ax{ 3,1}, mesh, d_c_ex, clim, colormap( H.Ax{ 3,1}));
    plot_mesh_data_c( H.Ax{ 3,2}, mesh, d_a_c, clim, colormap( H.Ax{ 3,2}));
    plot_mesh_data_c( H.Ax{ 3,3}, mesh, d_b_c, clim, colormap( H.Ax{ 3,3}));
    % plot_mesh_data_c( H.Ax{ 3,4}, mesh, d_c_c, clim, colormap( H.Ax{ 3,4}));
    plot_mesh_data_c( H.Ax{ 3,5}, mesh, err_a_c, clim_err, colormap( H.Ax{ 3,5}));
    plot_mesh_data_c( H.Ax{ 3,6}, mesh, err_b_c, clim_err, colormap( H.Ax{ 3,6}));
    % plot_mesh_data_c( H.Ax{ 3,7}, mesh, err_c_c, clim_err, colormap( H.Ax{ 3,7}));
  
    %% Add colorbar for errors
    pos = get( H.Ax{ 1,7},'position');
    x_left  = pos( 1) + pos( 3) + 5;
    y_top   = pos( 2) + pos( 4);
    pos = get( H.Ax{ 3,7},'position');
    y_bot   = pos( 2);
    ha = y_top - y_bot;
    pos = get( H.Fig,'position');
    x_right = pos( 3) - 5;
  
    wa = x_right - x_left;
    H.Ax_cbar = axes('parent',H.Fig,'units','pixels','position',[x_left,y_bot,wa,ha],...
      'fontsize',get( H.Ax{ 1,1},'fontsize'));
    colormap( H.Ax_cbar, colormap( H.Ax{ 1,5}));
    H.Ax_cbar.XAxis.Visible = 'off';
    H.Ax_cbar.YAxis.Visible = 'off';
    set( H.Ax_cbar,'clim',clim_err);
    H.cbar = colorbar( H.Ax_cbar,'location','west');
  
    ticklabels = get( H.cbar,'ticklabels');
    for i = 1: length( ticklabels)
      ticklabels{ i} = ['10^{' ticklabels{i} '}'];
    end
    set( H.cbar,'ticklabels',ticklabels);

  end

  %% Convert maximum error back from logarithm
  err_max = 10^err_max;

end

function [H, err_max] = analyse_mapping_derivatives_tests_mesh_function_ddxy( filename, ddxy)

  H = [];

  % Read the mesh
  mesh = read_mesh_from_file( filename);
  
  %% Set up the plot
  if do_print

    wa = 150;
    ha = 150;
    margins_hor = [50, 50, 5, 5, 50, 5, 5, 150];
    margins_ver = [75, 5, 5, 5];
    H = setup_multipanel_figure( wa, ha, margins_hor, margins_ver);
    
    for i = 1: size( H.Ax,1)
      for j = 1: size( H.Ax,2)
        set( H.Ax{ i,j},'xtick',[],'ytick',[],'clim',[-1,1]);
      end
    end
  
    for i = 1: 1
      for j = 1: size( H.Ax,2)
        set( H.Ax{ i,j},'xaxislocation','top')
      end
    end
  
    title( H.Ax{ 1,1},'Exact solution')
    xlabel( H.Ax{ 1,1},'-')
    title( H.Ax{ 1,3},'Discretised approximation')
    title( H.Ax{ 1,6},'Discretisation error')
  
    ylabel( H.Ax{ 1,1},'A')
    ylabel( H.Ax{ 2,1},'B')
    ylabel( H.Ax{ 3,1},'C')
  
    xlabel( H.Ax{ 1,2},'from A')
    xlabel( H.Ax{ 1,3},'from B')
    xlabel( H.Ax{ 1,4},'from C')
  
    ylabel( H.Ax{ 1,2},'to A')
    ylabel( H.Ax{ 2,2},'to B')
    ylabel( H.Ax{ 3,2},'to C')
  
    xlabel( H.Ax{ 1,5},'from A')
    xlabel( H.Ax{ 1,6},'from B')
    xlabel( H.Ax{ 1,7},'from C')
  
    ylabel( H.Ax{ 1,5},'to A')
    ylabel( H.Ax{ 2,5},'to B')
    ylabel( H.Ax{ 3,5},'to C')
  
    for i = 1:3
      for j = 5:7
        colormap( H.Ax{ i,j},turbo(256));
      end
    end

  end
  
  %% Read data
  d_a_ex = ncread( filename,[ddxy '_a_ex']);
  d_b_ex = ncread( filename,[ddxy '_b_ex']);
  d_c_ex = ncread( filename,[ddxy '_c_ex']);

  d_a_a  = ncread( filename,[ddxy '_a_a']);
  d_a_b  = ncread( filename,[ddxy '_a_b']);
  d_a_c  = ncread( filename,[ddxy '_a_c']);

  d_b_a  = ncread( filename,[ddxy '_b_a']);
  d_b_b  = ncread( filename,[ddxy '_b_b']);
  d_b_c  = ncread( filename,[ddxy '_b_c']);

  d_c_a  = ncread( filename,[ddxy '_c_a']);
  d_c_b  = ncread( filename,[ddxy '_c_b']);
  d_c_c  = ncread( filename,[ddxy '_c_c']);

  %% Calculate errors
  err_a_a = calc_log10_discretisation_error( d_a_ex, d_a_a);
  err_a_b = calc_log10_discretisation_error( d_b_ex, d_a_b);
  err_a_c = calc_log10_discretisation_error( d_c_ex, d_a_c);

  err_b_a = calc_log10_discretisation_error( d_a_ex, d_b_a);
  err_b_b = calc_log10_discretisation_error( d_b_ex, d_b_b);
  err_b_c = calc_log10_discretisation_error( d_c_ex, d_b_c);

  err_c_a = calc_log10_discretisation_error( d_a_ex, d_c_a);
  err_c_b = calc_log10_discretisation_error( d_b_ex, d_c_b);
  err_c_c = calc_log10_discretisation_error( d_c_ex, d_c_c);

  %% Calculate ranges

  % Calculate range for data
  clim = calc_color_limits( d_a_ex);

  % Calculate range for errors
  err_min = min([ ...
    min( err_a_a( err_a_a > -Inf & mesh.VBI   == 0))
    min( err_a_b( err_a_b > -Inf & mesh.TriBI == 0))
    min( err_a_c( err_a_c > -Inf & mesh.EBI   == 0))
    min( err_b_a( err_b_a > -Inf & mesh.VBI   == 0))
    min( err_b_b( err_b_b > -Inf & mesh.TriBI == 0))
    min( err_b_c( err_b_c > -Inf & mesh.EBI   == 0))
    min( err_c_a( err_c_a > -Inf & mesh.VBI   == 0))
    min( err_c_b( err_c_b > -Inf & mesh.TriBI == 0))
    min( err_c_c( err_c_c > -Inf & mesh.EBI   == 0))
    ]);
  err_max = max([ ...
    max( err_a_a( err_a_a > -Inf & mesh.VBI   == 0))
    max( err_a_b( err_a_b > -Inf & mesh.TriBI == 0))
    max( err_a_c( err_a_c > -Inf & mesh.EBI   == 0))
    max( err_b_a( err_b_a > -Inf & mesh.VBI   == 0))
    max( err_b_b( err_b_b > -Inf & mesh.TriBI == 0))
    max( err_b_c( err_b_c > -Inf & mesh.EBI   == 0))
    max( err_c_a( err_c_a > -Inf & mesh.VBI   == 0))
    max( err_c_b( err_c_b > -Inf & mesh.TriBI == 0))
    max( err_c_c( err_c_c > -Inf & mesh.EBI   == 0))
    ]);
  clim_err = [err_min, err_max];

  %% Plot data and errors
  
  if do_print

    plot_mesh_data_a( H.Ax{ 1,1}, mesh, d_a_ex, clim);
    plot_mesh_data_a( H.Ax{ 1,2}, mesh, d_a_a, clim);
    plot_mesh_data_a( H.Ax{ 1,3}, mesh, d_b_a, clim);
    plot_mesh_data_a( H.Ax{ 1,4}, mesh, d_c_a, clim);
    plot_mesh_data_a( H.Ax{ 1,5}, mesh, err_a_a, clim_err);
    plot_mesh_data_a( H.Ax{ 1,6}, mesh, err_b_a, clim_err);
    plot_mesh_data_a( H.Ax{ 1,7}, mesh, err_c_a, clim_err);
    
    plot_mesh_data_b( H.Ax{ 2,1}, mesh, d_b_ex, clim);
    plot_mesh_data_b( H.Ax{ 2,2}, mesh, d_a_b, clim);
    plot_mesh_data_b( H.Ax{ 2,3}, mesh, d_b_b, clim);
    plot_mesh_data_b( H.Ax{ 2,4}, mesh, d_c_b, clim);
    plot_mesh_data_b( H.Ax{ 2,5}, mesh, err_a_b, clim_err);
    plot_mesh_data_b( H.Ax{ 2,6}, mesh, err_b_b, clim_err);
    plot_mesh_data_b( H.Ax{ 2,7}, mesh, err_c_b, clim_err);
  
    plot_mesh_data_c( H.Ax{ 3,1}, mesh, d_c_ex, clim, colormap( H.Ax{ 3,1}));
    plot_mesh_data_c( H.Ax{ 3,2}, mesh, d_a_c, clim, colormap( H.Ax{ 3,2}));
    plot_mesh_data_c( H.Ax{ 3,3}, mesh, d_b_c, clim, colormap( H.Ax{ 3,3}));
    plot_mesh_data_c( H.Ax{ 3,4}, mesh, d_c_c, clim, colormap( H.Ax{ 3,4}));
    plot_mesh_data_c( H.Ax{ 3,5}, mesh, err_a_c, clim_err, colormap( H.Ax{ 3,5}));
    plot_mesh_data_c( H.Ax{ 3,6}, mesh, err_b_c, clim_err, colormap( H.Ax{ 3,6}));
    plot_mesh_data_c( H.Ax{ 3,7}, mesh, err_c_c, clim_err, colormap( H.Ax{ 3,7}));
  
    %% Add colorbar for errors
    pos = get( H.Ax{ 1,7},'position');
    x_left  = pos( 1) + pos( 3) + 5;
    y_top   = pos( 2) + pos( 4);
    pos = get( H.Ax{ 3,7},'position');
    y_bot   = pos( 2);
    ha = y_top - y_bot;
    pos = get( H.Fig,'position');
    x_right = pos( 3) - 5;
  
    wa = x_right - x_left;
    H.Ax_cbar = axes('parent',H.Fig,'units','pixels','position',[x_left,y_bot,wa,ha],...
      'fontsize',get( H.Ax{ 1,1},'fontsize'));
    colormap( H.Ax_cbar, colormap( H.Ax{ 1,5}));
    H.Ax_cbar.XAxis.Visible = 'off';
    H.Ax_cbar.YAxis.Visible = 'off';
    set( H.Ax_cbar,'clim',clim_err);
    H.cbar = colorbar( H.Ax_cbar,'location','west');
  
    ticklabels = get( H.cbar,'ticklabels');
    for i = 1: length( ticklabels)
      ticklabels{ i} = ['10^{' ticklabels{i} '}'];
    end
    set( H.cbar,'ticklabels',ticklabels);

  end

  %% Convert maximum error back from logarithm
  err_max = 10^err_max;

end

function [H, err_ddx_max, err_ddy_max, err_d2dx2_max, err_d2dxdy_max, err_d2dy2_max] = ...
    analyse_mapping_derivatives_tests_mesh_function_2nd( filename)

  H = [];

  % Read the mesh
  mesh = read_mesh_from_file( filename);
  
  %% Set up the plot
  if do_print

    wa = 150;
    ha = 150;
    margins_hor = [50, 5, 5, 5, 5, 5];
    margins_ver = [75, 5, 5, 75];
    H = setup_multipanel_figure( wa, ha, margins_hor, margins_ver);
    
    for i = 1: size( H.Ax,1)
      for j = 1: size( H.Ax,2)
        set( H.Ax{ i,j},'xtick',[],'ytick',[],'clim',[-1,1]);
      end
    end
  
    for i = 1: 1
      for j = 1: size( H.Ax,2)
        set( H.Ax{ i,j},'xaxislocation','top')
      end
    end
  
    ylabel( H.Ax{ 1,1},'Exact')
    ylabel( H.Ax{ 2,1},'Disc.')
    ylabel( H.Ax{ 3,1},'Err.')
  
    xlabel( H.Ax{ 1,1},'d/dx')
    xlabel( H.Ax{ 1,2},'d/dy')
    xlabel( H.Ax{ 1,3},'d^2/dx^2')
    xlabel( H.Ax{ 1,4},'d^2/dxdy')
    xlabel( H.Ax{ 1,5},'d^2/dy^2')
  
    for i = 3
      for j = 1:5
        colormap( H.Ax{ i,j},turbo(256));
      end
    end

  end
  
  %% Read data
  ddx_b_ex    = ncread( filename,'ddx_b_ex');
  ddy_b_ex    = ncread( filename,'ddy_b_ex');
  d2dx2_b_ex  = ncread( filename,'d2dx2_b_ex');
  d2dxdy_b_ex = ncread( filename,'d2dxdy_b_ex');
  d2dy2_b_ex  = ncread( filename,'d2dy2_b_ex');

  ddx_b_b_2nd    = ncread( filename,'ddx_b_b_2nd');
  ddy_b_b_2nd    = ncread( filename,'ddy_b_b_2nd');
  d2dx2_b_b_2nd  = ncread( filename,'d2dx2_b_b_2nd');
  d2dxdy_b_b_2nd = ncread( filename,'d2dxdy_b_b_2nd');
  d2dy2_b_b_2nd  = ncread( filename,'d2dy2_b_b_2nd');

  %% Calculate errors
  err_ddx_b_b_2nd    = calc_log10_discretisation_error( ddx_b_ex   , ddx_b_b_2nd);
  err_ddy_b_b_2nd    = calc_log10_discretisation_error( ddy_b_ex   , ddy_b_b_2nd);
  err_d2dx2_b_b_2nd  = calc_log10_discretisation_error( d2dx2_b_ex , d2dx2_b_b_2nd);
  err_d2dxdy_b_b_2nd = calc_log10_discretisation_error( d2dxdy_b_ex, d2dxdy_b_b_2nd);
  err_d2dy2_b_b_2nd  = calc_log10_discretisation_error( d2dy2_b_ex , d2dy2_b_b_2nd);

  %% Calculate ranges

  % Calculate range for data
  clim_ddx    = calc_color_limits( ddx_b_ex);
  clim_ddy    = calc_color_limits( ddy_b_ex);
  clim_d2dx2  = calc_color_limits( d2dx2_b_ex);
  clim_d2dxdy = calc_color_limits( d2dxdy_b_ex);
  clim_d2dy2  = calc_color_limits( d2dy2_b_ex);

  % Calculate range for errors
  clim_ddx_err    = calc_color_limits( err_ddx_b_b_2nd(    err_ddx_b_b_2nd    > -Inf & mesh.TriBI == 0));
  clim_ddy_err    = calc_color_limits( err_ddy_b_b_2nd(    err_ddy_b_b_2nd    > -Inf & mesh.TriBI == 0));
  clim_d2dx2_err  = calc_color_limits( err_d2dx2_b_b_2nd(  err_d2dx2_b_b_2nd  > -Inf & mesh.TriBI == 0));
  clim_d2dxdy_err = calc_color_limits( err_d2dxdy_b_b_2nd( err_d2dxdy_b_b_2nd > -Inf & mesh.TriBI == 0));
  clim_d2dy2_err  = calc_color_limits( err_d2dy2_b_b_2nd(  err_d2dy2_b_b_2nd  > -Inf & mesh.TriBI == 0));

  err_ddx_max    = clim_ddx_err(    2);
  err_ddy_max    = clim_ddy_err(    2);
  err_d2dx2_max  = clim_d2dx2_err(  2);
  err_d2dxdy_max = clim_d2dxdy_err( 2);
  err_d2dy2_max  = clim_d2dy2_err(  2);

  %% Plot data and errors
  
  if do_print

    plot_mesh_data_b( H.Ax{ 1,1}, mesh, ddx_b_ex   , clim_ddx);
    plot_mesh_data_b( H.Ax{ 1,2}, mesh, ddy_b_ex   , clim_ddy);
    plot_mesh_data_b( H.Ax{ 1,3}, mesh, d2dx2_b_ex , clim_d2dx2);
    plot_mesh_data_b( H.Ax{ 1,4}, mesh, d2dxdy_b_ex, clim_d2dxdy);
    plot_mesh_data_b( H.Ax{ 1,5}, mesh, d2dy2_b_ex , clim_d2dy2);
    
    plot_mesh_data_b( H.Ax{ 2,1}, mesh, ddx_b_b_2nd   , clim_ddx);
    plot_mesh_data_b( H.Ax{ 2,2}, mesh, ddy_b_b_2nd   , clim_ddy);
    plot_mesh_data_b( H.Ax{ 2,3}, mesh, d2dx2_b_b_2nd , clim_d2dx2);
    plot_mesh_data_b( H.Ax{ 2,4}, mesh, d2dxdy_b_b_2nd, clim_d2dxdy);
    plot_mesh_data_b( H.Ax{ 2,5}, mesh, d2dy2_b_b_2nd , clim_d2dy2);
    
    plot_mesh_data_b( H.Ax{ 3,1}, mesh, err_ddx_b_b_2nd   , clim_ddx_err);
    plot_mesh_data_b( H.Ax{ 3,2}, mesh, err_ddy_b_b_2nd   , clim_ddy_err);
    plot_mesh_data_b( H.Ax{ 3,3}, mesh, err_d2dx2_b_b_2nd , clim_d2dx2_err);
    plot_mesh_data_b( H.Ax{ 3,4}, mesh, err_d2dxdy_b_b_2nd, clim_d2dxdy_err);
    plot_mesh_data_b( H.Ax{ 3,5}, mesh, err_d2dy2_b_b_2nd , clim_d2dy2_err);
  
    %% Add colorbars for errors
  
    for j = 1: 5
  
      pos = get( H.Ax{ 3,j},'position');
      x_left  = pos( 1);
      x_right = pos( 1) + pos( 3);
      y_top   = pos( 2) - 5;
      y_bot   = 5;
      wa = x_right - x_left;
      ha = y_top   - y_bot;
    
      H.Ax_cbar{j} = axes('parent',H.Fig,'units','pixels','position',[x_left,y_bot,wa,ha],...
        'fontsize',get( H.Ax{ 3,j},'fontsize'));
      colormap( H.Ax_cbar{j}, colormap( H.Ax{ 3,1}));
      H.Ax_cbar{j}.XAxis.Visible = 'off';
      H.Ax_cbar{j}.YAxis.Visible = 'off';
      set( H.Ax_cbar{j},'clim',get(H.Ax{ 3,j},'clim'));
      H.cbar{ j} = colorbar( H.Ax_cbar{j},'location','north');
    
      ticklabels = get( H.cbar{ j},'ticklabels');
      for ii = 1: length( ticklabels)
        ticklabels{ ii} = ['10^{' ticklabels{ii} '}'];
      end
      set( H.cbar{ j},'ticklabels',ticklabels);
  
    end

  end

  %% Convert maximum error back from logarithm
  err_ddx_max    = 10^err_ddx_max;
  err_ddy_max    = 10^err_ddy_max;
  err_d2dx2_max  = 10^err_d2dx2_max;
  err_d2dxdy_max = 10^err_d2dxdy_max;
  err_d2dy2_max  = 10^err_d2dy2_max;

end

function clim = calc_color_limits( d)

  d_min = min( d);
  d_max = max( d);
  if d_min == d_max
    if d_min == 0
      d_min = min( d);
      d_max = max( d);
      if d_min == d_max
        d_min = -1e-15;
        d_max = 1e-15;
      end
    else
      d = 1e-5;
      d_min = d_min * (1 - d);
      d_max = d_max + (d_max - d_min);
    end
  end

  clim = [d_min, d_max];

end

function err = calc_log10_discretisation_error( d_ex, d_disc)
  if (max( d_ex) == 0)
    err = log( abs( d_disc - d_ex)) / log( 10);
  else
    err = log( abs( d_disc - d_ex) ./ max( abs( d_ex))) / log( 10);
  end
end

function plot_mesh_data_a( ax, mesh, d, clim)
  patch('parent',ax,'vertices',mesh.V(1:mesh.nV,:),'faces',mesh.Tri(1:mesh.nTri,:),...
    'facecolor','interp','facevertexcdata',d,'edgecolor','none');
  set(ax,'xlim',[mesh.xmin,mesh.xmax],'ylim',[mesh.ymin,mesh.ymax],'clim',clim);
end
function plot_mesh_data_b( ax, mesh, d, clim)
  patch('parent',ax,'vertices',mesh.V(1:mesh.nV,:),'faces',mesh.Tri(1:mesh.nTri,:),...
    'facecolor','flat','facevertexcdata',d,'edgecolor','none');
  set(ax,'xlim',[mesh.xmin,mesh.xmax],'ylim',[mesh.ymin,mesh.ymax],'clim',clim);
end
function plot_mesh_data_c( ax, mesh, d, clim, cmap)

  set(ax,'xlim',[mesh.xmin,mesh.xmax],'ylim',[mesh.ymin,mesh.ymax]);
  
  ncols = size( cmap,1);
  set(ax,'clim',clim)
  
  for ci = 1: ncols
    linedata(ci).x = [];
    linedata(ci).y = [];
  end
  for aci = 1: mesh.nE
    ci = (d( aci) - clim(1)) / (clim(2) - clim(1));
    ci = max(1,min(ncols,1 + round(ci*(ncols-1))));
    linedata(ci).x( end+1) = mesh.E( aci,1);
    linedata(ci).y( end+1) = mesh.E( aci,2);
  end
  for ci = 1: ncols
    line('parent',ax,'xdata',linedata( ci).x,'ydata',linedata( ci).y,'linestyle','none',...
      'marker','o','markerfacecolor',cmap( ci,:),'markeredgecolor',cmap( ci,:),'markersize',8);
  end
  
  % lastly NaN values
  line('parent',ax,'xdata',mesh.E( isnan(d),1),'ydata',mesh.E( isnan(d),2),'linestyle','none',...
      'marker','x','markerfacecolor','r','markeredgecolor','r','markersize',8);

end

end