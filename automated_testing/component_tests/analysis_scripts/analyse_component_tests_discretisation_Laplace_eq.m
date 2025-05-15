function analyse_component_tests_discretisation_Laplace_eq( foldername_automated_testing, do_print_figures)
% Analyse the results of all the Laplace eq solving discretisation
% component tests

disp('    Analysing Laplace equation solving component tests...')
disp('')

foldername_results = [foldername_automated_testing '/component_tests/results/discretisation/solve_Laplace_eq'];
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
  analyse_component_test_discretisation_Laplace_eq( filenames( fi).name)
end

disp('')

function analyse_component_test_discretisation_Laplace_eq( filename_short)
  % Analyse the results of the Laplace equation solving component
  % test for a single mesh

  disp(['      ' filename_short '...'])

  filename_full = [foldername_results '/' filename_short];
  
  [H, test_results] = analyse_Laplace_eq_test(  filename_full);
  
  write_to_scoreboard_file( filename_short, test_results);

  if do_print_figures
    filename_png = strrep( filename_short, '.nc', '.png');
    print( H.Fig, [foldername_figures '/' filename_png], '-dpng');
    close( H.Fig)
  end
  
end
  
function write_to_scoreboard_file( filename_short, test_results)

  % Set up a scoreboard results structure
  test_name = filename_short(5:end-3);
  single_run = initialise_single_test_run( test_name, ...
    'component_tests/discretisation/Laplace_eq_solving');

  single_run = add_cost_function_to_single_run( single_run, 'max_abs_err', 'max( abs( f - f_ex))'        , test_results.max_abs_err);
  single_run = add_cost_function_to_single_run( single_run, 'rmse'       , 'sqrt( mean( (f - f_ex).^2 ))', test_results.rmse);
 
  % Write to scoreboard file
  write_scoreboard_file( foldername_automated_testing, single_run)

end

function [H, test_results] = analyse_Laplace_eq_test( filename)

  H = [];

  % Read the mesh
  mesh = read_mesh_from_file( filename);
  
  %% Set up the plot
  if do_print_figures

    wa = 250;
    ha = 250;
    margins_hor = [150, 5, 25, 150];
    margins_ver = [75, 35];
    H = setup_multipanel_figure( wa, ha, margins_hor, margins_ver);
    
    for i = 1: size( H.Ax,1)
      for j = 1: size( H.Ax,2)
        set( H.Ax{ i,j},'xtick',[],'ytick',[]);
      end
    end
  
    for i = 1: 1
      for j = 1: size( H.Ax,2)
        set( H.Ax{ i,j},'xaxislocation','top')
      end
    end
  
    title( H.Ax{ 1,1},'Exact solution')
    title( H.Ax{ 1,2},'Discretised approximation')
    title( H.Ax{ 1,3},'Discretisation error')
  
    for i = 1:1
      for j = 1:3
        colormap( H.Ax{ i,j},turbo(256));
      end
    end
  
    % Add colorbars
    pos = get( H.Ax{1,1},'position');
    H.Cbar1 = colorbar( H.Ax{1,1},'location','westoutside');
    set( H.Ax{1,1},'position',pos);
  
    pos = get( H.Ax{1,3},'position');
    H.Cbar2 = colorbar( H.Ax{1,3},'location','eastoutside');
    set( H.Ax{1,3},'position',pos);

  end
  
  %% Read data
  f_ex = ncread( filename,'f_ex');
  f    = ncread( filename,'f');

  %% Calculate errors
  err = f - f_ex;

  % Calculate integrated values for output
  test_results.max_abs_err = max( abs( err));
  test_results.rmse = sqrt( mean( err.^2));

  %% Plot data and errors
  
  if do_print_figures

    plot_mesh_data_b( H.Ax{ 1,1}, mesh, f_ex, [min(f_ex), max(f_ex)]);
    plot_mesh_data_b( H.Ax{ 1,2}, mesh, f   , [min(f_ex), max(f_ex)]);
    plot_mesh_data_b( H.Ax{ 1,3}, mesh, err , [min(err ), max(err )]);

  end

end

function plot_mesh_data_b( ax, mesh, d, clim)
  patch('parent',ax,'vertices',mesh.V(1:mesh.nV,:),'faces',mesh.Tri(1:mesh.nTri,:),...
    'facecolor','flat','facevertexcdata',d,'edgecolor','none');
  set(ax,'xlim',[mesh.xmin,mesh.xmax],'ylim',[mesh.ymin,mesh.ymax],'clim',clim);
end

end