function analyse_component_tests_mass_conservation( foldername_automated_testing, do_print_figures)
% Analyse the results of all the mass conservation component tests

disp('    Analysing mass conservation component tests...')
disp('')

foldername_results = [foldername_automated_testing '/component_tests/results/mass_conservation'];
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
  % Analyse the results of the complete set of mesh-to-grid remapping component
  % tests for a single mesh-grid combination

  disp(['      ' filename_short '...']);

  filename_full = [foldername_results '/' filename_short];

  % Read test results
  mesh                = read_mesh_from_file( filename_full);
  Hi                  = ncread( filename_full, 'Hi');
  dHi_dt_ex           = ncread( filename_full, 'dHi_dt_ex');
  dHi_dt_explicit     = ncread( filename_full, 'dHi_dt_explicit');
  dHi_dt_semiimplicit = ncread( filename_full, 'dHi_dt_semiimplicit');
  dHi_dt_implicit     = ncread( filename_full, 'dHi_dt_implicit');
  dHi_dt_overimplicit = ncread( filename_full, 'dHi_dt_overimplicit');

  m = Hi > 100;
  m( mesh.VBI > 0) = false;

  % Plot error
  if do_print_figures

    %% Set up figure

    if contains( filename_full,'linear')
      clim_abs = -1 + [-1,1]*1e-3;
      clim_err = [-0.2,0.2];
    elseif contains( filename_full,'periodic')
      clim_abs = [-10,10];
      clim_err = [-2,2];
    elseif contains( filename_full,'Halfar')
      clim_abs = [-4,4];
      clim_err = [-1,1];
    else
      % error('unknown test function, need to define clim')
    end

    cmap_abs = crameri('bam',33);
    cmap_err = crameri('vik',33);

    wa = 180;
    ha = wa;
    H = setup_multipanel_figure( wa, ha, [100,5,5,5,5,5], [50,5,5]);

    for i = 1: size( H.Ax,1)
      for j = 1: size( H.Ax,2)
        ax = H.Ax{ i,j};
        set( ax,'xtick',[],'ytick',[],'xgrid','off','ygrid','off')
        ax.XAxis.Visible = 'off';
        ax.YAxis.Visible = 'off';
        H.Patch(i,j) = patch( 'parent',ax,'vertices',mesh.Vor,...
          'faces',changem(double(mesh.VVor),NaN),'facecolor','flat',...
          'facevertexcdata',zeros(mesh.nV,1),'edgecolor','none');
        if i==1
          colormap(ax,cmap_abs)
          set(ax,'clim',clim_abs)
        elseif i==2
          colormap(ax,cmap_err)
          set(ax,'clim',clim_err)
        end
        % set(ax,'xlim',[mesh.xmin,mesh.xmax]*0.55,'ylim',[mesh.ymin,mesh.ymax]*0.55)
      end
    end

    title( H.Ax{1,1},'exact')
    title( H.Ax{1,2},'explicit')
    title( H.Ax{1,3},'semi-implicit')
    title( H.Ax{1,4},'implicit')
    title( H.Ax{1,5},'over-implicit')
  
    pos = get( H.Ax{1,1},'position');
    H.Cbar_abs = colorbar(H.Ax{1,1},'location','westoutside');
    set( H.Ax{1,1},'position',pos)
    ylabel( H.Cbar_abs,'dH/dt')
  
    pos = get( H.Ax{2,1},'position');
    H.Cbar_err = colorbar(H.Ax{2,1},'location','westoutside');
    set( H.Ax{2,1},'position',pos)
    ylabel( H.Cbar_err,'err dH/dt')

    %% Plot results

    dHi_dt_ex(           ~m) = dHi_dt_ex( ~m);
    dHi_dt_explicit(     ~m) = dHi_dt_ex( ~m);
    dHi_dt_semiimplicit( ~m) = dHi_dt_ex( ~m);
    dHi_dt_implicit(     ~m) = dHi_dt_ex( ~m);
    dHi_dt_overimplicit( ~m) = dHi_dt_ex( ~m);

    set( H.Patch(1,1),'facevertexcdata',dHi_dt_ex)
    set( H.Patch(1,2),'facevertexcdata',dHi_dt_explicit)
    set( H.Patch(1,3),'facevertexcdata',dHi_dt_semiimplicit)
    set( H.Patch(1,4),'facevertexcdata',dHi_dt_implicit)
    set( H.Patch(1,5),'facevertexcdata',dHi_dt_overimplicit)

    err_ex           = dHi_dt_ex           - dHi_dt_ex;
    err_explicit     = dHi_dt_explicit     - dHi_dt_ex;
    err_semiimplicit = dHi_dt_semiimplicit - dHi_dt_ex;
    err_implicit     = dHi_dt_implicit     - dHi_dt_ex;
    err_overimplicit = dHi_dt_overimplicit - dHi_dt_ex;

    set( H.Patch(2,1),'facevertexcdata',err_ex          )
    set( H.Patch(2,2),'facevertexcdata',err_explicit    )
    set( H.Patch(2,3),'facevertexcdata',err_semiimplicit)
    set( H.Patch(2,4),'facevertexcdata',err_implicit    )
    set( H.Patch(2,5),'facevertexcdata',err_overimplicit)

    %% Save figure

    filename_png = strrep( filename_short, '.nc', '.png');
    print( H.Fig, [foldername_figures '/' filename_png], '-dpng');
    close( H.Fig)

  end

  write_to_scoreboard_file( filename_short, m, dHi_dt_ex, dHi_dt_explicit, ...
    dHi_dt_semiimplicit, dHi_dt_implicit, dHi_dt_overimplicit);
  
end
  
function write_to_scoreboard_file( filename_short, m, dHi_dt_ex, dHi_dt_explicit, ...
    dHi_dt_semiimplicit, dHi_dt_implicit, dHi_dt_overimplicit)

  % Set up a scoreboard results structure
  test_name = filename_short(5:end-3);
  single_run = initialise_single_test_run( test_name, ...
    'component_tests/mass_conservation/mesh_vertices_to_grid');

  % Calculate cost functions
  rmse_explicit     = sqrt( mean( (dHi_dt_explicit(    m) - dHi_dt_ex(m)).^2));
  rmse_semiimplicit = sqrt( mean( (dHi_dt_semiimplicit(m) - dHi_dt_ex(m)).^2));
  rmse_implicit     = sqrt( mean( (dHi_dt_implicit(    m) - dHi_dt_ex(m)).^2));
  rmse_overimplicit = sqrt( mean( (dHi_dt_overimplicit(m) - dHi_dt_ex(m)).^2));

  % Add cost functions to results structure
  single_run = add_cost_function_to_single_run( single_run, 'rmse_explicit'    , 'sqrt( mean( (dHi_dt_explicit - dHi_dt_ex).^2))'    , rmse_explicit);
  single_run = add_cost_function_to_single_run( single_run, 'rmse_semiimplicit', 'sqrt( mean( (dHi_dt_semiimplicit - dHi_dt_ex).^2))', rmse_semiimplicit);
  single_run = add_cost_function_to_single_run( single_run, 'rmse_implicit'    , 'sqrt( mean( (dHi_dt_implicit - dHi_dt_ex).^2))'    , rmse_implicit);
  single_run = add_cost_function_to_single_run( single_run, 'rmse_overimplicit', 'sqrt( mean( (dHi_dt_overimplicit - dHi_dt_ex).^2))', rmse_overimplicit);

  % Write to scoreboard file
  write_scoreboard_file( foldername_automated_testing, single_run);

end
  
end