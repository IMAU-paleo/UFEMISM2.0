clc
clear all
close all

% Read Pattyn et al. (2008) ensemble results
foldername_ensemble = '/Users/Beren017/Documents/GitHub/data/model_ensembles/ISMIP-HOM/ismip_all';
[HO,FS] = process_ISMIP_HOM_ensemble_experiment_B( foldername_ensemble);

%% Set up GUI

close all

c_HO      = [0.1 ,0.5 ,0.2];
c_FS      = [0.2 ,0.5 ,1.0];

colors = crameri('hawaii',3);

c_SIASSA = colors(1,:);
c_DIVA   = colors(2,:);
c_BPA    = colors(3,:);

wa = 250;
ha = 200;

margins_hor = [90,50,50,25];
margins_ver = [50,50,90];

H = setup_multipanel_figure( wa, ha, margins_hor, margins_ver);

set( H.Ax{ 1,1},'xtick',0:0.2:1,'xticklabels','');
set( H.Ax{ 1,2},'xtick',0:0.2:1,'xticklabels','');
set( H.Ax{ 1,3},'xtick',0:0.2:1,'xticklabels','');
set( H.Ax{ 2,1},'xtick',0:0.2:1);
set( H.Ax{ 2,2},'xtick',0:0.2:1);
set( H.Ax{ 2,3},'xtick',0:0.2:1);

xlabel( H.Ax{ 2,1},'x / L');
xlabel( H.Ax{ 2,2},'x / L');
xlabel( H.Ax{ 2,3},'x / L');

ylabel( H.Ax{ 1,1},'Surface x-velocity (m/yr)');
ylabel( H.Ax{ 2,1},'Surface x-velocity (m/yr)');

title( H.Ax{ 1,1},'5 km');
title( H.Ax{ 1,2},'10 km');
title( H.Ax{ 1,3},'20 km');
title( H.Ax{ 2,1},'40 km');
title( H.Ax{ 2,2},'80 km');
title( H.Ax{ 2,3},'160 km');

% Empty objects for legend
patch( 'parent',H.Ax{ 2,3},'vertices',[],'faces',[],'facecolor',c_HO,...
  'edgecolor','none','facealpha',0.3);
patch( 'parent',H.Ax{ 2,3},'vertices',[],'faces',[],'facecolor',c_FS,...
  'edgecolor','none','facealpha',0.7);
% line( 'parent',H.Ax{ 2,3},'xdata',[],'ydata',[],'color',c_HO,'linewidth',2);
% line( 'parent',H.Ax{ 2,3},'xdata',[],'ydata',[],'color',c_FS,'linewidth',2);
line( 'parent',H.Ax{ 2,3},'xdata',[],'ydata',[],'linewidth',2,'color',c_SIASSA);
line( 'parent',H.Ax{ 2,3},'xdata',[],'ydata',[],'linewidth',2,'color',c_DIVA);
line( 'parent',H.Ax{ 2,3},'xdata',[],'ydata',[],'linewidth',2,'color',c_BPA);

set( H.Ax{ 1,1},'ylim',[0,30]);
set( H.Ax{ 1,2},'ylim',[0,40]);
set( H.Ax{ 1,3},'ylim',[0,60]);
set( H.Ax{ 2,1},'ylim',[0,80]);
set( H.Ax{ 2,2},'ylim',[0,120]);
set( H.Ax{ 2,3},'ylim',[0,120]);

for Li = 1: 6

  if Li==1
    ex = 'L005';
    ax = H.Ax{ 1,1};
  elseif Li == 2
    ex = 'L010';
    ax = H.Ax{ 1,2};
  elseif Li == 3
    ex = 'L020';
    ax = H.Ax{ 1,3};
  elseif Li == 4
    ex = 'L040';
    ax = H.Ax{ 2,1};
  elseif Li == 5
    ex = 'L080';
    ax = H.Ax{ 2,2};
  elseif Li == 6
    ex = 'L160';
    ax = H.Ax{ 2,3};
  end

  % Ensemble ranges
  xdata = [HO.(ex).x    ; flipud( HO.(ex).x    )];
  ydata = [HO.(ex).u_min; flipud( HO.(ex).u_max)];
  patch('parent',ax,'xdata',xdata,'ydata',ydata,'facecolor',c_HO,...
    'edgecolor','none','facealpha',0.3);

  xdata = [FS.(ex).x    ; flipud( FS.(ex).x    )];
  ydata = [FS.(ex).u_min; flipud( FS.(ex).u_max)];
  patch('parent',ax,'xdata',xdata,'ydata',ydata,'facecolor',c_FS,...
    'edgecolor','none','facealpha',0.7);

  % % Ensemble means
  % line('parent',ax,'xdata',HO.(ex).x,'ydata',HO.(ex).u_av,...
  %   'color',c_HO,'linewidth',2)
  % line('parent',ax,'xdata',FS.(ex).x,'ydata',FS.(ex).u_av,...
  %   'color',c_FS,'linewidth',2)
end

%% Read and plot UFEMISM results


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

  if     L == 5
    ax = H.Ax{ 1,1};
  elseif L == 10
    ax = H.Ax{ 1,2};
  elseif L == 20
    ax = H.Ax{ 1,3};
  elseif L == 40
    ax = H.Ax{ 2,1};
  elseif L == 80
    ax = H.Ax{ 2,2};
  elseif L == 160
    ax = H.Ax{ 2,3};
  end

  for approxi = 1:3
    approx = approxs{ approxi};

    if strcmpi( approx,'SIASSA')
      color = c_SIASSA;
    elseif strcmpi( approx,'DIVA')
      color = c_DIVA;
    elseif strcmpi( approx,'BPA')
      color = c_BPA;
    end

    foldername = ['results_ISMIP_HOM_B_' num2str(L) '_' approx];
    filename = [foldername '/transect_ISMIP-HOM.nc'];
    V = ncread( filename,'V');
    xt = fliplr(linspace(0,1,size(V,1)));
    u_3D = ncread( filename,'u_par');
    u_surf = u_3D(:,1);

    % Plot results
    line('parent',ax,'xdata',xt,'ydata',u_surf,'color',color,...
      'linewidth',2)

  end
end

%% Legend

legend( H.Ax{ 2,3},'Full-Stokes','Higher-Order','SIA/SSA','DIVA','BPA','location','northwest')
