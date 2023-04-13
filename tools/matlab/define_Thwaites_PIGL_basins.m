clc
clear all
close all

% Plot range
xmid = -1300e3;
ymid =  -450e3;
w    = 350e3;

xmin = xmid - w;
xmax = xmid + w;
ymin = ymid - w;
ymax = ymid + w;

%% Set up GUI
wa = 700;
ha = 700;

margins_hor = [125,125];
margins_ver = [125,50];

wf = margins_hor( 1) + wa + margins_hor( 2);
hf = margins_ver( 1) + ha + margins_ver( 2);

H.Fig = figure( 'position',[100,100,wf,hf],'color','w');
H.Ax  = axes('parent',H.Fig,'units','pixels','position',[margins_hor(1),margins_ver(1),wa,ha],...
  'xlim',[xmin,xmax],'ylim',[ymin,ymax],'fontsize',24,'xgrid','on','ygrid','on');

%% Plot velocity

filename = '/Users/berends/Documents/Datasets/Ant_velocity_Rignot/Ant_velocity_Rignot2011_10km.nc';

x = ncread( filename,'x');
y = ncread( filename,'y');
u = ncread( filename,'u_surf');
v = ncread( filename,'v_surf');
uabs = sqrt( u.^2 + v.^2);

image('parent',H.Ax,'xdata',x,'ydata',y,'cdata',uabs','cdatamapping','scaled');
set( H.Ax,'colorscale','log');
colormap( H.Ax, sentinelmap(256));

%% Plot PIGL basin
V   = [
  -1.64e6, -3.4e5;
  -1.60e6, -3.5e5;
  -1.55e6, -3.4e5;
  -1.50e6, -3.2e5;
  -1.45e6, -2.9e5;
  -1.40e6, -2.5e5;
  -1.37e6, -2.0e5;
  -1.34e6, -1.7e5;
  -1.30e6, -1.6e5;
  -1.26e6, -1.6e5;
  -1.22e6, -1.7e5;
  -1.18e6, -1.75e5;
  -1.14e6, -1.75e5;
  -1.11e6, -1.72e5;
  -1.09e6, -1.6e5;
  -1.085e6, -1.4e5;
  -1.09e6, -1.2e5;
  -1.1e6, -1.0e5;
  -1.13e6, -0.7e5;
  -1.17e6, -0.4e5;
  -1.21e6, -0.2e5;
  -1.26e6, -0.0e5;
  -1.32e6, 0.1e5;
  -1.45e6, 0.1e5;
  -1.48e6, 0.15e5;
  -1.51e6, 0.35e5;
  -1.53e6, 0.75e5;
  -1.55e6, 0.95e5;
  -1.58e6, 0.1e6;
  -1.62e6, 0.11e6;
  -1.65e6, 0.12e6;
  -1.67e6, 0.10e6;
  -1.69e6, 0.9e5;
  -1.71e6, 0.5e5;
  -1.74e6, 0.1e5;
  -1.75e6, -0.5e5;
  -1.75e6, -0.15e6;
  -1.71e6, -0.19e6;
  -1.66e6, -0.2e6;
  -1.64e6, -0.21e6;
  -1.63e6, -0.23e6;
  -1.63e6, -0.29e6
  ];
Tri = 1:length(V);
patch('parent',H.Ax,'vertices',V,'faces',Tri,'facecolor','none','edgecolor','k','linewidth',3,'marker','o');


%% Plot Thwaites basin
V   = [
  -1.6e6, -5.4e5;
  -1.55e6, -5.4e5;
  -1.50e6, -5.5e5;
  -1.45e6, -5.6e5;
  -1.40e6, -5.65e5;
  -1.37e6, -5.75e5;
  -1.35e6, -6e5;
  -1.35e6, -6.5e5;
  -1.34e6, -6.9e5;
  -1.32e6, -7.3e5;
  -1.29e6, -7.6e5;
  -1.25e6, -7.8e5;
  -1.22e6, -7.8e5;
  -1.20e6, -7.6e5;
  -1.18e6, -7.4e5;
  -1.15e6, -6.9e5;
  -1.14e6, -6.4e5;
  -1.14e6, -5.9e5;
  -1.11e6, -5.6e5;
  -1.08e6, -5.5e5;
  -1.04e6, -5.4e5;
  -1.01e6, -5.2e5;
  -0.99e6, -5.0e5;
  -0.99e6, -4.6e5;
  -1.02e6, -4.4e5;
  -1.04e6, -4.2e5;
  -1.06e6, -3.9e5;
  -1.07e6, -3.5e5;
  -1.07e6, -3.2e5;
  -1.09e6, -2.8e5;
  -1.12e6, -2.5e5;
  -1.15e6, -2.2e5;
  -1.18e6, -1.9e5;
  -1.22e6, -1.7e5;
  -1.26e6, -1.6e5;
  -1.30e6, -1.6e5;
  -1.34e6, -1.7e5;
  -1.37e6, -2.0e5;
  -1.40e6, -2.5e5;
  -1.45e6, -2.9e5;
  -1.50e6, -3.2e5;
  -1.55e6, -3.4e5;
  -1.60e6, -3.5e5;
  -1.64e6, -3.4e5;
  ];
Tri = 1:length(V);
patch('parent',H.Ax,'vertices',V,'faces',Tri,'facecolor','none','edgecolor','b','linewidth',3,'marker','o');