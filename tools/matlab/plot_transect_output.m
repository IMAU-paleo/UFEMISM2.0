clc
clear all
close all

foldername = ['../../automated_testing/integrated_tests/realistic/Antarctica/' ...
  'initialisation/Ant_init_20kyr_invBMB_invpwf_40km/results/'];

% filename = [foldername '/transect_PineIsland_centralflowline.nc'];
% filename = [foldername '/transect_PineIsland_groundingline.nc'];
filename = [foldername '/transect_Thwaites_centralflowline.nc'];
% filename = [foldername '/transect_Thwaites_groundingline.nc'];

data_to_plot = 'u_par';

% Read transect geometry
transect.V    = ncread( filename,'V');
transect.d    = calc_distance_along_transect_in_km( transect.V);
transect.zeta = ncread( filename,'zeta');

time = ncread( filename,'time');
ti = length( time);

% Read ice geometry
Hi   = ncread( filename,'Hi',   [1,  ti],[Inf,    1]);
Hb   = ncread( filename,'Hb',   [1,  ti],[Inf,    1]);
Hs   = ncread( filename,'Hs',   [1,  ti],[Inf,    1]);
Hib  = ncread( filename,'Hib',  [1,  ti],[Inf,    1]);
u_3D = ncread( filename, data_to_plot, [1,1,ti],[Inf,Inf,1]);

% Calculate transect vertices and faces
[V,F,C] = calc_transect_patch( transect.d, transect.zeta, Hi, Hs, u_3D);

% Plot
H = setup_multipanel_figure( 600, 400, [100,150],[25,100]); H.Ax = H.Ax{1,1};
xlabel( H.Ax,'Distance along transect (km)');
ylabel( H.Ax,'z (m)');
pos = get( H.Ax,'position');
H.Cbar = colorbar( H.Ax,'location','eastoutside');
set( H.Ax,'position',pos);
ylabel( H.Cbar,'Velocity (m/yr)');
colormap( H.Ax, turbo(256));
set( H.Ax,'colorscale','log','clim',[50,5000]);
set( H.Cbar,'ticks',[50,100,200,500,1000,2000,5000]);

H.Patch = patch('parent',H.Ax,'vertices',V,'faces',F,'facevertexcdata',C,...
  'facecolor','interp','edgecolor','none');
H.Hs  = line('parent',H.Ax,'xdata',transect.d,'ydata',Hs, 'linewidth',2);
H.Hb  = line('parent',H.Ax,'xdata',transect.d,'ydata',Hb, 'linewidth',2);
H.Hib = line('parent',H.Ax,'xdata',transect.d,'ydata',Hib,'linewidth',2);

function d = calc_distance_along_transect_in_km( V)
  d = zeros( size( V,1),1);
  for i = 2: size( V,1)
    d(i) = d(i-1) + norm( V(i,:) - V(i-1,:)) / 1e3;
  end
end

function [V,F,C] = calc_transect_patch( d, zeta, Hi, Hs, cc)

  nx = length(d);
  nz = length(zeta);

  nV = nx * nz;
  nF = (nx-1) * (nz-1);

  ik2vi = zeros( nx,nz);
  V = zeros( nV,2);
  F = zeros( nF,4);
  C = zeros( nV,1);

  vi = 0;
  for i = 1: nx
    for k = 1: nz
      vi = vi + 1;
      ik2vi( i,k) = vi;
      x = d( i);
      z = Hs( i) - zeta( k) * Hi( i);
      V( vi,:) = [x,z];
      C( vi  ) = cc( i,k);
    end
  end

  fi = 0;
  for i = 1: nx-1
    for k = 1: nz-1
      fi = fi + 1;
      F( fi,:) = [ik2vi( i,k), ik2vi( i+1,k), ik2vi( i+1,k+1), ik2vi( i,k+1)];
    end
  end

end

