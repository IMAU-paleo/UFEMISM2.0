clc
clear all
close all

filename_src       = '../../exp_II_target_5km/restart_ANT.nc';
filename_src_topo  = 'exp_II_topography_5km.nc';
filename_dst_hi    = 'exp_II_topography_hi_5km.nc';
filename_dst_lo    = 'exp_II_topography_lo_5km.nc';

% Topography parameters
f_TAF = 0.1;

% Read data
x    = ncread( filename_src,'x'); nx = length( x);
y    = ncread( filename_src,'y'); ny = length( y);
time = ncread( filename_src,'time');
ti   = length( time);

topo.Hi  = ncread( filename_src,'Hi',[1,1,ti],[Inf,Inf,1]);
topo.Hb  = ncread( filename_src,'Hb',[1,1,ti],[Inf,Inf,1]);
topo.Hs  = ncread( filename_src,'Hs',[1,1,ti],[Inf,Inf,1]);

ice_density                      =  910.0;
seawater_density                 = 1028.0;
topo.TAF = topo.Hi - max(0, (-topo.Hb) * (seawater_density / ice_density));
m = topo.TAF > 0;

% Create perturbed versions
topo_hi = topo;
topo_hi.Hi( m) = topo_hi.Hi( m) + topo.TAF( m) * f_TAF;
topo_hi.Hb( m) = topo_hi.Hb( m) - topo.TAF( m) * f_TAF;

topo_lo = topo;
topo_lo.Hi( m) = topo_lo.Hi( m) - topo.TAF( m) * f_TAF;
topo_lo.Hb( m) = topo_lo.Hb( m) + topo.TAF( m) * f_TAF;

% Delete files if necessary
if exist( filename_dst_hi,'file')
  delete( filename_dst_hi)
end
if exist( filename_dst_lo,'file')
  delete( filename_dst_lo)
end

copyfile( filename_src_topo,filename_dst_hi);
copyfile( filename_src_topo,filename_dst_lo);

% Write perturbed data to new files
ncwrite( filename_dst_hi,'Hi',topo_hi.Hi);
ncwrite( filename_dst_hi,'Hb',topo_hi.Hb);
ncwrite( filename_dst_hi,'Hs',topo_hi.Hs);

ncwrite( filename_dst_lo,'Hi',topo_lo.Hi);
ncwrite( filename_dst_lo,'Hb',topo_lo.Hb);
ncwrite( filename_dst_lo,'Hs',topo_lo.Hs);
