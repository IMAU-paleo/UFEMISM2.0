clc
clear all
close all

filename_src       = '../../exp_I_target_20km/restart_ANT.nc';
filename_src_topo  = 'exp_I_topography_20km.nc';
filename_dst_hi    = 'exp_I_topography_hi_20km.nc';
filename_dst_lo    = 'exp_I_topography_lo_20km.nc';

% Topography parameters
f_Hi_hi = 1.1;
f_Hi_lo = 0.9;

% Read data
x    = ncread( filename_src,'x');
y    = ncread( filename_src,'y');
time = ncread( filename_src,'time');
ti   = length( time);

topo.Hi = ncread( filename_src,'Hi',[1,1,ti],[Inf,Inf,1]);
topo.Hb = ncread( filename_src,'Hb',[1,1,ti],[Inf,Inf,1]);
topo.Hs = ncread( filename_src,'Hs',[1,1,ti],[Inf,Inf,1]);

% Create perturbed versions
topo_hi = topo;
topo_hi.Hi = topo.Hi * f_Hi_hi;
topo_hi.Hb = topo_hi.Hs - topo_hi.Hi;

topo_lo = topo;
topo_lo.Hi = topo.Hi * f_Hi_lo;
topo_lo.Hb = topo_lo.Hs - topo_lo.Hi;

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