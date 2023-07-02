clc
clear all
close all

filename = '../../results_benchmark_MISMIPplus_spinup_with_spatial_roughness/main_output_ANT_00001.nc';

mesh = read_mesh_from_file( filename);

trans.n = 100;
trans.x = linspace( mesh.xmin, mesh.xmax, trans.n)';
trans.y = zeros( trans.n,1);
trans.M = calc_transect_matrix( mesh, trans.x, trans.y);

time = ncread( filename,'time');
ti = length( time);

Hi  = ncread( filename,'Hi',[1,ti],[Inf,1]);
Hb  = ncread( filename,'Hb',[1,ti],[Inf,1]);
Hs  = ncread( filename,'Hs',[1,ti],[Inf,1]);
Hib = Hs - Hi;

trans.Hi  = trans.M * Hi;
trans.Hb  = trans.M * Hb;
trans.Hs  = trans.M * Hs;
trans.Hib = trans.M * Hib;

wa = 600;
ha = 400;

margins_hor = [110,25];
margins_ver = [25,80];

H = setup_multipanel_figure( wa, ha, margins_hor, margins_ver);

xlabel( H.Ax{ 1,1},'x (km)')
ylabel( H.Ax{ 1,1},'z (m)')

line( 'parent',H.Ax{ 1,1},'xdata',trans.x / 1e3,'ydata',trans.Hb ,'linewidth',3,'color','k');
line( 'parent',H.Ax{ 1,1},'xdata',trans.x / 1e3,'ydata',trans.Hs ,'linewidth',3,'color','k');
line( 'parent',H.Ax{ 1,1},'xdata',trans.x / 1e3,'ydata',trans.Hib,'linewidth',3,'color','k');