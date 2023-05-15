clc
clear all
close all

filename = 'testfile.nc';
filename = ['../../' filename];

mesh = read_mesh_from_file( filename);
R = ncread( filename,'R');
plot_mesh( mesh);
H = plot_mesh_data_a( mesh, R);
set( H.Ax,'colorscale','log','clim',[1e3,500e3]);
colormap( H.Ax, flipud( parula( 32)));
colorbar( H.Ax)