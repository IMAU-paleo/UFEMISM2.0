clc
clear all
close all

% filename = '../../test_mesh.txt';
% mesh = read_mesh_from_text_file( filename);
% figure; patch('vertices',mesh.V,'faces',mesh.Tri,'facecolor','none');

filename = '../../test_mesh_netcdf.nc';
mesh = read_mesh_from_file( filename);
figure; patch('vertices',mesh.V,'faces',mesh.Tri,'facecolor','interp',...
  'facevertexcdata',mesh.lon,'edgecolor','none');