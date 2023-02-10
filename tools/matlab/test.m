clc
clear all
close all

filename = '../../test_mesh.txt';
mesh = read_mesh_from_text_file( filename);
figure; patch('vertices',mesh.V,'faces',mesh.Tri,'facecolor','none');