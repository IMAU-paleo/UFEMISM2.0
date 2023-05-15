function H = plot_mesh_data( mesh, d)

if length(d) == mesh.nV
  H = plot_mesh_data_a( mesh,d);
elseif length(d) == mesh.nTri
  H = plot_mesh_data_b( mesh,d);
elseif length(d) == mesh.nE
  H = plot_mesh_data_c( mesh,d);
else
  error('unexpected vector length!')
end

end