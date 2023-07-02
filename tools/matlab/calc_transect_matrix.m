function A = calc_transect_matrix( mesh, xt, yt)

n = length(xt);

A = sparse( n, mesh.nV);

mesh.tol_dist  = 1e-7;
mesh.TriMap    = zeros( mesh.nTri,1);
mesh.TriStack1 = zeros( mesh.nTri,1);
mesh.TriStack2 = zeros( mesh.nTri,1);

ti = 1;
for i = 1:n
  
  % The point of the transect
  p = [xt( i), yt( i)];
  
  % The triangle containing this point
  ti = find_containing_triangle( mesh, p, ti);
  
  % The three vertices spanning the triangle
  via = mesh.Tri( ti,1);
  vib = mesh.Tri( ti,2);
  vic = mesh.Tri( ti,3);
  pa  = mesh.V( via,:);
  pb  = mesh.V( vib,:);
  pc  = mesh.V( vic,:);
  
  % The three interpolation weights
  Atri_abp = calc_triangle_area( pa, pb, p);
  Atri_bcp = calc_triangle_area( pb, pc, p);
  Atri_cap = calc_triangle_area( pc, pa, p);
  Atri_abc = Atri_abp + Atri_bcp + Atri_cap;

  wa = Atri_bcp / Atri_abc;
  wb = Atri_cap / Atri_abc;
  wc = Atri_abp / Atri_abc;
  
  % Fill them into the sparse matrix
  A( i, via) = wa;
  A( i, vib) = wb;
  A( i, vic) = wc;
  
end
  
end