function [Vor, Vor_vi, Vor_ti, nVor] = calc_Voronoi_cell_free( mesh, vi)
% % Find the points spanning the Voronoi cell of vertex vi of the mesh.
% % Sorted counted-clockwise, with no double entries.
% !
% % Point Vor( i) corresponds to the circumcentre of triangle Vor_ti( i).
% !
% % The line connecting Vor( i) and Vor( i+1) is shared with the Voronoi cell
% % of vertex Vor_vi( i).
% 
% IMPLICIT NONE
% 
% % In/output variables:
% TYPE(type_mesh),                     INTENT(IN)    :: mesh
% INTEGER,                             INTENT(IN)    :: vi
% REAL(dp), DIMENSION(mesh.nC_mem,2),  INTENT(OUT)   :: Vor
% INTEGER,  DIMENSION(mesh.nC_mem  ),  INTENT(OUT)   :: Vor_vi
% INTEGER,  DIMENSION(mesh.nC_mem  ),  INTENT(OUT)   :: Vor_ti
% INTEGER,                             INTENT(OUT)   :: nVor
% 
% % Local variables:
% INTEGER                                            :: iti,ti,vj,n1,n2
% 
% % Initialise
Vor    = zeros( mesh.nC_mem,2);
Vor_vi = zeros( mesh.nC_mem,1);
Vor_ti = zeros( mesh.nC_mem,1);
nVor   = 0;

% List the circumcentres of all triangles surrounding vi
% as points spanning the Voronoi cell

for iti = 1: mesh.niTri( vi)

  ti  = mesh.iTri( vi,iti);

  % Find vertex vj such that triangle ti lies to the right of the line vi-vj
  vj = 0;
  for n1 = 1: 3
    n2 = n1 + 1;
    if (n2 == 4); n2 = 1; end
    if (mesh.Tri( ti,n2) == vi)
      vj = mesh.Tri( ti,n1);
      break
    end
  end
  % Safety
  if (vj == 0); error('calc_Voronoi_cell_free - couldnt find vertex vj in triangle ti!'); end

  % Safety
  if (mesh.Tricc( ti,1) < mesh.xmin || mesh.Tricc( ti,1) > mesh.xmax || ...
      mesh.Tricc( ti,2) < mesh.ymin || mesh.Tricc( ti,2) > mesh.ymax)
    error('calc_Voronoi_cell_free - found triangle circumcentre outside the mesh domain!')
  end

  % List the new Voronoi cell-spanning point
  nVor = nVor + 1;
  Vor(     nVor,:) = mesh.Tricc( ti,:);
  Vor_vi(  nVor  ) = vj;
  Vor_ti(  nVor  ) = ti;

end

end