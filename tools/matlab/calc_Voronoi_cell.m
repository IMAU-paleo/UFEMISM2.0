function [Vor, Vor_vi, Vor_ti, nVor] = calc_Voronoi_cell( mesh, vi, dx)
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
% REAL(dp),                            INTENT(IN)    :: dx
% REAL(dp), DIMENSION(mesh.nC_mem,2),  INTENT(OUT)   :: Vor
% INTEGER,  DIMENSION(mesh.nC_mem  ),  INTENT(OUT)   :: Vor_vi
% INTEGER,  DIMENSION(mesh.nC_mem  ),  INTENT(OUT)   :: Vor_ti
% INTEGER,                             INTENT(OUT)   :: nVor

% Local variables:

% Vertices lying on the domain border or corners are best treated separately
if (mesh.VBI( vi) == 0)
  % Free vertex
  [Vor, Vor_vi, Vor_ti, nVor] = calc_Voronoi_cell_free(   mesh, vi);
else
  % Border vertex
  [Vor, Vor_vi, Vor_ti, nVor] = calc_Voronoi_cell_border( mesh, vi, dx);
end

end