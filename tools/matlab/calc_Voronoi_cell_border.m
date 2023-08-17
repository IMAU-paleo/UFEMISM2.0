function [Vor, Vor_vi, Vor_ti, nVor] = calc_Voronoi_cell_border( mesh, vi, dx)
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
% 
% % Local variables:
% INTEGER                                            :: ci,vj,ti,iti,tj,n,n2

% Initialise
Vor    = zeros( mesh.nC_mem,2);
Vor_vi = zeros( mesh.nC_mem,1);
Vor_ti = zeros( mesh.nC_mem,1);
nVor   = 0;

% For each neighbouring vertex vj, list the circumcentre of the triangle
% lying right of the line vi-vj. Except, we skip the first neighbouring
% vertex, as DO that one the area to the right of the line vi-vj is
% outside of the mesh domain

for ci = 2: mesh.nC( vi)

  vj  = mesh.C( vi,ci);

  % Find the triangle ti lying right of the line vi-vj
  ti = 0;
  for iti = 1: mesh.niTri( vi)
    tj = mesh.iTri( vi,iti);
    for n = 1: 3
      n2 = n + 1;
      if (n2 == 4); n2 = 1; end
      if (mesh.Tri( tj,n) == vj && mesh.Tri( tj,n2) == vi)
        ti = tj;
        break
      end
    end
    if (ti > 0); break; end
  end

  % Safety
  if (ti == 0); error('couldnt find triangle right of the line vi-vj!'); end

  % Safety
  if (mesh.Tricc( ti,1) < mesh.xmin || mesh.Tricc( ti,1) > mesh.xmax || ...
      mesh.Tricc( ti,2) < mesh.ymin || mesh.Tricc( ti,2) > mesh.ymax)
    error('calc_Voronoi_cell_free - found triangle circumcentre outside the mesh domain!')
  end

  % List the circumcentre of this triangle
  nVor = nVor + 1;
  Vor(     nVor,:) = mesh.Tricc( ti,:);
  Vor_vi(  nVor  ) = vj;
  Vor_ti(  nVor  ) = ti;

end % DO ci = 2, mesh.nC( vi)

% == Add the projection of the first Voronoi vertex on the domain boundary as an additional point

nVor = nVor + 1;
Vor(     2:nVor,:) = Vor(     1:nVor-1,:);
Vor_vi(  2:nVor  ) = Vor_vi(  1:nVor-1  );
Vor_ti(  2:nVor  ) = Vor_ti(  1:nVor-1  );

Vor_vi( 1) = mesh.C(    vi,1);
Vor_ti( 1) = mesh.iTri( vi,1);

if     (mesh.VBI( vi) == 1 || mesh.VBI( vi) == 2)
  % vi is on the northern border, or on the northeast corner; its first neighbour lies on the northern border
  Vor(    1,:) = [Vor( 2,1), mesh.ymax + dx];
elseif (mesh.VBI( vi) == 3 || mesh.VBI( vi) == 4)
  % vi is on the eastern border, or on the southeast corner; its first neighbour lies on the eastern border
  Vor(    1,:) = [mesh.xmax + dx, Vor( 2,2)];
elseif (mesh.VBI( vi) == 5 || mesh.VBI( vi) == 6)
  % vi is on the southern border, or on the southwest corner; its first neighbour lies on the southern border
  Vor(    1,:) = [Vor( 2,1), mesh.ymin - dx];
elseif (mesh.VBI( vi) == 7 || mesh.VBI( vi) == 8)
  % vi is on the western border, or on the northwest corner; its first neighbour lies on the western border
  Vor(    1,:) = [mesh.xmin - dx, Vor( 2,2)];
end

% == Add the projection of the last Voronoi vertex on the domain boundary as an additional point

nVor = nVor + 1;

Vor_vi( nVor) = mesh.C(    vi, mesh.nC( vi));
Vor_ti( nVor) = mesh.iTri( vi, mesh.nC( vi));

if     (mesh.VBI( vi) == 2 || mesh.VBI( vi) == 3)
  % vi is on the eastern border, or on the northeast corner; its last neighbour lies on the eastern border
  Vor(    nVor,:) = [mesh.xmax + dx, Vor( nVor-1,2)];
elseif (mesh.VBI( vi) == 4 || mesh.VBI( vi) == 5)
  % vi is on the southern border, or on the southeast corner; its last neighbour lies on the southern border
  Vor(    nVor,:) = [Vor( nVor-1,1), mesh.ymin - dx];
elseif (mesh.VBI( vi) == 6 || mesh.VBI( vi) == 7)
  % vi is on the western border, or on the southwest corner; its last neighbour lies on the western border
  Vor(    nVor,:) = [mesh.xmin - dx, Vor( nVor-1,2)];
elseif (mesh.VBI( vi) == 8 || mesh.VBI( vi) == 1)
  % vi is on the northern border, or on the northwest corner; its last neighbour lies on the northern border
  Vor(    nVor,:) = [Vor( nVor-1,1), mesh.ymax + dx];
end

% In the case of the four corners, add the corners themselves as well
if     (mesh.VBI( vi) == 2)
  % Northeast corner
  nVor = nVor + 1;
  Vor( nVor,:) = [mesh.xmax + dx, mesh.ymax + dx];
elseif (mesh.VBI( vi) == 4)
  % Southeast corner
  nVor = nVor + 1;
  Vor( nVor,:) = [mesh.xmax + dx, mesh.ymin - dx];
elseif (mesh.VBI( vi) == 6)
  % Southwest corner
  nVor = nVor + 1;
  Vor( nVor,:) = [mesh.xmin - dx, mesh.ymin - dx];
elseif (mesh.VBI( vi) == 8)
  % Northwest corner
  nVor = nVor + 1;
  Vor( nVor,:) = [mesh.xmin - dx, mesh.ymax + dx];
end

end