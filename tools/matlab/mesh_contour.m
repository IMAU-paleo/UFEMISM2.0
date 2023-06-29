function C = mesh_contour( mesh,d,level)
% Given a data field d on the mesh, return a list of x,y points describing
% the contour lines at level.

C    = [];
Csub = zeros( mesh.nTri*2,2);
tri_checked = false( mesh.nTri,1);

for ti = 1:mesh.nTri
  
  if tri_checked( ti); continue; end
  
  via = mesh.Tri( ti,1);
  vib = mesh.Tri( ti,2);
  vic = mesh.Tri( ti,3);
  
  dmin = min( min( d( via),d( vib)),d( vic));
  dmax = max( max( d( via),d( vib)),d( vic));
  
  % If the level is not within the data range of this triangle, continue;
  if dmin > level || dmax < level
    continue
  end
  
  % Follow the contour line until we're back in the current triangle
  Csub(:,:) = 0;
  nC        = 0;
  ti_cur    = ti;
  ti_prev   = 0;  
  Finished  = false;
  while (~Finished)
    % Find the next triangle
    for n1 = 1: 3
      n2 = n1+1; if (n2 == 4); n2 = 1; end
      n3 = n2+1; if (n3 == 4); n3 = 1; end
      if (d( mesh.Tri( ti_cur,n1)) >= level && d( mesh.Tri( ti_cur,n2)) <= level) || ...
         (d( mesh.Tri( ti_cur,n1)) <= level && d( mesh.Tri( ti_cur,n2)) >= level)
        ti_next = mesh.TriC( ti_cur,n3);
        if (ti_next ~= ti_prev && ti_next ~= 0)
          % Interpolate position between the two vertices
          w = (level - d( mesh.Tri( ti_cur,n1))) / (d( mesh.Tri( ti_cur,n2)) - d( mesh.Tri( ti_cur,n1)));
          p = mesh.V( mesh.Tri( ti_cur,n1),:) * (1 - w) + mesh.V( mesh.Tri( ti_cur,n2),:) * (w);
          nC = nC+1;
          Csub( nC,:) = p;
          break
        end
      end
    end
    % Cycle
    ti_prev = ti_cur;
    ti_cur  = ti_next;
    tri_checked( ti_cur) = true;
    % Check if finished
    if (ti_cur == ti); Finished = true; end
  end
  % Repeat last point to make it a closed loop
  nC = nC+1;
  Csub( nC,:) = Csub(1,:);
  
  % Add to total line
  if ~isempty(C); C = [C; [NaN,NaN]]; end
  C = [C; Csub(1:nC,:)];
  
end

end