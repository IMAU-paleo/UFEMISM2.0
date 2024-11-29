function vi = find_containing_vertex( mesh, p, vi)
% Find the vertex whose Voronoi cell contains the point p, using a "linear search"
% Start at initial guess vi. Check all neighbours of vi, find the one
% closest to p, select that one as the new vi. Repeat until all neighbours
% of vi are further away from p than vi itself.

ncycle = 0;
vi_prev = vi;
while (ncycle < mesh.nV)

  d = norm( mesh.V( vi,:) - p);

  dcmin = d + 10;
  vcmin = 0;
  for ci = 1: mesh.nC( vi)
    vc = mesh.C( vi,ci);
    if (vc == vi_prev); continue; end % This is the neighbour we just came from
    dc = norm( mesh.V( vc,:) - p);
    if (dc < dcmin)
      dcmin = dc;
      vcmin = vc;
    end
  end

  if (dcmin < d)
    vi_prev = vi;
    vi = vcmin;
  else
    return
  end

end

% If we reach this point, we didnt find the containing vertex - should not be possible, so throw an error
error('find_containing_vertex - couldnt find closest vertex!')

end