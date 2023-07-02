function t = find_containing_triangle( mesh, p, t)
  % Start at initial guess t, search outward from there using a
  % flood-fill algorithm.

  % If p lies outside the mesh domain, throw an error
  if (p(1) < mesh.xmin || p(1) > mesh.xmax || p(2) < mesh.ymin || p(2) > mesh.ymax) 
    error('FindContainingTriangle - ERROR: point lies outside mesh domain%')
  end

  % See if the initial guess is correct.
  q = mesh.V(mesh.Tri(t,1),:);
  r = mesh.V(mesh.Tri(t,2),:);
  s = mesh.V(mesh.Tri(t,3),:);
  if (is_in_triangle(q, r, s, p)); return; end

  % If not, start with a linear search.
  ncycle = 0;
  t_prev = t;
  while (ncycle < mesh.nTri)

    gcti = (mesh.V(mesh.Tri(t,1),:) + mesh.V(mesh.Tri(t,2),:) + mesh.V(mesh.Tri(t,3),:)) / 3;
    d = norm(gcti - p);

    dcmin = d + 10;
    tcmin = 0;
    for n = 1: 3
      tc   = mesh.TriC(t,n);
      if (tc==0);      continue; end % This triangle neighbour doesn't exist
      if (tc==t_prev); continue; end % This triangle neighbour is the one we just came from
      gctc = (mesh.V(mesh.Tri(tc,1),:) + mesh.V(mesh.Tri(tc,2),:) + mesh.V(mesh.Tri(tc,3),:)) / 3;
      dc = norm( gctc - p);
      if (dc < dcmin) 
        dcmin = dc;
        tcmin = tc;
      end
    end

    if (dcmin < d) 
      t_prev = t;
      t = tcmin;
    else
      break
    end

  end % while (ncycle < mesh.nTri)

  % Check if the result from the linear search is correct.
  q = mesh.V(mesh.Tri(t,1),:);
  r = mesh.V(mesh.Tri(t,2),:);
  s = mesh.V(mesh.Tri(t,3),:);
  if (is_in_triangle(q, r, s, p)); return; end
  if (lies_on_line_segment( q, r, p, mesh.tol_dist)); return; end
  if (lies_on_line_segment( r, s, p, mesh.tol_dist)); return; end
  if (lies_on_line_segment( s, q, p, mesh.tol_dist)); return; end

  % It's not. Perform a flood-fill style outward search.

  % Initialise map and stack.
  mesh.TriMap(:)    = 0;
  mesh.TriStack1(:) = 0;
  mesh.TriStack2(:) = 0;
  mesh.TriStackN1   = 0;
  mesh.TriStackN2   = 0;
  mesh.TriMap(t)    = 1; % We checked that one.

  % Add t' neighbours to the stack.
  for n = 1: 3
    if (mesh.TriC(t,n) > 0) 
      mesh.TriStackN1 = mesh.TriStackN1+1;
      mesh.TriStack1( mesh.TriStackN1) = mesh.TriC(t,n);
    end
  end

  FoundIt = false;
  while (~FoundIt)
    % Check all triangles in the stack. If they're not it, add their
    % non-checked neighbours to the new stack.

    mesh.TriStack2  = 0;
    mesh.TriStackN2 = 0;

    for n = 1: mesh.TriStackN1
      ti = mesh.TriStack1(n);
      q = mesh.V(mesh.Tri(ti,1),:);
      r = mesh.V(mesh.Tri(ti,2),:);
      s = mesh.V(mesh.Tri(ti,3),:);
      if (is_in_triangle(q, r, s, p))

        % Found it!
        t = ti;

        return

      else % if (is_in_triangle(ti,p))
        % Did not find it. And add this triangle's non-checked neighbours to the new stack.

        for n2 = 1: 3
          tin = mesh.TriC(ti,n2);
          if (tin==0);              continue; end % This neighbour doesn't exist.
          if (mesh.TriMap(tin)==1); continue; end % This neighbour has already been checked or is already in the stack.
          mesh.TriStackN2 = mesh.TriStackN2 + 1;
          mesh.TriStack2( mesh.TriStackN2) = tin;
          mesh.TriMap(tin) = 1;
        end

      end % if (is_in_triangle(q, r, s, p, tol)) 
    end % for n = 1: mesh.triStackN1

  % Cycle stacks.
  mesh.TriStack1  = mesh.TriStack2;
  mesh.TriStackN1 = mesh.TriStackN2;

  % If no more non-checked neighbours could be found, terminate and throw an error.
  if (mesh.TriStackN2==0) 
    error('FindContainingTriangle - ERROR: couldnt find triangle containing this point!')
  end

  end % while (~FoundIt)


end