function C = calc_mesh_contour( mesh, d, level)
% Return a contour in the same format as Matlab's own 'contour' function
% See: https://nl.mathworks.com/help/matlab/ref/contour.html

C  = zeros( mesh.nE,2);
nC = 0;

d_scaled = scale_d( d, level);

% Mark all edges that cross the contour, and
% count number of contoour-crossing edges per triangle
E_cross_C = false( mesh.nE,1);
nT_cross_C = zeros( mesh.nTri,1);
for ei = 1: mesh.nE
  vi = mesh.EV( ei,1);
  vj = mesh.EV( ei,2);
  if (d_scaled( vi) * d_scaled( vj) < 0)
    E_cross_C( ei) = true;
    ti = mesh.ETri( ei,1);
    tj = mesh.ETri( ei,2);
    if (ti > 0); nT_cross_C( ti) = nT_cross_C( ti) + 1; end
    if (tj > 0); nT_cross_C( tj) = nT_cross_C( tj) + 1; end
  end
end

% Mark edges where a contour ends
E_end = false( mesh.nE,1);
for ei = 1: mesh.nE
  if (E_cross_C( ei))
    if (mesh.EBI( ei) > 0)
      E_end( ei) = true;
    else
      is_next_to_single_edge_tri = false;
      ti = mesh.ETri( ei,1);
      if (ti > 0)
        if (nT_cross_C( ti) == 1)
          is_next_to_single_edge_tri = true;
        end
      end
      tj = mesh.ETri( ei,2);
      if (tj > 0)
        if (nT_cross_C( tj) == 1)
          is_next_to_single_edge_tri = true;
        end
      end
      if (is_next_to_single_edge_tri)
        E_end( ei) = true;
      end
    end
  end
end

% Trace linear contours
for ei = 1: mesh.nE
  if (E_end( ei))
    [C_sub, n_sub, E_end, E_cross_C] = trace_linear_contour( ...
      mesh, d_scaled, E_end, E_cross_C, ei);
    % Add to total contour list (native Matlab format)
    C( nC+1,:) = [n_sub, NaN];
    C( nC+2: nC+n_sub+1,:) = C_sub( 1:n_sub,:);
    nC = nC + n_sub + 1;
  end
end

% Trace circular contours
for ei = 1: mesh.nE
  if (E_cross_C( ei))
    [C_sub, n_sub, E_cross_C] = trace_circular_contour( ...
      mesh, d_scaled, E_cross_C, ei);
    % Add to total contour list (native Matlab format)
    C( nC+1,:) = [n_sub, NaN];
    C( nC+2: nC+n_sub+1,:) = C_sub( 1:n_sub,:);
    nC = nC + n_sub + 1;
  end
end

C = C( 1:nC,:);
  
  function d_scaled = scale_d( d, level)
    % Scale d so that it runs from -1 to 1, with 0 corresponding to level
    d_scaled = d - level;
    d_scaled( d_scaled > 0) =  d_scaled( d_scaled > 0) / max( d_scaled( d_scaled > 0));
    d_scaled( d_scaled < 0) = -d_scaled( d_scaled < 0) / min( d_scaled( d_scaled < 0));
  end
  
  function [C_sub, n_sub, E_end, E_cross_C] = trace_linear_contour( ...
    mesh, d_scaled, E_end, E_cross_C, ei)
  
    E_C = zeros( mesh.nE,1);
    n_C = 0;
    
    ei_prev = 0;
    
    nit = 0;
    while (nit < mesh.nE)
      nit = nit+1;
    
      % Add current edge to the contour
      n_C = n_C + 1;
      E_C( n_C) = ei;
    
      % Find next edge (if any)
      ei_next = 0;
    
      % Try both adjacent triangles
      ti = mesh.ETri( ei,1);
      if (ti > 0)
        for n = 1: 3
          ej = mesh.TriE( ti,n);
          if (ej ~= ei && ej ~= ei_prev && E_cross_C( ej))
            ei_next = ej;
          end
        end
      end
      tj = mesh.ETri( ei,2);
      if (tj > 0)
        for n = 1: 3
          ej = mesh.TriE( tj,n);
          if (ej ~= ei && ej ~= ei_prev && E_cross_C( ej))
            ei_next = ej;
          end
        end
      end
    
      % If no next edge could be found, we probably 
      % reached the end of the linear contour
      reached_end = false;
      if (ei_next == 0)
        if (~E_end( ei))
          error('whaa!')
        elseif (ei_prev == 0)
          error('whaa!')
        else
          reached_end = true;
        end
      end
    
      if (reached_end)
        break
      else
        % Move a step along the contour
        % line('xdata',[mesh.E(ei,1),mesh.E(ei_next,1)],...
        %   'ydata',[mesh.E(ei,2),mesh.E(ei_next,2)],'color','r','linewidth',2)
        % drawnow('update')
        ei_prev = ei;
        ei = ei_next;
      end
    
    end
    
    % Unmark the edges of this contour on the maps
    for i = 1: n_C
      ei = E_C( i);
      E_cross_C( ei) = false;
      E_end( ei) = false;
    end
    
    % Add precise contour coordinates to the list
    n_sub = n_C;
    C_sub = zeros( n_sub,2);
    for i = 1: n_sub
      ei = E_C( i);
      vi = mesh.EV( ei,1);
      vj = mesh.EV( ei,2);
      C_sub( i,:) = linint_points( mesh.V( vi,:), mesh.V( vj,:), ...
        d_scaled( vi), d_scaled( vj), 0);
    end
  
  end
    
  function [C_sub, n_sub, E_cross_C] = trace_circular_contour( ...
    mesh, d_scaled, E_cross_C, ei_start)
  
    E_C = zeros( mesh.nE,1);
    n_C = 0;
    
    ei = ei_start;
    ei_prev = 0;
    
    nit = 0;
    while (nit < mesh.nE)
      nit = nit+1;
    
      % Add current edge to the contour
      n_C = n_C + 1;
      E_C( n_C) = ei;
    
      % Find next edge (if any)
      ei_next = 0;
    
      % Try both adjacent triangles
      ti = mesh.ETri( ei,1);
      if (ti > 0)
        for n = 1: 3
          ej = mesh.TriE( ti,n);
          if (ej ~= ei && ej ~= ei_prev && E_cross_C( ej))
            ei_next = ej;
          end
        end
      end
      tj = mesh.ETri( ei,2);
      if (tj > 0)
        for n = 1: 3
          ej = mesh.TriE( tj,n);
          if (ej ~= ei && ej ~= ei_prev && E_cross_C( ej))
            ei_next = ej;
          end
        end
      end
    
      % Safety
      if (ei_next == 0)
        error('whaa!')
      end
    
      if (ei_next == ei_start)
        n_C = n_C + 1;
        E_C( n_C) = ei_next;
        break
      else
        % Move a step along the contour
        ei_prev = ei;
        ei = ei_next;
      end
    
    end
    
    % Umark the edges of this contour on the maps
    for i = 1: n_C
      ei = E_C( i);
      E_cross_C( ei) = false;
    end
    
    % Add precise contour coordinates to the list
    n_sub = n_C;
    C_sub = zeros( n_sub,2);
    for i = 1: n_sub
      ei = E_C( i);
      vi = mesh.EV( ei,1);
      vj = mesh.EV( ei,2);
      C_sub( i,:) = linint_points( mesh.V( vi,:), mesh.V( vj,:), ...
        d_scaled( vi), d_scaled( vj), 0);
    end
  
  end

  function r = linint_points( p, q, fp, fq, f0)
    % Given a function f( x) and points p, q such that f( p) = fp, f( q) = fq,
    % interpolate f linearly to find the point r such that f( r) = f0
  
    % Safety - if fp == fq, then f = fp = fq everywhere
    if (abs( 1 - fp/fq) < 1E-9)
      r = (p + q) / 2;
      return
    end
  
    x1 = 0;
    x2 = 1;
  
    lambda = (fq - fp) / (x2 - x1);
    w = x1 + (f0 - fp) / lambda;
  
    r = w * q + (1-w) * p;
  
  end

end