clc
clear all
close all

filename = 'Ant_Bedmachine_mesh1.nc';

mesh = read_mesh_from_file( filename);
mesh.tol_dist = max( mesh.xmax-mesh.xmin, mesh.ymax-mesh.ymin) * 1e-9;

% Create test mask
mask = false( mesh.nV,1);
for vi = 1: mesh.nV
  if (norm( mesh.V( vi,:) - [mesh.xmin,mesh.ymin]) < mesh.xmax/2)
    mask( vi) = true;
  end
end

% Plot
plot_mesh_data_a( mesh, double( mask));

vi0 = 1;
[V_poly, n_poly] = calc_mesh_mask_as_polygon( mesh, mask, vi0);

function [V_poly, n_poly] = calc_mesh_mask_as_polygon( mesh, mask, vi0)
  % Calculate a polygon enveloping the set of TRUE-valued mask vertices around vi0,
  % and remove that set of vertices from the mask

%     IMPLICIT NONE
% 
%     % In/output variables
%     TYPE(type_mesh),                         INTENT(IN)    :: mesh
%     LOGICAL,  DIMENSION(:    ),              INTENT(INOUT) :: mask
%     INTEGER,                                 INTENT(IN)    :: vi0
%     REAL(dp), DIMENSION(:,:  ),              INTENT(OUT)   :: poly
%     INTEGER,                                 INTENT(OUT)   :: n_poly
% 
%     % Local variables:
%     CHARACTER(LEN=256), PARAMETER                          :: routine_name = 'calc_mesh_mask_as_polygon'
%     INTEGER,  DIMENSION(:    ), ALLOCATABLE                :: map, stack
%     INTEGER                                                :: stackN
%     INTEGER                                                :: vi, ci, vj, ei0, ei, ei_next, it, cii, ci2, vk
% 
%     % Add routine to path
%     CALL init_routine( routine_name)
    
    V_poly = zeros( mesh.nE,2);

    % Safety
    if (size( mask,1) ~= mesh.nV); error('incorrect data dimensions!'); end
    if (~mask( vi0)); error('seed at vi0 is not TRUE!'); end

    % Use a flood-fill algorithm to find the map of same-valued grid cells around mask cell i0,j0
    map   = zeros( mesh.nV,1);
    stack = zeros( mesh.nV,1);

    map( vi0) = 1;
    stackN = 1;
    stack( 1) = vi0;

    while (stackN > 0)

      % Take the last element from the stack
      vi = stack( stackN);
      stackN = stackN - 1;

      % Mark it as mapped
      map( vi) = 2;

      % Remove it from the input mask
      mask( vi) = false;

      % Add all non-mapped, non-stacked, TRUE-valued neighbours to the stack
      for ci = 1: mesh.nC( vi)
        vj = mesh.C( vi,ci);
        % If this neighbour lies outside the TRUE region, store the connection
        % as a starting point for the outline tracer
        if (map( vj) == 0 && ~mask( vj) && ~(mesh.VBI( vi) > 0 && mesh.VBI( vj) > 0)); ei0 = mesh.VE( vi,ci); end
        % If this neighbour lies inside the TRUE region and isn't
        % marked yet, add it to the stack and mark it
        if (map( vj) == 0 && mask( vj))
          % Add this neighbour to the stack
          stackN = stackN + 1;
          stack( stackN) = vj;
          % Mark this neighbour on the map as stacked
          map( vj) = 1;
        end
      end

    end % DO WHILE (stackN > 0)
    % Safety
    if (ei0 == 0); error('couldnt find starting edge!'); end

    % Starting at the edge we found earlier, trace the outline of the TRUE region
    ei     = ei0;
    n_poly = 0;
    it     = 0;
    
    while true

      % Safety
      it = it + 1;
      if (it > mesh.nE); error('outline tracer got stuck!'); end

      % Find the vertices vi and vj spanning the current edge ei,
      % sorted such that vi lies inside the TRUE region and vj lies outside of it.
      if     (map( mesh.EV( ei,1)) == 2 && map( mesh.EV( ei,2)) == 0)
        vi  = mesh.EV(   ei,1);
        vj  = mesh.EV(   ei,2);
        til = mesh.ETri( ei,1);
        tir = mesh.ETri( ei,2);
      elseif (map( mesh.EV( ei,1)) == 0 && map( mesh.EV( ei,2)) == 2)
        vi  = mesh.EV(   ei,2);
        vj  = mesh.EV(   ei,1);
        til = mesh.ETri( ei,2);
        tir = mesh.ETri( ei,1);
      else
        % Apparently this edge doesn't cross the border of the TRUE region
        error('found non-border edge!')
      end
      
      % Add this edge to the polygon
      n_poly = n_poly + 1;
      if (tir == 0)
        % vi-vj is a border edge?
        if (mesh.VBI( vi) > 0 && mesh.VBI( vj) > 0)
          V_poly( n_poly,:) = (mesh.V( vi,:) + mesh.V( vj,:)) / 2.0;
        else
          error('expected vi-vj to be a border edge!')
        end
      else
        V_poly( n_poly,:) = mesh.Tricc( tir,:);
      end
    
      % DENK DROM
      if it > 1
        line('xdata',[V_poly( n_poly-1,1),V_poly( n_poly,1)],'ydata',[V_poly( n_poly-1,2),V_poly( n_poly,2)],...
          'color','r','linewidth',2)
        drawnow('update')
%         pause
      end

      % Find the index ci such that VE( vi,ci) = ei
      ci = 0;
      for cii = 1: mesh.nC( vi)
        if (mesh.VE( vi,cii) == ei)
          ci = cii;
          break
        end
      end
      % Safety
      if (ci == 0); error('couldnt find connection index ci such that VE( vi,ci) = ei!'); end
      
      if (mesh.VBI( vi) > 0 && ci == mesh.nC( vi))
        % Special case for when the tracer reaches the domain border
        
        % Move along the border until we find the other end of the TRUE region
        
        % Add the last Voronoi cell boundary section
        n_poly = n_poly + 1;
        V_poly( n_poly,:) = (mesh.V( vi,:) + mesh.V( vj,:)) / 2.0;
        % DENK DROM
        line('xdata',[V_poly( n_poly-1,1),V_poly( n_poly,1)],'ydata',[V_poly( n_poly-1,2),V_poly( n_poly,2)],...
          'color','r','linewidth',2)
        
        % Add the section from the Voronoi cell boundary to vi
        n_poly = n_poly + 1;
        V_poly( n_poly,:) = mesh.V( vi,:);
        % DENK DROM
        line('xdata',[V_poly( n_poly-1,1),V_poly( n_poly,1)],'ydata',[V_poly( n_poly-1,2),V_poly( n_poly,2)],...
          'color','r','linewidth',2)
        
        % Add all border sections until we find the other end of the TRUE region
        it2 = 0;
        while (map( mesh.C( vi,1)) == 2)
          % Safety
          it2 = it2 + 1;
          if (it2 > mesh.nV); error('outline tracer for mesh border got stuck!'); end
          % Move to next border vertex
          vi = mesh.C( vi,1);
          % Add section to polygon
          n_poly = n_poly + 1;
          V_poly( n_poly,:) = mesh.V( vi,:);
          % DENK DROM
          line('xdata',[V_poly( n_poly-1,1),V_poly( n_poly,1)],'ydata',[V_poly( n_poly-1,2),V_poly( n_poly,2)],...
            'color','r','linewidth',2)
        end
        
        % The next edge
        ei = mesh.VE( vi,1);
        
      else % if (mesh.VBI( vi) > 0 && ci == mesh.nC( vi))
        % Regular case for the domain interior

        % The index ci2 of the next edge originating from vi, counterclockwise from ci
        ci2 = ci + 1;
        if (ci2 > mesh.nC( vi)); ci2 = 1; end
        
        % The vertex that this edge points to
        vk = mesh.C( vi,ci2);
        
        if (map( vk) == 2)
          % If vk also lies inside the TRUE region, move to its Voronoi cell boundary
          
          % Find the connection index ck such that C( vk,ci) = vi
          ck = 0;
          for cii = 1: mesh.nC( vk)
            if (mesh.C( vk,cii) == vi)
              ck = cii;
              break
            end
          end
          % Safety
          if (ck == 0); error('couldnt find connection index ck such that C( vk,ci) = vi'); end

          % The index ck2 of the next edge originating from vk, counterclockwise from ck
          ck2 = ck + 1;
          if (ck2 > mesh.nC( vk)); ck2 = 1; end
          
          % The next edge
          ei = mesh.VE( vk,ck2);
          
        else % if (map( vk) == 2)
          % If vk lies outside the TRUE region, move to the next edge along vi
          
          % The next edge
          ei = mesh.VE( vi,ci2);
          
        end % if (map( vk) == 2)
        
      
      end % if (mesh.VBI( vi) > 0 && ci == mesh.nC( vi))

      % If we've reached the starting point again, stop
      if (it > 1)
        if (norm( V_poly( n_poly,:) - V_poly( 1,:)) < mesh.tol_dist)
          break
        end
      end

    end % DO WHILE (.TRUE.)
% 
%     % Finalise routine path
%     CALL finalise_routine( routine_name)

end % calc_mesh_mask_as_polygon