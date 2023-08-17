function H = plot_mesh_data_Voronoi( mesh, d, varargin)

if (size(d,1) ~= mesh.nV); error('data field is the wrong size!'); end

V_Vor  = zeros( mesh.nV * mean( mesh.nC), 2);
nV_Vor = 0;
F_Vor  = NaN( mesh.nV, mesh.nC_mem);

for vi = 1: mesh.nV
  
  [Vor, ~, ~, nVor] = calc_Voronoi_cell( mesh, vi, 0.0);
  
  V_Vor( nV_Vor + 1 : nV_Vor + nVor,:) = Vor( 1:nVor,:);
  F_Vor( vi,1:nVor) = nV_Vor + 1 : nV_Vor + nVor;
  nV_Vor = nV_Vor + nVor;
  
  if nV_Vor > size( V_Vor,1)-mesh.nC_mem
    V_Vor = [V_Vor; zeros( 100 * mesh.nC_mem,2)];
  end
  
end

V_Vor = V_Vor( 1:nV_Vor,:);

% Axes and figure size
xw = mesh.xmax - mesh.xmin;
yw = mesh.ymax - mesh.ymin;
if xw >= yw
  wa = 800;
  ha = wa * yw / xw;
else
  ha = 800;
  wa = ha * xw / yw;
end
wf = 25 + wa + 100;
hf = 25 + ha + 50;

H.Fig = figure('position',[200,200,wf,hf],'color','w');
H.Ax  = axes('parent',H.Fig,'units','pixels','position',[25,25,wa,ha],'fontsize',24,...
  'xtick',[],'ytick',[],'xlim',[mesh.xmin,mesh.xmax],'ylim',[mesh.ymin,mesh.ymax]);

if ~isempty( varargin)
  H.Patch = patch('vertices',V_Vor,'faces',F_Vor,'facecolor','flat','facevertexcdata',d,'edgecolor',varargin{1});
else
  H.Patch = patch('vertices',V_Vor,'faces',F_Vor,'facecolor','flat','facevertexcdata',d,'edgecolor','none');
end

pos = get( H.Ax,'position');
H.Cbar = colorbar( H.Ax,'location','eastoutside');
set( H.Ax,'position',pos);

set( H.Ax,'units','normalized');

end