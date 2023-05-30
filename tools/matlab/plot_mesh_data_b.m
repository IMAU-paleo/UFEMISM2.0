function H = plot_mesh_data_b( mesh, d)

edgecolor = 'none';
% edgecolor = 'k';

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
wf = 25 + wa + 50;
hf = 25 + ha + 50;

H.Fig = figure('position',[200,200,wf,hf],'color','w');
H.Ax  = axes('parent',H.Fig,'units','pixels','position',[25,25,wa,ha],'fontsize',24,...
  'xtick',[],'ytick',[],'xlim',[mesh.xmin,mesh.xmax],'ylim',[mesh.ymin,mesh.ymax]);
H.Patch = patch('vertices',mesh.V( 1:mesh.nV,:),'faces',mesh.Tri( 1:mesh.nTri,:),...
  'facecolor','flat','facevertexcdata',d,'edgecolor',edgecolor);

end