function H = plot_mesh_data_c( mesh, d)

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
pos = get( H.Ax,'position');
H.Cbar = colorbar( H.Ax,'location','eastoutside');
set( H.Ax,'position',pos);

ncols = 256;
cmap = parula(ncols);
clim = [min(d), max(d)];
if (clim(1)==clim(2))
  clim(1) = clim(1)-1;
  clim(2) = clim(2)+1;
end
set( H.Ax,'clim',clim)

for ci = 1: ncols
  linedata(ci).x = [];
  linedata(ci).y = [];
end
for aci = 1: mesh.nE
  ci = (d( aci) - clim(1)) / (clim(2) - clim(1));
  ci = max(1,min(ncols,1 + round(ci*(ncols-1))));
  linedata(ci).x( end+1) = mesh.E( aci,1);
  linedata(ci).y( end+1) = mesh.E( aci,2);
end
for ci = 1: ncols
  line('parent',H.Ax,'xdata',linedata( ci).x,'ydata',linedata( ci).y,'linestyle','none',...
    'marker','o','markerfacecolor',cmap( ci,:),'markeredgecolor',cmap( ci,:),'markersize',8);
end

% lastly NaN values
line('parent',H.Ax,'xdata',mesh.E( isnan(d),1),'ydata',mesh.E( isnan(d),2),'linestyle','none',...
    'marker','x','markerfacecolor','r','markeredgecolor','r','markersize',8);

set( H.Ax,'units','normalized');

end