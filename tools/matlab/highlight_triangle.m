function H = highlight_triangle( mesh, ti)

via = mesh.Tri( ti,1);
vib = mesh.Tri( ti,2);
vic = mesh.Tri( ti,3);

va = mesh.V( via,:);
vb = mesh.V( vib,:);
vc = mesh.V( vic,:);

H = patch('vertices',[va;vb;vc],'faces',[1,2,3],'facecolor','none','edgecolor','r','linewidth',2);

end