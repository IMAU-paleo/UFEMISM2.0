function mesh = read_mesh_from_text_file( filename)

fid = fopen(filename);
temp = textscan(fid,'%s %f',7,'headerlines',1,'delimiter',{'=','\n'});

mesh.xmin    = temp{2}(1);
mesh.xmax    = temp{2}(2);
mesh.ymin    = temp{2}(3);
mesh.ymax    = temp{2}(4);
mesh.nC_mem  = temp{2}(5);
mesh.nV      = temp{2}(6);
mesh.nTri    = temp{2}(7);

temp = textscan(fid,'%s',4,'delimiter','\n');

% Vertex data
ncols = 2 + 1 + mesh.nC_mem + 1 + mesh.nC_mem + 1;
formatspec = '%f %f';
for n=3:ncols
  formatspec = [formatspec ' %d'];
end
temp = textscan(fid,formatspec,mesh.nV);

nc = 1;
mesh.V = [temp{1} temp{2}];
nc = nc+2;

mesh.nC = temp{nc};
nc = nc+1;
mesh.C = temp{nc}; for c = 2:mesh.nC_mem; nc = nc+1; mesh.C = [mesh.C temp{nc}]; end
nc = nc+1;

mesh.niTri = temp{nc};
nc = nc+1;
mesh.iTri = temp{nc}; for c = 2:mesh.nC_mem; nc = nc+1; mesh.iTri = [mesh.iTri temp{nc}]; end
nc = nc+1;

mesh.edge_index = temp{nc};

temp = textscan(fid,'%s',3,'delimiter','\n');

% Triangle data
temp = textscan(fid,'%d %d %d %d %d %d %d',mesh.nTri);

mesh.Tri  = [temp{1} temp{2} temp{3}];
mesh.TriC = [temp{4} temp{5} temp{6}];
mesh.Tri_edge_index = temp{7};

fclose(fid);

end