function mesh = read_mesh_from_file( filename)

  mesh.V              = ncread( filename,'V');
  mesh.nC             = ncread( filename,'nC');
  mesh.C              = ncread( filename,'C');
  mesh.niTri          = ncread( filename,'niTri');
  mesh.iTri           = ncread( filename,'iTri');
  mesh.VBI            = ncread( filename,'VBI');
  
  mesh.Tri            = ncread( filename,'Tri');
  mesh.TriC           = ncread( filename,'TriC');
  mesh.Tricc          = ncread( filename,'Tricc');
  mesh.TriBI          = ncread( filename,'TriBI');
  
  mesh.nV             = size( mesh.V,1);
  mesh.nTri           = size( mesh.Tri,1);
  
  mesh.xmin           = mesh.V( 1,1);
  mesh.xmax           = mesh.V( 2,1);
  mesh.ymin           = mesh.V( 2,2);
  mesh.ymax           = mesh.V( 3,2);
  
  mesh.E              = ncread(filename,'E');
  mesh.VE             = ncread(filename,'VE');
  mesh.EV             = ncread(filename,'EV');
  mesh.ETri           = ncread(filename,'ETri');
  mesh.nE             = size( mesh.E,1);
  
  mesh.R              = ncread( filename,'R');
  mesh.A              = ncread( filename,'A');
  mesh.lon            = ncread( filename,'lon');
  mesh.lat            = ncread( filename,'lat');

end