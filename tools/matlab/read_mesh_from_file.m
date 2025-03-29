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
  mesh.nC_mem         = size( mesh.C,2);
  
  mesh.xmin           = min( mesh.V( :,1));
  mesh.xmax           = max( mesh.V( :,1));
  mesh.ymin           = min( mesh.V( :,2));
  mesh.ymax           = max( mesh.V( :,2));
  
  mesh.E              = ncread(filename,'E');
  mesh.VE             = ncread(filename,'VE');
  mesh.EV             = ncread(filename,'EV');
  mesh.ETri           = ncread(filename,'ETri');
  mesh.TriE           = ncread(filename,'TriE');
  mesh.EBI            = ncread(filename,'EBI');
  mesh.nE             = size( mesh.E,1);

  mesh.vi2vori        = ncread( filename, 'vi2vori');
  mesh.ti2vori        = ncread( filename, 'ti2vori');
  mesh.ei2vori        = ncread( filename, 'ei2vori');
  mesh.vori2vi        = ncread( filename, 'vori2vi');
  mesh.vori2ti        = ncread( filename, 'vori2ti');
  mesh.vori2ei        = ncread( filename, 'vori2ei');
  mesh.Vor            = ncread( filename, 'Vor');
  mesh.VornC          = ncread( filename, 'VornC');
  mesh.VorC           = ncread( filename, 'VorC');
  mesh.nVVor          = ncread( filename, 'nVVor');
  mesh.VVor           = ncread( filename, 'VVor');
  
  mesh.TriGC          = ncread( filename,'TriGC');
  mesh.TriA           = ncread( filename,'TriA');
  mesh.R              = ncread( filename,'R');
  mesh.A              = ncread( filename,'A');
  mesh.lon            = ncread( filename,'lon');
  mesh.lat            = ncread( filename,'lat');

end