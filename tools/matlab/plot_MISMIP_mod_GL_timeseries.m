clc
clear all
close all

foldernames = {
  '../../results_benchmark_MISMIP_mod_spinup/'
  '../../results_benchmark_MISMIP_mod_advance/'
  '../../results_benchmark_MISMIP_mod_retreat/'
  };

wa = 600;
ha = 400;

margins_hor = [90,25];
margins_ver = [25,80];

H = setup_multipanel_figure( wa, ha, margins_hor, margins_ver);

xlabel( H.Ax{ 1,1},'Time (kyr)')
ylabel( H.Ax{ 1,1},'GL position (km)')

for fi = 1: length( foldernames)
  
  filename = [foldernames{ fi} '/main_output_ANT_00001.nc'];
  
  mesh = read_mesh_from_file( filename);
  
  time = ncread( filename,'time') / 1e3;
  r_GL_av = zeros( size( time));
  
  for ti = 1: length( time)
    
    Hi = ncread( filename,'Hi',[1,ti],[Inf,1]);
    Hb = ncread( filename,'Hb',[1,ti],[Inf,1]);
    
    seawater_density = 1028;
    ice_density      = 910;
    TAF = Hi - max( 0.0, (0 - Hb) * (seawater_density / ice_density));
    
    C_GL = mesh_contour( mesh,TAF,0);
    r_GL = sqrt( C_GL( :,1).^2 + C_GL( :,2).^2) / 1e3;
    r_GL_av( ti) = mean( r_GL);
    
  end
  
  line('parent',H.Ax{ 1,1},'xdata',time,'ydata',r_GL_av,'linewidth',3,'color','b');
end