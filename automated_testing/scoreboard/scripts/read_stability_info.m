function stab = read_stability_info( filename, nskip)

stab.n_dt_ice   = length( ncread( filename,'dt_ice'    ,nskip,Inf));
stab.n_visc_its = sum(    ncread( filename,'n_visc_its',nskip,Inf));
stab.n_Axb_its  = sum(    ncread( filename,'n_Axb_its' ,nskip,Inf));
  
end