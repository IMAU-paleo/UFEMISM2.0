function stab = read_stability_info( filename, nskip)

stab.n_dt_ice   = length( ncread( filename,'dt_ice'    ,nskip+1,Inf));
stab.n_visc_its = sum(    ncread( filename,'n_visc_its',nskip+1,Inf));
stab.n_Axb_its  = sum(    ncread( filename,'n_Axb_its' ,nskip+1,Inf));
  
end