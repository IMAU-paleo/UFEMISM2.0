function stab = read_stability_info( filename, nskip)
  stab.n_dt_ice                 = sum(  ncread( filename,'n_dt_ice'                ,nskip,Inf));
  stab.min_dt_ice               = min(  ncread( filename,'min_dt_ice'              ,nskip,Inf));
  stab.max_dt_ice               = max(  ncread( filename,'max_dt_ice'              ,nskip,Inf));
  stab.mean_dt_ice              = mean( ncread( filename,'mean_dt_ice'             ,nskip,Inf));
  
  stab.n_visc_its               = sum(  ncread( filename,'n_visc_its'              ,nskip,Inf));
  stab.min_visc_its_per_dt      = min(  ncread( filename,'min_visc_its_per_dt'     ,nskip,Inf));
  stab.max_visc_its_per_dt      = max(  ncread( filename,'max_visc_its_per_dt'     ,nskip,Inf));
  stab.mean_visc_its_per_dt     = mean( ncread( filename,'mean_visc_its_per_dt'    ,nskip,Inf));
  
  stab.n_Axb_its                = sum(  ncread( filename,'n_Axb_its'               ,nskip,Inf));
  stab.min_Axb_its_per_visc_it  = min(  ncread( filename,'min_Axb_its_per_visc_it' ,nskip,Inf));
  stab.max_Axb_its_per_visc_it  = max(  ncread( filename,'max_Axb_its_per_visc_it' ,nskip,Inf));
  stab.mean_Axb_its_per_visc_it = mean( ncread( filename,'mean_Axb_its_per_visc_it',nskip,Inf));
end