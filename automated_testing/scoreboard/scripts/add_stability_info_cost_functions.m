function single_run = add_stability_info_cost_functions( single_run, stab)

single_run = add_cost_function_to_single_run( single_run, ...
  'n_dt_ice', 'total number of ice-dynamical time steps', stab.n_dt_ice);
single_run = add_cost_function_to_single_run( single_run, ...
  'n_visc_its', 'total number of viscosity iterations', stab.n_visc_its);
single_run = add_cost_function_to_single_run( single_run, ...
  'n_Axb_its', 'total number of linear solver iterations', stab.n_Axb_its);

end