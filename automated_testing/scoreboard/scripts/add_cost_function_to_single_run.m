function single_run = add_cost_function_to_single_run( single_run, name, definition, value)
% Add a new cost function to the results of a single run of a
% component/integrated test
%
% e.g.:
% 
% name       = 'max_abs_err'
% definition = 'max( abs( d_mesh - d_mesh_ex))'
% value       = 1e-3

cost_function.name       = name;
cost_function.definition = definition;
cost_function.value      = value;

if strcmp( single_run.cost_functions(1).name,'') && ...
   strcmp( single_run.cost_functions(1).definition,'') && ...
   single_run.cost_functions(1).value == 0
  single_run.cost_functions( 1) = cost_function;
else
  single_run.cost_functions( end+1) = cost_function;
end

end