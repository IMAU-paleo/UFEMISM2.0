function res = add_result_to_test_results( res, name, description, cost_function)
% Add a single result to the results of a single component/integrated test
%
% e.g.:
% 
% name          = 'max_abs_err'
% description   = 'max( abs( d_mesh - d_mesh_ex))'
% cost_function = 1e-3
%
% Note that the 'cost_function' is completely user-defined for each
% individual test, and should be minimised; any cost functions that have
% increased with respect to the previous run of the test will be
% interpreted by the scoreboard as a decrease in model performance.

result.name          = name;
result.description   = description;
result.cost_function = cost_function;

if strcmp( res.results(1).name,'') && ...
   strcmp( res.results(1).description,'') && ...
   res.results(1).cost_function == 0
  res.results( 1) = result;
else
  res.results( end+1) = result;
end

end