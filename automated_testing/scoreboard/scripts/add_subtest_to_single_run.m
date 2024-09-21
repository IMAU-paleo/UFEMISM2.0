function single_run = add_subtest_to_single_run( single_run, name, description, cost_function)
% Add a single sub-result to the results of a single component/integrated
% test run
%
% e.g.:
% 
% name          = 'max_abs_err'
% description   = 'max( abs( d_mesh - d_mesh_ex))'
% cost_function = 1e-3
%
% Note that the 'cost_function' is completely user-defined for each
% individual sub-test, and should be minimised; any cost functions that
% have increased with respect to the previous run of the test will be
% interpreted by the scoreboard as a decrease in model performance.

subtest.name          = name;
subtest.description   = description;
subtest.cost_function = cost_function;

if strcmp( single_run.subtests(1).name,'') && ...
   strcmp( single_run.subtests(1).description,'') && ...
   single_run.subtests(1).cost_function == 0
  single_run.subtests( 1) = subtest;
else
  single_run.subtests( end+1) = subtest;
end

end