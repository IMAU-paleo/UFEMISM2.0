function single_run = initialise_single_test_run( test_name, test_category)
% Initialise the results of a single run of a component/integrated test

single_run.name            = test_name;
single_run.category        = test_category;
single_run.date_and_time   = string( datetime);
single_run.git_hash_string = git_hash_string;

cost_function_dummy.name       = '';
cost_function_dummy.definition = '';
cost_function_dummy.value      = 0;

single_run.cost_functions(1) = cost_function_dummy;

end