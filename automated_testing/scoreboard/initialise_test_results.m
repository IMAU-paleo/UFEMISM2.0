function res = initialise_test_results( test_name, test_category)
% Initialise the results of a single component/integrated test

res.name            = test_name;
res.category        = test_category;
res.date_and_time   = string( datetime);
res.git_hash_string = git_hash_string;

result_dummy.name          = '';
result_dummy.description   = '';
result_dummy.cost_function = 0;

res.results(1) = result_dummy;

end