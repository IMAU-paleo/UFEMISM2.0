function single_run = initialise_single_test_run( test_name, test_category)
% Initialise the results of a single component/integrated test run

single_run.name            = test_name;
single_run.category        = test_category;
single_run.date_and_time   = string( datetime);
single_run.git_hash_string = git_hash_string;

subtest_dummy.name          = '';
subtest_dummy.description   = '';
subtest_dummy.cost_function = 0;

single_run.subtests(1) = subtest_dummy;

end