function single_run = read_single_run_from_temporary_scoreboard_file( scoreboard_filename)

all_runs_of_test = read_all_runs_of_test_from_scoreboard_file( scoreboard_filename);
% Safety
if all_runs_of_test.n_runs > 1
  error('whaa!')
end

% Convert to single test run structure
single_run = initialise_single_test_run( all_runs_of_test.name, all_runs_of_test.category);
single_run.date_and_time = all_runs_of_test.dates_and_times{ 1};
single_run.git_hash_string = all_runs_of_test.git_hash_strings{ 1};

for i = 1: length( all_runs_of_test.subtests)
  single_run = add_subtest_to_single_run( single_run, ...
    all_runs_of_test.subtests( i).name, ...
    all_runs_of_test.subtests( i).description, ...
    all_runs_of_test.subtests( i).cost_function);
end

end