function append_test_results_to_main_scoreboard_file( all_runs_new, filename_main)

% Safety
if ~isscalar( all_runs_new.single_run)
  error('Temporary scoreboard file contains results of more than one test run!')
end

single_run_new = all_runs_new.single_run;

% If the main scoreboard file doesn't exist yet, create it.
if ~exist( filename_main, 'file')
  writestruct( all_runs_new, filename_main)
  return
end

% Read the main file too
all_runs_main = read_scoreboard_file( filename_main);

% If the git hash string of the last entry in the main file is identical to
% that of the new single run, overwrite the last run
do_overwrite_last_run = false;
if strcmpi( all_runs_main.single_run( end).git_hash_string, single_run_new.git_hash_string)
  do_overwrite_last_run = true;
end

if do_overwrite_last_run
  all_runs_main.single_run( end  ) = single_run_new;
else
  all_runs_main.single_run( end+1) = single_run_new;
end

% Delete the old main scoreboard file
delete( filename_main);

% Create a new main scoreboard file with the new single run included
write_scoreboard_file( all_runs_main, filename_main);

end