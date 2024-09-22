function write_scoreboard_file( all_runs, filename)
% Write a structure containing the results of all the runs of a
% component/integrated test to a scoreboard file

% Only keep the last 50 test runs
n = length( all_runs.single_run);
i = max( 1, n - 49);
all_runs.single_run = all_runs.single_run( i:n);

% Write
writestruct( all_runs, filename, 'StructNodeName', 'all_runs');

end