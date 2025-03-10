function delete_old_scoreboard_files( varargin)
% Analyse all the scoreboard files and delete all but the 10 most recent
% ones

%%

% In the GitHub Workflow, provide the automated_testing folder as
% input; but retain the option of running without input (i.e. as a script)
% locally, with user-defined folders

input_args = varargin;
if isempty( input_args)
  % Assume this is a local run
  foldername_automated_testing = '/Users/Beren017/Documents/GitHub/UFEMISM2.0/automated_testing';
elseif isscalar( input_args)
  % Assume this is a GitHub Workflow run
  foldername_automated_testing = varargin{1};
else
  error('need either foldername_automated_testing, or nothing as input!')
end

%%

list_of_scoreboard_files = dir( [foldername_automated_testing '/scoreboard/scoreboard_files']);

i = 1;
while i <= length( list_of_scoreboard_files)
  if ~endsWith( list_of_scoreboard_files(i).name, '.xml')
    i = i+1;
    continue
  end

  run_i = read_scoreboard_file( list_of_scoreboard_files(i).name);
  run_i.filename = list_of_scoreboard_files(i).name;

  disp(['Deleting old scoreboard files for test ' char(run_i.name) '...'])

  all_runs_of_this_test = run_i;

  j = i+1;
  while j <= length( list_of_scoreboard_files)
    if ~endsWith( list_of_scoreboard_files(j).name, '.xml')
      break
    end
    run_j = read_scoreboard_file( list_of_scoreboard_files(j).name);
    run_j.filename = list_of_scoreboard_files(j).name;
    if strcmpi( run_i.name, run_j.name) && strcmpi( run_i.category, run_j.category)
      all_runs_of_this_test( end+1) = run_j;
      j = j+1;
    else
      break
    end
  end

  % Sort by date and time
  for ii = 1: length( all_runs_of_this_test)-1
    for jj = ii+1: length( all_runs_of_this_test)
      if all_runs_of_this_test(ii).date_and_time > all_runs_of_this_test(jj).date_and_time
        all_runs_of_this_test([ii,jj]) = all_runs_of_this_test([jj,ii]);
      end
    end
  end

  % Delete all but 10 most recent runs of this test
  for ii = 1: length(all_runs_of_this_test)-9
    filename = [foldername_automated_testing '/scoreboard/scoreboard_files/' all_runs_of_this_test(ii).filename];
    delete( filename)
  end

  i = j;

end

end