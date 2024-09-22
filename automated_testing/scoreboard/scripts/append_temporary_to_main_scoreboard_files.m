function append_temporary_to_main_scoreboard_files( varargin)
% We want to be able to run each individual integrated test as a separate
% job in the GitHub workflow, so that they will be executed (and analysed)
% in parallel. However, that means that the test results (both the NetCDF
% files produced by UFEMISM, and the scoreboard files created by the Matlab
% analysis script) will exist only within that job. Committing the changes
% in the scoreboard files from those jobs would mean creating a huge amount
% of commits, which will mess up the repo.
% So instead, we create a small, temporary scoreboard file within the job,
% and upload that file as a job artifact. Once all integrated tests have
% finished, the main job will download all those artifacts, add the test
% results to the "real" scoreboard files, and then commit all the changes
% together.
% A little bit convoluted, perhaps, but much cleaner overall!

disp('Appending results from temporary scoreboard files to main scoreboard files...')
disp('')

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

foldername_integrated_tests = [foldername_automated_testing '/integrated_tests'];
foldername_scoreboard       = [foldername_automated_testing '/scoreboard'];

addpath([foldername_scoreboard '/scripts'])

%% Component tests

list_of_temporary_scoreboard_files = dir( [foldername_automated_testing ...
  '/component_tests/temporary_scoreboard_files']);

for i = 1: length( list_of_temporary_scoreboard_files)

  filename_short = list_of_temporary_scoreboard_files( i).name;
  
  if contains( filename_short,'.xml')

    disp(['  Adding results of ' filename_short ' to main scoreboard file...'])

    filename_temporary_scoreboard_file = [foldername_automated_testing ...
      '/component_tests/temporary_scoreboard_files/' filename_short];

    % Read the single new test run from the temporary scoreboard file
    all_runs_temp = read_scoreboard_file( filename_temporary_scoreboard_file);

    % Append the results of the single run from the temporary scoreboard file to the main scoreboard file
    filename_main_scoreboard_file = [foldername_scoreboard '/scoreboard_files/' filename_short];
    append_test_results_to_main_scoreboard_file( all_runs_temp, filename_main_scoreboard_file);

  end
end

%% Integrated tests

list_of_tests = list_all_integrated_tests( foldername_integrated_tests);

for i = 1: length( list_of_tests)
  test_path = list_of_tests{ i};
  add_integrated_test_result_to_main_scoreboard_file( test_path);
end

  function add_integrated_test_result_to_main_scoreboard_file( test_path)

    % Check if a temporary scoreboard file exists in this integrated test folder
    scoreboard_file_exists = false;
    henk = dir( test_path);
    for ii = 1: length( henk)
      if contains( henk( ii).name,'scoreboard_temp_') && contains( henk( ii).name,'.xml')
        scoreboard_file_exists = true;
        filename_temporary_scoreboard_file = henk( ii).name;
      end
    end
    % If not, do nothing
    if ~scoreboard_file_exists
      return
    end

    ii = strfind( test_path,'/');ii = ii( end);
    test_path_firstname = test_path( ii+1:end);
    disp(['  Adding results of ' test_path_firstname ' to main scoreboard file...'])

    % Read the single new test run from the temporary scoreboard file
    all_runs_temp = read_scoreboard_file( [test_path '/' filename_temporary_scoreboard_file]);

    % Append the results of the single run from the temporary scoreboard file to the main scoreboard file
    filename_main_scoreboard_file = [strrep( filename_temporary_scoreboard_file, ...
      'scoreboard_temp_', '')];
    append_test_results_to_main_scoreboard_file( ...
      all_runs_temp, [foldername_scoreboard '/scoreboard_files/' filename_main_scoreboard_file]);

  end

end