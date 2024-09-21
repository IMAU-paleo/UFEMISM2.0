function add_all_integrated_test_results_to_main_scoreboard_files( varargin)
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

disp('Adding all integrated test results to the main scoreboard files...')
disp('')

%%

% In the GitHub Workflow, provide the component test results folders as
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

%%

list_of_tests = list_all_integrated_tests( foldername_integrated_tests);

for i = 1: length( list_of_tests)
  test_path = list_of_tests{ i};
  ii = strfind( test_path,'/');ii = ii( end);
  test_path_firstname = test_path( ii+1:end);
  disp(['  Adding results of ' test_path_firstname ' to main scoreboard file...'])
  add_integrated_test_result_to_main_scoreboard_file( test_path);
end

  function add_integrated_test_result_to_main_scoreboard_file( test_path)

    % Check if a temporary scoreboard file exists in this integrated test folder
    scoreboard_file_exists = false;
    henk = dir( test_path);
    for ii = 1: length( henk)
      if contains( henk( ii).name,'scoreboard_') && contains( henk( ii).name,'.txt')
        scoreboard_file_exists = true;
        scoreboard_filename = [test_path '/' henk( ii).name];
      end
    end
    % If not, do nothing
    if ~scoreboard_file_exists
      return
    end

    % Read single test run results from temporary scoreboard file
    single_run = read_single_run_from_temporary_scoreboard_file( scoreboard_filename);

    % Write single test run results to main scoreboard file
    write_single_test_run_to_scoreboard_file( single_run, foldername_scoreboard);

  end

end