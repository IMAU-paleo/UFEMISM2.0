function verify_unit_tests_results( varargin)
% Verify that all the unit tests were passed.
% If not, throw an error (which, when done in a GitHub Workflow,
% will cause the job to fail, so the merge cannot be completed).

disp('Verifying unit test results...')

%%

% In the GitHub Workflow, provide the scoreboard folder as
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

filename = [foldername_automated_testing '/unit_tests/results/unit_tests_output.txt'];

fid = fopen(filename);
temp = textscan(fid,'%s','delimiter','\n'); temp = temp{1};
fclose(fid);

all_unit_tests_passed = true;
for i = 1: length( temp)
  if contains(temp{i},'Unit test passed:')
    % This unit test was passed succesfully
  elseif contains(temp{i},'Unit test failed:')
    % This unit test was failed
    all_unit_tests_passed = false;
    warning(temp{i})
  end
end

if ~all_unit_tests_passed
  disp('')
  disp('===================================================')
  disp('===== ERROR - not all unit tests were passed! =====')
  disp('===================================================')
  disp('')
  error('')
else
  disp('Unit tests passed!')
end

end