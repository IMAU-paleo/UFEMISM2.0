function verify_unit_tests_results(foldername)
% Verify that all the unit tests were passed.
% If not, throw an error (which, when done in a GitHub Workflow,
% will cause the job to fail, so the merge cannot be completed).

disp('Verifying unit test results...')

% foldername = '../results_unit_tests';
filename = [foldername '/unit_tests_output.txt'];

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
end

end