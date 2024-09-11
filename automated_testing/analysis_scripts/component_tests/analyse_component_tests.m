function analyse_component_tests( varargin)
% Analyse the results of all component tests and add the results to the
% scoreboard files.

%%

% In the GitHub Workflow, provide the component test results folders as
% input; but retain the option of running without input (i.e. as a script)
% locally, with user-defined folders

input_args = varargin;
if isempty( input_args)
  % Assume this is a local run

  foldername_component_tests = '/Users/Beren017/Documents/GitHub/UFEMISM2.0/results_component_tests';
  addpath('/Users/Beren017/Documents/GitHub/UFEMISM2.0/tools/matlab/')
  addpath('/Users/Beren017/Documents/GitHub/UFEMISM2.0/automated_testing/analysis_scripts/basic')
  do_print_figures = false;
  foldername_automated_testing = '.';

elseif length( input_args) == 1
  % Assume this is a GitHub Workflow run

  foldername_component_tests = varargin{1};
  addpath('tools/matlab/')
  addpath('automated_testing/analysis_scripts/basic')
  do_print_figures = true;
  foldername_automated_testing = 'automated_testing';

else
  error('need either foldername_component_tests, or nothing as input!')
end

%%

analyse_component_tests_discretisation( ...
  [foldername_component_tests '/discretisation'], foldername_automated_testing, do_print_figures);

analyse_component_tests_remapping( ...
  [foldername_component_tests '/remapping'], foldername_automated_testing, do_print_figures);

end