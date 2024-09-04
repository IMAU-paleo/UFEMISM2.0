function analyse_component_and_integrated_tests( varargin)
% Analyse the results of all component tests and integrated tests, and add
% the results to the scoreboard files

% In the GitHub Workflow, provide the component test and integrated test
% results folders as input; but retain the option of running without input
% (i.e. as a script) locally, with user-defined folders

input_args = varargin;
if isempty( input_args)
  % Assume this is a local run

  foldername_component_tests = '/Users/Beren017/Documents/GitHub/UFEMISM2.0/results_component_tests';
  % foldername_integrated_tests = '';
  addpath('/Users/Beren017/Documents/GitHub/UFEMISM2.0/tools/matlab/')
  do_print_figures = false;
  foldername_automated_testing = '.';

elseif length( input_args) == 2
  % Assume this is a GitHub Workflow run

  foldername_component_tests = varargin{1};
  % foldername_integrated_tests = varargin{2};
  addpath('tools/matlab/')
  do_print_figures = true;
  foldername_automated_testing = 'automated_testing';

else
  error('need either foldername_component_tests and foldername_integrated_tests, or nothing as input!')
end

addpath('analysis_scripts')
addpath('analysis_scripts/basic')
addpath('analysis_scripts/component_tests')
addpath('analysis_scripts/integrated_tests')

analyse_component_tests( foldername_component_tests, foldername_automated_testing, do_print_figures);

end