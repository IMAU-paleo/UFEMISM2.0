function analyse_component_and_integrated_tests( foldername_component_tests, foldername_integrated_tests)
% Analyse the results of all component tests and integrated tests, and add
% the results to the scoreboard files

foldername_component_tests = '/Users/Beren017/Documents/GitHub/UFEMISM2.0/results_component_tests';
% foldername_integrated_tests = '';

addpath('tools/matlab/')
addpath('/Users/Beren017/Documents/GitHub/UFEMISM2.0/tools/matlab/')
addpath('analysis_scripts')
addpath('analysis_scripts/basic')
addpath('analysis_scripts/component_tests')
addpath('analysis_scripts/integrated_tests')

do_print_figures = false;

analyse_component_tests( foldername_component_tests, do_print_figures);

end