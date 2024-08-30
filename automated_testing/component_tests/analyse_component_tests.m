function analyse_component_tests( foldername)

% foldername = '../../results_component_tests';
% addpath('../../tools/matlab/')

addpath('tools/matlab/')

do_print = false;

% analyse_discretisation_tests( [foldername '/discretisation'], do_print);
analyse_remapping_tests(      [foldername '/remapping']     , do_print);

end