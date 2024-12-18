function write_scoreboard_file( foldername_automated_testing, single_run)
% Write a structure containing the results of all the runs of a
% component/integrated test to a scoreboard file

% Abbreviations to shorten the filenames. These do not affect the
% scoreboard itself in any way.
string_replacements_test_category = {...
  '/',                        '_';...
  'component_tests',          'ct';...
  'integrated_tests',         'it';...
  'discretisation',           'disc';...
  'mapping_and_derivatives',  'map_deriv';...
  'remapping',                'remap';...
  'mesh_to_grid',             'm2g';...
  'grid_to_mesh',             'g2m';...
  'mesh_to_mesh',             'm2m';...
  'idealised',                'ideal';...
  'Halfar',                   'Hlf'};

replace_these = cell(1,size(string_replacements_test_category,1));
replace_with  = cell(1,size(string_replacements_test_category,1));
for i = 1: length( string_replacements_test_category)
  replace_these{i} = string_replacements_test_category{i,1};
  replace_with{ i} = string_replacements_test_category{i,2};
end

filename = [foldername_automated_testing ...
  '/scoreboard/temporary_scoreboard_files/' ...
  replace( single_run.category, replace_these, replace_with) ...
  '_' single_run.name '_' single_run.git_hash_string '.xml'];

% Write
writestruct( single_run, filename, 'StructNodeName', 'single_run');

end