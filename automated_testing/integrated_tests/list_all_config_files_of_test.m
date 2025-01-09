function config_files = list_all_config_files_of_test( test_path)

henk = dir( test_path);

config_files = {};
for i = 1: length( henk)
  if contains( henk( i).name,'config') && contains( henk( i).name,'.cfg')
    config_files{ end+1} = henk( i).name;
  end
end

end