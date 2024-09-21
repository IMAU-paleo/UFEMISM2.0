function list_of_tests = list_all_integrated_tests( foldername_integrated_tests)

list_of_tests = list_all_integrated_tests_recursive( {}, foldername_integrated_tests);

  function list_of_tests = list_all_integrated_tests_recursive( list_of_tests, test_path)
  
    henk = dir( test_path);
  
    is_test = false;
    for i = 1: length( henk)
      if strcmpi( henk( i).name,'config.cfg')
        is_test = true;
      end
    end
    if is_test
      list_of_tests{ end+1} = test_path;%( 3:end);
      return
    end
  
    for i = 1: length( henk)
      if strcmpi( henk( i).name,'.') || strcmpi( henk( i).name, '..')
        continue
      end
      if henk( i).isdir
        list_of_tests = list_all_integrated_tests_recursive( list_of_tests, [test_path '/' henk( i).name]);
      end
    end
  
  end

end