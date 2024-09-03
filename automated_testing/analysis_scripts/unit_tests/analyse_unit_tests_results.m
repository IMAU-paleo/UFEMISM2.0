function analyse_unit_tests_results(foldername)
% Process the UFEMISM unit tests results into a nice, semi-interactive html report.

disp('Analysing unit test results...')

% foldername = '../../results_unit_tests';
filename = [foldername '/unit_tests_output.txt'];

R = read_unit_tests_structure( filename);

filename_html = [foldername '/unit_tests_report.html'];
create_unit_tests_report_html( R, filename_html);

  function R = read_unit_tests_structure( filename)
    % Read the UFEMISM unit tests output file and collect
    % its results into a single nested Matlab structure.

    disp(['Reading unit tests results from file "' filename '"...'])

    fid = fopen(filename);
    temp = textscan(fid,'%s','delimiter','\n'); temp = temp{1};
    fclose(fid);

    R.UFEMISM = [];
    
    for i = 1: length( temp)
      R = process_line( R, temp{i});
    end

    R = gather_test_results_up_tree( R);

    function R = process_line( R, line_full)

      thing = 'Unit test passed: ';
      test_result_str = line_full( 1:length(thing));
      line            = line_full( length(thing)+1:end);

      if strcmpi( test_result_str, 'Unit test passed: ')
        test_result = true;
      elseif strcmpi( test_result_str, 'Unit test failed: ')
        test_result = false;
      else
        error(['Unrecognised unit test result: "' test_result_str '"'])
      end

      if contains( line, '/')
        % Move down the unit test family tree

        ii = strfind( line, '/'); ii = ii(1);
        name     = line( 1   :ii-1);
        children = line( ii+1:end );

        % Check if this name is already listed in R. If so, move down a level.
        found_match = false;
        if ~isempty(R)
          f = fields(R);
          for fi = 1: length( f)
            if strcmpi( f{fi},name)
              % Found a match; move down a level
              found_match = true;
              R.(f{fi}) = process_line( R.(f{fi}), [test_result_str children]);
            end
          end
        end
        % If not, add it.
        if ~found_match
          R.(name) = [];
          R.(name) = process_line( R.(name), [test_result_str children]);
        end

      else
        % We've reached an end of the tree

        R.(line) = test_result;

      end

    end
    function R = gather_test_results_up_tree( R)

      test_result = true;

      f = fields(R);
      for fi = 1: length( f)
        if islogical( R.(f{fi}))
          test_result = test_result && R.(f{fi});
        elseif isstruct( R.(f{fi}))
          R.(f{fi}) = gather_test_results_up_tree( R.(f{fi}));
          test_result = test_result && R.(f{fi}).test_result;
        end
      end

      R.test_result = test_result;

    end

  end

  function create_unit_tests_report_html( R, filename)
    % Create a nice, semi-interactive html file where the
    % user can easily inspect the UFEMISM unit test results.

    disp(['Creating unit tests report in file "' filename '"...'])

    fid = fopen(filename,'w');

    fprintf(fid,'<!DOCTYPE html>\n');
    fprintf(fid,'<html>\n');
    fprintf(fid,'\n');
    fprintf(fid,'<!-- -->\n');
    fprintf(fid,'<!-- UFEMISM unit tests report -->\n');
    fprintf(fid,'<!-- -->\n');
    fprintf(fid,'\n');

    print_html_head( fid)

    fprintf(fid,'<body>\n');
    fprintf(fid,'\n');
    fprintf(fid,'<div class="container">\n');
    fprintf(fid,'<h1>UFEMISM unit tests report</h1>\n');

    process_unit_tests_tree( R.UFEMISM, fid, 0);

    fprintf(fid,'</div>\n');
    fprintf(fid,'\n');

    print_html_tail_script( fid);

    fprintf(fid,'\n');
    fprintf(fid,'</body>\n');
    fprintf(fid,'</html>\n');

    fclose(fid);

    function print_html_head( fid)
      fprintf(fid,'<head>\n');
      fprintf(fid,'<meta name="viewport" content="width=device-width, initial-scale=1">\n');
      fprintf(fid,'<style>\n');
      fprintf(fid,'.collapsible_pass {\n');
      fprintf(fid,'  background-color: #00FF00;\n');
      fprintf(fid,'  color: black;\n');
      fprintf(fid,'  cursor: pointer;\n');
      fprintf(fid,'  padding: 8px;\n');
      fprintf(fid,'  width: 100%;\n');
      fprintf(fid,'  border: none;\n');
      fprintf(fid,'  text-align: left;\n');
      fprintf(fid,'  outline: none;\n');
      fprintf(fid,'  font-size: 20px;\n');
      fprintf(fid,'}\n');
      fprintf(fid,'\n');
      fprintf(fid,'.collapsible_fail {\n');
      fprintf(fid,'  background-color: #FF6666;\n');
      fprintf(fid,'  color: black;\n');
      fprintf(fid,'  cursor: pointer;\n');
      fprintf(fid,'  padding: 8px;\n');
      fprintf(fid,'  width: 100%;\n');
      fprintf(fid,'  border: none;\n');
      fprintf(fid,'  text-align: left;\n');
      fprintf(fid,'  outline: none;\n');
      fprintf(fid,'  font-size: 20px;\n');
      fprintf(fid,'}\n');
      fprintf(fid,'\n');
      fprintf(fid,'.active {\n');
      fprintf(fid,'  border: solid 3px black;\n');
      fprintf(fid,'}\n');
      fprintf(fid,'\n');
      fprintf(fid,'.content {\n');
      fprintf(fid,'  padding: 0 0 0 50px;\n');
      fprintf(fid,'  display: none;\n');
      fprintf(fid,'  overflow: hidden;\n');
      fprintf(fid,'  background-color: white;\n');
      fprintf(fid,'  border: solid 2px black;\n');
      fprintf(fid,'  line-height: 1.0;\n');
      fprintf(fid,'}\n');
      fprintf(fid,'</style>\n');
      fprintf(fid,'</head>\n');
    end
    function print_html_tail_script( fid)
      fprintf(fid,'<script>\n');
      fprintf(fid,'\n');
      fprintf(fid,'var coll = document.getElementsByClassName("collapsible_pass");\n');
      fprintf(fid,'var i;\n');
      fprintf(fid,'\n');
      fprintf(fid,'for (i = 0; i < coll.length; i++) {\n');
      fprintf(fid,'  coll[i].addEventListener("click", function() {\n');
      fprintf(fid,'    this.classList.toggle("active");\n');
      fprintf(fid,'    var content = this.nextElementSibling;\n');
      fprintf(fid,'    if (content.style.display === "block") {\n');
      fprintf(fid,'      content.style.display = "none";\n');
      fprintf(fid,'    } else {\n');
      fprintf(fid,'      content.style.display = "block";\n');
      fprintf(fid,'    }\n');
      fprintf(fid,'  });\n');
      fprintf(fid,'}\n');
      fprintf(fid,'var coll = document.getElementsByClassName("collapsible_fail");\n');
      fprintf(fid,'var i;\n');
      fprintf(fid,'\n');
      fprintf(fid,'for (i = 0; i < coll.length; i++) {\n');
      fprintf(fid,'  coll[i].addEventListener("click", function() {\n');
      fprintf(fid,'    this.classList.toggle("active");\n');
      fprintf(fid,'    var content = this.nextElementSibling;\n');
      fprintf(fid,'    if (content.style.display === "block") {\n');
      fprintf(fid,'      content.style.display = "none";\n');
      fprintf(fid,'    } else {\n');
      fprintf(fid,'      content.style.display = "block";\n');
      fprintf(fid,'    }\n');
      fprintf(fid,'  });\n');
      fprintf(fid,'}\n');
      fprintf(fid,'\n');
      fprintf(fid,'\n');
      fprintf(fid,'</script>\n');
    end
    function process_unit_tests_tree( R, fid, depth)

      is_leaf = true;
      f = fields(R);
      for fi = 1: length(f)
        if isstruct(R.(f{fi}))
          is_leaf = false;
        end
      end

      if ~is_leaf
        process_unit_tests_tree_branch( R, fid, depth)
      else
        process_unit_tests_tree_leaf( R, fid, depth)
      end

    end
    function process_unit_tests_tree_branch( R, fid, depth)

      f = fields(R);
      for fi = 1: length(f)
        if strcmpi(f{fi},'test_result'); continue; end

        if R.(f{fi}).test_result
          div_class = 'collapsible_pass';
        else
          div_class = 'collapsible_fail';
        end

        for i = 1: depth
          fprintf(fid,' ');
        end
        fprintf(fid,'<div>\n');
        for i = 1: depth
          fprintf(fid,' ');
        end
        fprintf(fid,['<button type="button" class="' div_class '">' f{fi} '</button>\n']);
        for i = 1: depth
          fprintf(fid,' ');
        end
        if R.(f{fi}).test_result
          fprintf(fid,'<div class="content">\n');
        else
          fprintf(fid,'<div class="content" style="display:block">\n');
        end
        process_unit_tests_tree( R.(f{fi}), fid, depth+2)
        for i = 1: depth
          fprintf(fid,' ');
        end
        fprintf(fid,'</div>\n');
        for i = 1: depth
          fprintf(fid,' ');
        end
        fprintf(fid,'</div>\n');

      end

    end
    function process_unit_tests_tree_leaf( R, fid, depth)
      f = fields(R);
      for fi = 1: length(f)
        if strcmpi(f{fi},'test_result'); continue; end
        for i = 1: depth
          fprintf(fid,' ');
        end
        if R.(f{fi})
          % fprintf(fid,['<p style="color: #006600; font-weight: bold; font-size: 14pt">' f{fi} '</p>\n']);
          fprintf(fid,['<div style="background-color: #00FF00; font-weight: bold; font-size: 16pt">' f{fi} '</div>\n']);
        else
          fprintf(fid,['<div style="background-color: #FF6666; font-weight: bold; font-size: 16pt">' f{fi} '</div>\n']);
        end
      end
    end

  end

end