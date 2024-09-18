function analyse_unit_tests_results( varargin)
% Process the UFEMISM unit tests results into a nice, semi-interactive html report.

disp('Analysing unit test results...')

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
R = read_unit_tests_structure( filename);

filename_html = [foldername_automated_testing '/test_reports/unit_tests_report.html'];
create_unit_tests_report_html( R, filename_html);

%%

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

    fprintf( fid, '%s\n', '<!DOCTYPE html>');
    fprintf( fid, '%s\n', '<html>');
    fprintf( fid, '%s\n', '');
    fprintf( fid, '%s\n', '<!-- -->');
    fprintf( fid, '%s\n', '<!-- UFEMISM unit tests report -->');
    fprintf( fid, '%s\n', '<!-- -->');
    fprintf( fid, '%s\n', ['<!-- Created: ' char(datetime) '-->']);
    fprintf( fid, '%s\n', '<!-- -->');
    fprintf( fid, '%s\n', '');

    print_html_head( fid)

    fprintf( fid, '%s\n', '<body>');
    fprintf( fid, '%s\n', '');
    fprintf( fid, '%s\n', '<div class="container">');
    fprintf( fid, '%s\n', '<h1>UFEMISM unit tests report</h1>');
    fprintf( fid, '%s\n', ['  <p style="font-size: 18pt;">Created: ' char(datetime) '</p>']);
    fprintf( fid, '%s\n', '  <p style="font-size: 18pt;">&#128994 = pass, &#128992 = fail</p>');

    process_unit_tests_tree( R.UFEMISM, fid, 0);

    fprintf( fid, '%s\n', '</div>');
    fprintf( fid, '%s\n', '');

    print_html_tail_script( fid);

    fprintf( fid, '%s\n', '');
    fprintf( fid, '%s\n', '</body>');
    fprintf( fid, '%s\n', '</html>');

    fclose(fid);

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

    p = '';
    for ii = 1: depth
      p = [p ' '];
    end

    f = fields(R);
    for fi = 1: length(f)
      if strcmpi(f{fi},'test_result'); continue; end

      fprintf( fid, [p '<div>\n']);
      
      if R.(f{fi}).test_result
        str = ['&#128994 - ' f{fi}];
      else
        str = ['&#128992 - ' f{fi}];
      end

      fprintf( fid, [p ' <button type="button" class="collapsible" style="font-size:18pt;">' str '</button>\n']);
      fprintf( fid, [p ' <div class="content">\n']);

      process_unit_tests_tree( R.(f{fi}), fid, depth+2)

      fprintf( fid, [p ' </div>\n']);
      fprintf( fid, [p '</div>\n']);

    end

  end
  function process_unit_tests_tree_leaf( R, fid, depth)

    p = '';
    for ii = 1: depth
      p = [p ' '];
    end

    f = fields(R);
    for fi = 1: length(f)
      if strcmpi(f{fi},'test_result'); continue; end
      
      if R.(f{fi})
        str = ['&#128994 - ' f{fi}];
      else
        str = ['&#128992 - ' f{fi}];
      end

      fprintf(fid,[p '<div style="font-size: 16pt">' str '</div>\n']);
    end
  end

  function print_html_head( fid)
    fprintf( fid, '%s\n', '<head>');
    fprintf( fid, '%s\n', '<meta name="viewport" content="width=device-width, initial-scale=1">');
    fprintf( fid, '%s\n', '<style>');
    fprintf( fid, '%s\n', '.content {');
    fprintf( fid, '%s\n', '  padding: 0 0 0 50px;');
    fprintf( fid, '%s\n', '  display: none;');
    fprintf( fid, '%s\n', '  overflow: hidden;');
    fprintf( fid, '%s\n', '  background-color: white;');
    fprintf( fid, '%s\n', '  border: solid 2px black;');
    fprintf( fid, '%s\n', '  line-height: 1.0;');
    fprintf( fid, '%s\n', '}');
    fprintf( fid, '%s\n', '</style>');
    fprintf( fid, '%s\n', '</head>');
  end
  function print_html_tail_script( fid)
    fprintf( fid, '%s\n', '<script>');
    fprintf( fid, '%s\n', '');
    fprintf( fid, '%s\n', 'var coll = document.getElementsByClassName("collapsible");');
    fprintf( fid, '%s\n', 'var i;');
    fprintf( fid, '%s\n', '');
    fprintf( fid, '%s\n', 'for (i = 0; i < coll.length; i++) {');
    fprintf( fid, '%s\n', '  coll[i].addEventListener("click", function() {');
    fprintf( fid, '%s\n', '    this.classList.toggle("active");');
    fprintf( fid, '%s\n', '    var content = this.nextElementSibling;');
    fprintf( fid, '%s\n', '    if (content.style.display === "block") {');
    fprintf( fid, '%s\n', '      content.style.display = "none";');
    fprintf( fid, '%s\n', '    } else {');
    fprintf( fid, '%s\n', '      content.style.display = "block";');
    fprintf( fid, '%s\n', '    }');
    fprintf( fid, '%s\n', '  });');
    fprintf( fid, '%s\n', '}');
    fprintf( fid, '%s\n', '');
    fprintf( fid, '%s\n', '</script>');
  end

end