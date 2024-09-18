function create_scoreboard_html( varargin)
% Analyse all the scoreboard files and create the big, interactive
% scoreboard HTML file.

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

scoreboard = read_scoreboard_files( [foldername_automated_testing '/scoreboard']);

%% Write to HTML

filename = [foldername_automated_testing '/test_reports/scoreboard.html'];
if exist( filename,'file')
  delete( filename)
end

disp(['Creating scoreboard visualisation in file "' filename '"...'])

fid = fopen( filename,'w');

fprintf( fid, '%s\n', '<!DOCTYPE html>');
fprintf( fid, '%s\n', '<html>');
fprintf( fid, '%s\n', '');
fprintf( fid, '%s\n', '<!-- -->');
fprintf( fid, '%s\n', '<!-- UFEMISM automated testing scoreboard -->');
fprintf( fid, '%s\n', '<!-- -->');
fprintf( fid, '%s\n', '');

print_html_head( fid)

fprintf( fid, '%s\n', '<body >');
fprintf( fid, '%s\n', '');
fprintf( fid, '%s\n', ' <div>');
fprintf( fid, '%s\n', '  <h1>UFEMISM automated testing scoreboard</h1>');
fprintf( fid, '%s\n', ['  <p style="font-size: 18pt;">Created: ' char(datetime) '</p>']);
fprintf( fid, '%s\n', '  <p style="font-size: 18pt;">&#128993 = unchanged, &#128994 = better, &#128992 = worse</p>');

fprintf( fid, '%s\n', '  <div style="border: solid 0px #f00;">');
process_scoreboard_branch( scoreboard, fid, 0, '')
fprintf( fid, '%s\n', '  </div>');

fprintf( fid, '%s\n', ' </div>');
fprintf( fid, '%s\n', '');

print_html_tail_script( fid);

fprintf( fid, '%s\n', '');
fprintf( fid, '%s\n', '</body>');
fprintf( fid, '%s\n', '</html>');

fclose( fid);

%%

  function scoreboard = read_scoreboard_files( foldername_scoreboard)
    % Read all the scoreboard files

    list_of_files = dir( foldername_scoreboard);
    i = 1;
    while i <= length( list_of_files)
      if contains( list_of_files( i).name,'scoreboard_') && ...
         contains( list_of_files( i).name,'.txt')
        i = i+1;
      else
        list_of_files( i) = [];
      end
    end

    for i = 1: length( list_of_files)
      single_test_results(i) = read_scoreboard_file( [foldername_scoreboard '/' list_of_files( i).name]);
    end

    scoreboard = categorise_test_results( single_test_results);

  end
  function single_test_results = read_scoreboard_file( filename)
    % Read the results of a single test from its scoreboard file

    k = strfind( filename, '/');
    filename_short = filename( k(end)+1:end);
    disp(['  Reading scoreboard file ' filename_short '...'])

    % Read text from file
    fid = fopen( filename);
    complete_text = textscan( fid,'%s','delimiter','\n'); complete_text = complete_text{1};
    fclose( fid);

    test_runs_text = separate_complete_text_into_individual_runs( complete_text);
    test_runs_cellstruct = test_runs_text2cellstruct( test_runs_text);
    single_test_results = test_runs_cellstruct2struct( test_runs_cellstruct);

  end
  function test_runs_text = separate_complete_text_into_individual_runs( complete_text)
    % Separate cell array with the complete text of a scoreboard file
    % into separate cell arrays for each individual test run.

    % Determine number of test runs
    n_runs = 0;
    for i = 1: length( complete_text)
      single_line = complete_text{ i};
      if contains( single_line,'<test>')
        n_runs = n_runs + 1;
      end
    end

    test_runs_text = cell( n_runs,1);

    i_run = 0;
    while ~isempty( complete_text)

      i_start = 1;
      found_end = false;
      for i_end = 2: length( complete_text)
        single_line = complete_text{ i_end};
        if contains( single_line,'</test>')
          found_end = true;
          break
        end
      end
      % Safety
      if ~found_end
        error('whaa!')
      end

      i_run = i_run + 1;
      test_runs_text{ i_run} = complete_text( i_start: i_end);
      complete_text( i_start: i_end) = [];

    end

  end
  function test_runs_cellstruct = test_runs_text2cellstruct( test_runs_text)
    % Create a cell array with in each cell a struct containing the results of a single test run
    test_runs_cellstruct = cell(0);
    for i_run = 1: length( test_runs_text)
      test_runs_cellstruct{ end+1} = test_run_cell2struct( test_runs_text{ i_run});
    end
  end
  function test_run_struct = test_run_cell2struct( test_run_text)
    % Create a struct containing the results of a single test run

    % Safety
    if ~contains( test_run_text{ 1},'<test>') || ...
       ~contains( test_run_text{ end},'</test>')
      error('whaa!')
    end

    test_run_struct.results = [];

    is_header = true;
    for i = 1: length( test_run_text)

      single_line = test_run_text{ i};

      if contains( single_line,'<result>')
        is_header = false;
      end

      if is_header
        % Header info
        if contains( single_line,'<test>')
        elseif contains( single_line,'name')
          test_run_struct.name            = single_line( strfind( single_line, ': ')+2 : length( single_line));
        elseif contains( single_line,'category')
          test_run_struct.category        = single_line( strfind( single_line, ': ')+2 : length( single_line));
        elseif contains( single_line,'git hash string')
          test_run_struct.git_hash_string = single_line( strfind( single_line, ': ')+2 : length( single_line));
        elseif contains( single_line,'date, time')
          test_run_struct.date_and_time   = single_line( strfind( single_line, ': ')+2 : length( single_line));
        else
          error(['unknown test run header line: ' single_line])
        end
      end

      if contains( single_line,'<result>')
        % Result info

        is_header = false;

        found_end = false;
        for i_end = i+1: length( test_run_text)
          single_line = test_run_text{ i_end};
          if contains( single_line,'</result>')
            found_end = true;
            break
          end
        end
        % Safety
        if ~found_end
          error('whaa!')
        end

        res = [];

        for j = i+1: i_end-1
          single_line = test_run_text{ j};
          if contains( single_line,'name')
            res.name = single_line( strfind( single_line, ': ')+2 : length( single_line));
          elseif contains( single_line,'description')
            res.description = single_line( strfind( single_line, ': ')+2 : length( single_line));
          elseif contains( single_line,'cost function')
            res.cost_function = str2double( single_line( strfind( single_line, ': ')+2 : length( single_line)));
          else
            error(['unknown test run result line: ' single_line])
          end
        end

        if isempty( test_run_struct.results)
          test_run_struct.results = res;
        else
          test_run_struct.results( end+1) = res;
        end

      end

    end

  end
  function test_runs_struct = test_runs_cellstruct2struct( test_runs_cellstruct)
    % Create a single struct with the cost functions for each sub-test for
    % each commit

    %% Header info

    n_runs = length( test_runs_cellstruct);
    test_runs_struct.n_runs = n_runs;

    test_runs_struct.name = test_runs_cellstruct{1}.name;
    % Safety
    for i_run = 2: n_runs
      if ~strcmpi( test_runs_cellstruct{ i_run}.name, test_runs_struct.name)
        error(['Corrupt scoreboard file - test name changed between commits ' ...
          test_runs_cellstruct{ i_run-1}.git_hash_string ' and ' ...
          test_runs_cellstruct{ i_run  }.git_hash_string])
      end
    end

    test_runs_struct.category = test_runs_cellstruct{1}.category;
    % Safety
    for i_run = 2: n_runs
      if ~strcmpi( test_runs_cellstruct{ i_run}.category, test_runs_struct.category)
        error(['Corrupt scoreboard file - test category changed between commits ' ...
          test_runs_cellstruct{ i_run-1}.git_hash_string ' and ' ...
          test_runs_cellstruct{ i_run  }.git_hash_string])
      end
    end

    %% Date-and-times
    test_runs_struct.dates_and_times = cell( n_runs,1);
    for i_run = 1: n_runs
      test_runs_struct.dates_and_times{ i_run} = test_runs_cellstruct{ i_run}.date_and_time;
    end
    % Safety
    for i_run = 2: n_runs
      if datetime( test_runs_struct.dates_and_times{ i_run}) <= ...
         datetime( test_runs_struct.dates_and_times{ i_run-1})
        error(['Corrupt scoreboard file - commit ' test_runs_cellstruct{ i_run}.git_hash_string ...
          ' has a date-and-time before previous commit ' test_runs_cellstruct{ i_run-1}.git_hash_string])
      end
    end

    %% Commits
    test_runs_struct.git_hash_strings = cell( n_runs,1);
    for i_run = 1: n_runs
      test_runs_struct.git_hash_strings{ i_run} = test_runs_cellstruct{ i_run}.git_hash_string;
    end

    %% Results
    test_runs_struct.results = [];

    for i_run = 1: n_runs
      single_run = test_runs_cellstruct{ i_run};
      for i_res = 1: length( single_run.results)
        single_result = single_run.results( i_res);

        % If this result is already listed, add its value. If not, add a
        % new result
        result_exists = false;
        for j_res = 1: length( test_runs_struct.results)
          if strcmpi( test_runs_struct.results( j_res).name, single_result.name)
            % This result is already listed; add its value
            result_exists = true;
            test_runs_struct.results( j_res).cost_function( i_run) = single_result.cost_function;
          end
        end
        if ~result_exists
          % This result is not yet listed; add it
          new_res.name                  = single_result.name;
          new_res.description           = single_result.description;
          new_res.cost_function         = NaN( n_runs,1);
          new_res.cost_function( i_run) = single_result.cost_function;
          if isempty( test_runs_struct.results)
            test_runs_struct.results = new_res;
          else
            test_runs_struct.results( end+1) = new_res;
          end
        end

      end
    end

  end
  function scoreboard = categorise_test_results( single_test_results)
    % Order the test results by category

    scoreboard.subs       = {};
    scoreboard.sub_names  = {};
    scoreboard.has_better = false;
    scoreboard.has_same   = false;
    scoreboard.has_worse  = false;

    for i = 1: length( single_test_results)
      for ri = 1: length( single_test_results( i).results)
        full_name = [single_test_results( i).category '/' single_test_results( i).name '/' single_test_results( i).results( ri).name];
        scoreboard = add_single_test_result( scoreboard, full_name, single_test_results( i).results( ri), ...
          single_test_results( i).git_hash_strings, single_test_results( i).dates_and_times);
      end
    end

  end
  function scoreboard = add_single_test_result( scoreboard, full_name, single_subtest_result, git_hash_strings, dates_and_times)

    if ~isempty( full_name)
      % Move down the category tree

      k = strfind( full_name,'/');

      if isempty( k)
        trunk = full_name;
        branch = [];
      else
        k = k(1);
        trunk = full_name( 1:k-1);
        branch = full_name( k+1:end);
      end
  
      j = find( strcmpi( scoreboard.sub_names, trunk));
      if isempty( j)
        % This category is not yet listed on the scoreboard; add it
        scoreboard.sub_names{ end+1} = trunk;
        j = length( scoreboard.sub_names);
        scoreboard.subs{ j}.subs = [];
        scoreboard.subs{ j}.sub_names = {};
        scoreboard.subs{ j}.has_better = false;
        scoreboard.subs{ j}.has_same   = false;
        scoreboard.subs{ j}.has_worse  = false;
      elseif ~isscalar( j)
        error('whaa!')
      end

      scoreboard.subs{ j} = add_single_test_result( ...
        scoreboard.subs{ j}, branch, single_subtest_result, git_hash_strings, dates_and_times);

      for j = 1: length( scoreboard.subs)
        scoreboard.has_better = scoreboard.has_better || scoreboard.subs{ j}.has_better;
        scoreboard.has_same   = scoreboard.has_same   || scoreboard.subs{ j}.has_same;
        scoreboard.has_worse  = scoreboard.has_worse  || scoreboard.subs{ j}.has_worse;
      end

    else
      % Reached the end of the category tree; add result
      
      scoreboard = rmfield( scoreboard,'subs');
      scoreboard = rmfield( scoreboard,'sub_names');

      scoreboard.description      = single_subtest_result.description;
      scoreboard.git_hash_strings = git_hash_strings;
      scoreboard.dates_and_times  = dates_and_times;
      scoreboard.cost_function    = single_subtest_result.cost_function;

      scoreboard.has_better = false;
      scoreboard.has_same   = false;
      scoreboard.has_worse  = false;

      if isscalar( scoreboard.cost_function)
        scoreboard.has_same = true;
      else
        if scoreboard.cost_function( end) == scoreboard.cost_function( end-1)
          scoreboard.has_same = true;
        elseif scoreboard.cost_function( end) > scoreboard.cost_function( end-1)
          scoreboard.has_worse = true;
        else
          scoreboard.has_better = true;
        end
      end

    end

  end
 
  function process_scoreboard_branch( scoreboard, fid, depth, name)

    p = '';
    for ii = 1: depth
      p = [p ' '];
    end

    if ~isfield( scoreboard,'subs')

      process_scoreboard_leaf( scoreboard, fid, depth);

    else
      for fi = 1: length( scoreboard.subs)

        button_text = scoreboard.sub_names{ fi};
        if scoreboard.subs{ fi}.has_same
          button_text = [button_text ' &#128993'];
        end
        if scoreboard.subs{ fi}.has_better
          button_text = [button_text ' &#128994'];
        end
        if scoreboard.subs{ fi}.has_worse
          button_text = [button_text ' &#128992'];
        end

        fprintf( fid, '%s\n', [p '<div>']);
        fprintf( fid, '%s\n', [p ' <button type="button" class="collapsible" style="font-size:18pt;">' button_text '</button>']);
        fprintf( fid, '%s\n', [p ' <div class="content" style="font-size:18pt;">']);

        process_scoreboard_branch( scoreboard.subs{ fi}, fid, depth+2, scoreboard.sub_names{ fi});

        fprintf( fid, '%s\n', [p ' </div>']);
        fprintf( fid, '%s\n', [p '</div>']);

      end
    end
  end
  function process_scoreboard_leaf( scoreboard, fid, depth, name)

    p = '';
    for ii = 1: depth
      p = [p ' '];
    end

    fprintf(fid,'%s\n',[p '<div>']);
    fprintf(fid,'%s\n',[p ' <p style="font-size:14pt">Description: ' scoreboard.description '</p>']); 

    fprintf(fid,'%s\n',[p ' <div style="height:200px;">']);
    fprintf(fid,'%s\n',[p '  <table>']);
    fprintf(fid,'%s\n',[p '   <tr style="border-bottom: solid 1px #000;">']);
    fprintf(fid,'%s\n',[p '     <th style="width:160px; text-align:left;">Time</th>']);
    fprintf(fid,'%s\n',[p '     <th style="width:160px; text-align:left;">Cost function</th>']);
    fprintf(fid,'%s\n',[p '     <th style="width:500px; text-align:left;">Commit</th>']);
    fprintf(fid,'%s\n',[p '   </tr>']);

    for ii = length( scoreboard.cost_function): -1: 1
      fprintf(fid,'%s\n',[p '   <tr style="border-bottom: solid 1px #000; font-size: 12pt;">']);
      fprintf(fid,'%s\n',[p '    <td>' scoreboard.dates_and_times{ ii} '</td>']);
      fprintf(fid,'%s\n',[p '    <td>' num2str( scoreboard.cost_function( ii), '%14.4e') '</td>']);
      fprintf(fid,'%s\n',[p '    <td>' scoreboard.git_hash_strings{ ii} '</td>']);
      fprintf(fid,'%s\n',[p '   </tr>']);
    end

    fprintf(fid,'%s\n',[p '  </table>']);
    fprintf(fid,'%s\n',[p ' </div>']);
    
    fprintf(fid,'%s\n',[p '</div>']);
      
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