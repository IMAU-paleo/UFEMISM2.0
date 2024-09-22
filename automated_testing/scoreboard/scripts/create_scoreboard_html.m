function create_scoreboard_html( varargin)
% Analyse all the scoreboard files and create the big, interactive
% scoreboard HTML file.

%%

% In the GitHub Workflow, provide the automated_testing folder as
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

scoreboard = read_scoreboard_files( [foldername_automated_testing '/scoreboard/scoreboard_files']);

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
process_scoreboard_branch( scoreboard, fid, 0)
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
      if contains( list_of_files( i).name,'.xml')
        i = i+1;
      else
        list_of_files( i) = [];
      end
    end

    for i = 1: length( list_of_files)
      % For each individual component/integrated tests, the results of all
      % the runs of that test
      disp(['  Reading scoreboard file ' list_of_files( i).name '...'])
      all_tests(i) = read_scoreboard_file( [foldername_scoreboard '/' list_of_files( i).name]);
    end

    scoreboard = categorise_test_results( all_tests);

  end

  function scoreboard = categorise_test_results( all_tests)
    % Order the tests by category

    scoreboard = create_category_tree( all_tests);
    scoreboard = add_results_to_category_tree( scoreboard, all_tests);
    scoreboard = check_results_of_last_run( scoreboard);

  end

  function scoreboard = create_category_tree( all_tests)
    % Create the test category tree

    scoreboard = initialise_empty_category_branch;

    for ti = 1: length( all_tests)
      single_run = all_tests( ti).single_run(1);
      for cfi = 1: length( single_run.cost_functions)
        cost_function = single_run.cost_functions( cfi);
        cost_function_path = [char(single_run.category) '/' char(single_run.name) '/' char(cost_function.name)];
        scoreboard = add_single_cost_function_to_category_tree( scoreboard, cost_function_path, cost_function);
      end
    end

  end
  function empty_branch = initialise_empty_category_branch

    empty_branch.subs       = {};
    empty_branch.sub_names  = {};
    empty_branch.has_better = false;
    empty_branch.has_same   = false;
    empty_branch.has_worse  = false;

  end
  function scoreboard = add_single_cost_function_to_category_tree( scoreboard, cost_function_path, cost_function)
    % Add a single cost function to the test category tree

    if ~contains( cost_function_path, '/')
      % We've reached the end of this cost function's path

      % Safety - check that this cost function is not yet listed
      is_listed = false;
      for ib  = 1: length( scoreboard.subs)
        if strcmpi( scoreboard.sub_names{ ib}, cost_function.name)
          is_listed = true;
          break
        end
      end
      if is_listed
        error('This cost function is already listed!')
      end

      % Add this cost function
      ib = length( scoreboard.subs) + 1;

      scoreboard.sub_names{ ib} = cost_function.name;
      scoreboard.subs{ ib} = initialise_empty_cost_function( cost_function);

    else
      % Move down this cost function's path

      ii = strfind( cost_function_path,'/'); ii = ii(1);
      trunk  = cost_function_path( 1:ii-1);
      branch = cost_function_path( ii+1:end);

      % Check if the trunk is already listed as a branch
      is_listed = false;
      for ib  = 1: length( scoreboard.subs)
        if strcmpi( scoreboard.sub_names{ ib}, trunk)
          is_listed = true;
          break
        end
      end

      if is_listed
        % This category is already listed; move into it.
      else
        % This category is not yet listed; add it.
        ib = length( scoreboard.subs)+1;
        scoreboard.sub_names{ ib} = trunk;
        scoreboard.subs{ ib} = initialise_empty_category_branch;
      end

      scoreboard.subs{ ib} = add_single_cost_function_to_category_tree( ...
        scoreboard.subs{ ib}, branch, cost_function);

    end
  end
  function empty_cost_function = initialise_empty_cost_function( cost_function)

    empty_cost_function.has_better = false;
    empty_cost_function.has_same   = false;
    empty_cost_function.has_worse  = false;

    empty_cost_function.definition       = cost_function.definition;
    empty_cost_function.n_runs           = 0;
    empty_cost_function.git_hash_strings = cell(0);
    empty_cost_function.date_and_times   = cell(0);
    empty_cost_function.values           = [];

  end

  function scoreboard = add_results_to_category_tree( scoreboard, all_tests)
    % Add the individual test results to the ordered category tree

    for ti = 1: length( all_tests)
      single_test = all_tests( ti);
      for ri = 1: length( single_test.single_run)
        single_run = single_test.single_run( ri);
        for cfi = 1: length( single_run.cost_functions)
          cost_function = single_run.cost_functions( cfi);
          cost_function_path = [char(single_run.category) '/' char(single_run.name) '/' char(cost_function.name)];
          scoreboard = add_single_cost_function_value_to_scoreboard( ...
            scoreboard, cost_function_path, cost_function, ...
            single_run.git_hash_string, single_run.date_and_time);
        end
      end
    end

  end
  function scoreboard = add_single_cost_function_value_to_scoreboard( ...
      scoreboard, cost_function_path, cost_function, git_hash_string_of_run, date_and_time_of_run)
    
    if ~contains( cost_function_path, '/')
      % We've reached the end of this cost function's path

      % Add the value to the corresponding branch of the scoreboard
      found_it = false;
      for ib = 1: length( scoreboard.subs)
        if strcmpi( scoreboard.sub_names{ ib}, cost_function_path)
          found_it = true;
          scoreboard.subs{ ib} = add_single_cost_function_value_to_scoreboard_sub( ...
            scoreboard.subs{ ib}, cost_function, git_hash_string_of_run, date_and_time_of_run);
        end
      end
      % Safety
      if ~found_it
        error('Couldnt find corresponding branch on the scoreboard!')
      end

    else
      % Move down this cost function's path

      ii = strfind( cost_function_path,'/'); ii = ii(1);
      trunk  = cost_function_path( 1:ii-1);
      branch = cost_function_path( ii+1:end);

      % Move into the corresponding branch of the scoreboard
      found_it = false;
      for ib = 1: length( scoreboard.subs)
        if strcmpi( scoreboard.sub_names{ ib}, trunk)
          found_it = true;
          scoreboard.subs{ ib} = add_single_cost_function_value_to_scoreboard( ...
            scoreboard.subs{ ib}, branch, cost_function, ...
            git_hash_string_of_run, date_and_time_of_run);
        end
      end
      % Safety
      if ~found_it
        error('Couldnt find corresponding branch on the scoreboard!')
      end

    end

  end
  function scoreboard = add_single_cost_function_value_to_scoreboard_sub( ...
      scoreboard, cost_function, git_hash_string_of_run, date_and_time_of_run)

    scoreboard.n_runs = scoreboard.n_runs + 1;
    n = scoreboard.n_runs;

    scoreboard.git_hash_strings{ n} = git_hash_string_of_run;
    scoreboard.date_and_times{   n} = date_and_time_of_run;
    scoreboard.values(           n) = cost_function.value;

  end

  function scoreboard = check_results_of_last_run( scoreboard)
    % Fill in the has_better, has_worse and has_same flags

    if isfield( scoreboard,'subs')
      % We havent reached the end of the category tree yet

      for ib = 1: length( scoreboard.subs)
        scoreboard.subs{ ib} = check_results_of_last_run( scoreboard.subs{ ib});
        scoreboard.has_better = scoreboard.has_better || scoreboard.subs{ ib}.has_better;
        scoreboard.has_worse  = scoreboard.has_worse  || scoreboard.subs{ ib}.has_worse;
        scoreboard.has_same   = scoreboard.has_same   || scoreboard.subs{ ib}.has_same;
      end

    else
      % We've reached the end of the category tree
      
      if scoreboard.n_runs == 1
        % Trivial answer
        scoreboard.has_same = true;
      else
        % Complicated answer
        if scoreboard.values( end) < scoreboard.values( end-1)
          scoreboard.has_better = true;
        elseif scoreboard.values( end) > scoreboard.values( end-1)
          scoreboard.has_worse = true;
        else
          scoreboard.has_same = true;
        end
      end 
    end
  end

  function process_scoreboard_branch( scoreboard, fid, depth)

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

        process_scoreboard_branch( scoreboard.subs{ fi}, fid, depth+2);

        fprintf( fid, '%s\n', [p ' </div>']);
        fprintf( fid, '%s\n', [p '</div>']);

      end
    end
  end
  function process_scoreboard_leaf( scoreboard, fid, depth)

    p = '';
    for ii = 1: depth
      p = [p ' '];
    end

    fprintf(fid,'%s\n',[p '<div>']);
    fprintf(fid,'%s\n',[p ' <p style="font-size:14pt"><b>Definition: </b>' scoreboard.definition '</p>']); 

    fprintf(fid,'%s\n',[p ' <div style="height:200px;">']);
    fprintf(fid,'%s\n',[p '  <table>']);
    fprintf(fid,'%s\n',[p '   <tr style="border-bottom: solid 1px #000;">']);
    fprintf(fid,'%s\n',[p '     <th style="width:160px; text-align:left;">Time</th>']);
    fprintf(fid,'%s\n',[p '     <th style="width:160px; text-align:left;">Value</th>']);
    fprintf(fid,'%s\n',[p '     <th style="width:500px; text-align:left;">Commit</th>']);
    fprintf(fid,'%s\n',[p '   </tr>']);

    for ii = length( scoreboard.values): -1: 1
      fprintf(fid,'%s\n',[p '   <tr style="border-bottom: solid 1px #000; font-size: 12pt;">']);
      fprintf(fid,'%s\n',[p '    <td>' char(scoreboard.date_and_times{ ii}) '</td>']);
      fprintf(fid,'%s\n',[p '    <td>' num2str( scoreboard.values( ii), '%14.4e') '</td>']);
      fprintf(fid,'%s\n',[p '    <td>' char(scoreboard.git_hash_strings{ ii}) '</td>']);
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