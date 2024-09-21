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
      % For each individual component/integrated tests, the results of all
      % the runs of that test
      disp(['  Reading scoreboard file ' list_of_files( i).name '...'])
      all_runs_of_test(i) = read_all_runs_of_test_from_scoreboard_file( [foldername_scoreboard '/' list_of_files( i).name]);
    end

    scoreboard = categorise_test_results( all_runs_of_test);

  end

  function scoreboard = categorise_test_results( all_runs_of_test)
    % Order the test results by category

    scoreboard.subs       = {};
    scoreboard.sub_names  = {};
    scoreboard.has_better = false;
    scoreboard.has_same   = false;
    scoreboard.has_worse  = false;

    for i = 1: length( all_runs_of_test)
      for ri = 1: length( all_runs_of_test( i).subtests)
        full_name = [all_runs_of_test( i).category '/' all_runs_of_test( i).name '/' all_runs_of_test( i).subtests( ri).name];
        scoreboard = add_single_test_result( scoreboard, full_name, all_runs_of_test( i).subtests( ri), ...
          all_runs_of_test( i).git_hash_strings, all_runs_of_test( i).dates_and_times);
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