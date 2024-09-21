function all_runs_of_test = read_all_runs_of_test_from_scoreboard_file( filename)
% Read the results of all runs of a single component/integrated test
% from its scoreboard file

% Read text from file
fid = fopen( filename);
complete_text = textscan( fid,'%s','delimiter','\n'); complete_text = complete_text{1};
fclose( fid);

all_runs_of_test_as_text = separate_complete_text_into_individual_runs( complete_text);
all_runs_of_test_as_cellstruct = all_runs_of_test_text2cellstruct( all_runs_of_test_as_text);
all_runs_of_test = all_runs_of_test_cellstruct2struct( all_runs_of_test_as_cellstruct);
  
  function test_runs_text = separate_complete_text_into_individual_runs( complete_text)
    % Separate cell array with the complete text of a scoreboard file
    % into separate cell arrays for each individual test run.
  
    % Determine number of test runs
    n_runs = 0;
    for i = 1: length( complete_text)
      single_line = complete_text{ i};
      if contains( single_line,'<test_run>')
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
        if contains( single_line,'</test_run>')
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
  function test_runs_cellstruct = all_runs_of_test_text2cellstruct( test_runs_text)
    % Create a cell array with in each cell a struct containing the results of a single test run
    test_runs_cellstruct = cell(0);
    for i_run = 1: length( test_runs_text)
      test_runs_cellstruct{ end+1} = test_run_cell2struct( test_runs_text{ i_run});
    end
  end
  function test_run_struct = test_run_cell2struct( test_run_text)
    % Create a struct containing the results of a single test run
  
    % Safety
    if ~contains( test_run_text{ 1},'<test_run>') || ...
       ~contains( test_run_text{ end},'</test_run>')
      error('whaa!')
    end
  
    test_run_struct.subtests = [];
  
    is_header = true;
    for i = 1: length( test_run_text)
  
      single_line = test_run_text{ i};
  
      if contains( single_line,'<subtest>')
        is_header = false;
      end
  
      if is_header
        % Header info
        if contains( single_line,'<test_run>')
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
  
      if contains( single_line,'<subtest>')
        % Result info
  
        is_header = false;
  
        found_end = false;
        for i_end = i+1: length( test_run_text)
          single_line = test_run_text{ i_end};
          if contains( single_line,'</subtest>')
            found_end = true;
            break
          end
        end
        % Safety
        if ~found_end
          error('whaa!')
        end
  
        subtest = [];
  
        for j = i+1: i_end-1
          single_line = test_run_text{ j};
          if contains( single_line,'name')
            subtest.name = single_line( strfind( single_line, ': ')+2 : length( single_line));
          elseif contains( single_line,'description')
            subtest.description = single_line( strfind( single_line, ': ')+2 : length( single_line));
          elseif contains( single_line,'cost function')
            subtest.cost_function = str2double( single_line( strfind( single_line, ': ')+2 : length( single_line)));
          else
            error(['unknown subtest line: ' single_line])
          end
        end
  
        if isempty( test_run_struct.subtests)
          test_run_struct.subtests = subtest;
        else
          test_run_struct.subtests( end+1) = subtest;
        end
  
      end
  
    end
  
  end
  function test_runs_struct = all_runs_of_test_cellstruct2struct( test_runs_cellstruct)
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
  
    %% Sub-tests
    test_runs_struct.subtests = [];
  
    for i_run = 1: n_runs
      single_run = test_runs_cellstruct{ i_run};
      for i_sub = 1: length( single_run.subtests)
        single_result = single_run.subtests( i_sub);
  
        % If this subtest is already listed, add its value. If not, add a
        % new subtest
        subtest_exists = false;
        for j_res = 1: length( test_runs_struct.subtests)
          if strcmpi( test_runs_struct.subtests( j_res).name, single_result.name)
            % This subtest is already listed; add its value
            subtest_exists = true;
            test_runs_struct.subtests( j_res).cost_function( i_run) = single_result.cost_function;
          end
        end
        if ~subtest_exists
          % This subtest is not yet listed; add it
          new_subtest.name                  = single_result.name;
          new_subtest.description           = single_result.description;
          new_subtest.cost_function         = NaN( n_runs,1);
          new_subtest.cost_function( i_run) = single_result.cost_function;
          if isempty( test_runs_struct.subtests)
            test_runs_struct.subtests = new_subtest;
          else
            test_runs_struct.subtests( end+1) = new_subtest;
          end
        end
  
      end
    end
  
  end

end