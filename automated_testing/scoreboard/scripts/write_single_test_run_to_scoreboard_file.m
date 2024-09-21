function write_single_test_run_to_scoreboard_file( single_run, foldername_scoreboard)
% Write a structure containing the results of a single
% component/integrated test run to the corresponding scoreboard file
% (and create that file if it doesn't exist yet)

filename = [foldername_scoreboard '/scoreboard_' single_run.name '.txt'];

% Remove earlier entry for the current commit if it exists
remove_earlier_run_from_this_commit( single_run, filename)

% Open the (new) scoreboard file
fid = fopen( filename,'a');

% Write the test results to the scoreboard
fprintf( fid, '%s\n', '<test_run>');
fprintf( fid, '%s%s\n', '  name           : ', single_run.name);
fprintf( fid, '%s%s\n', '  category       : ', single_run.category);
fprintf( fid, '%s%s\n', '  git hash string: ', single_run.git_hash_string);
fprintf( fid, '%s%s\n', '  date, time     : ', single_run.date_and_time);

for ri = 1: length( single_run.subtests)
  fprintf( fid, '%s\n'      , '  <subtest>');
  fprintf( fid, '%s%s\n'    , '    name         : ', single_run.subtests( ri).name);
  fprintf( fid, '%s%s\n'    , '    description  : ', single_run.subtests( ri).description);
  fprintf( fid, '%s%14.4e\n', '    cost function: ', single_run.subtests( ri).cost_function);
  fprintf( fid, '%s\n'      , '  </subtest>');
end

fprintf( fid, '%s\n', '</test_run>');

% Close the scoreboard file
fclose( fid);

  function remove_earlier_run_from_this_commit( single_run, filename)
    % Remove earlier run from the current commit if it exists

    if ~exist( filename,'file')
      % No scoreboard file exists yet)
      return
    end

    % Open the scoreboard file
    fid = fopen( filename,'r');
    complete_text = textscan( fid, '%s', 'delimiter', '\n', 'whitespace', ''); complete_text = complete_text{1};
    fclose( fid);

    % Find the end of the last run in the file
    i_end = length( complete_text);
    % Safety
    if ~contains( complete_text{ i_end},'</test_run>')
      error('whaa!')
    end

    % Find the start of the last run in the file
    foundit = false;
    i_start = -1;
    for i = i_end-1:-1:1
      if contains( complete_text{i},'<test_run>')
        foundit = true;
        i_start = i;
        break
      end
    end
    % Safety
    if ~foundit
      error('whaa!')
    end

    % Determine the git commit of the last run in the scoreboard file
    foundit = false;
    for i = i_start: i_end
      if contains( complete_text{i},'git hash string')
        foundit = true;
        k = strfind( complete_text{i},': ');
        git_hash_string_file = complete_text{i}(k+2:end);
        break
      end
    end
    % Safety
    if ~foundit
      error('whaa!')
    end

    if strcmp( git_hash_string_file, single_run.git_hash_string)
      % The run in the file is from the current commit; remove it

      if i_start == 1
        % The file only contains this single run; simply remove the entire file
        delete( filename)
        return
      end

      % Remove this particular run from the text
      complete_text( i_start:i_end) = [];

      % Write the reduced text to the file
      fid = fopen( filename,'w');
      for i = 1: length( complete_text)
        fprintf( fid, '%s\n', complete_text{i});
      end
      fclose( fid);

    end

  end

end