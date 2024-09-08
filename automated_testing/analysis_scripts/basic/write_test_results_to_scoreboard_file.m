function write_test_results_to_scoreboard_file( res, foldername_scoreboard)
% Write a structure containing the results of a single
% component/integrated test to the corresponding scoreboard file
% (and create that file if it doesn't exist yet)

% Remove earlier entry for the current commit if it exists
remove_earlier_entry_for_this_commit( res, foldername_scoreboard)

% Open the (new) scoreboard file
filename_scoreboardfile = [foldername_scoreboard '/scoreboard_' strrep( res.name, '.', 'p') '.txt']
fid_scoreboardfile = fopen( filename_scoreboardfile,'a');

% % DENK DROM
% fid = 1;
% 
% % Write the test results to the scoreboard
% fprintf( fid, '%s\n', '<test>');
% fprintf( fid, '%s%s\n', '  category       : ', res.category);
% fprintf( fid, '%s%s\n', '  git hash string: ', res.git_hash_string);
% fprintf( fid, '%s%s\n', '  date, time     : ', res.date_and_time);
% 
% for ri = 1: length( res.results)
%   fprintf( fid, '%s\n'      , '  <result>');
%   fprintf( fid, '%s%s\n'    , '    name         : ', res.results( ri).name);
%   fprintf( fid, '%s%s\n'    , '    description  : ', res.results( ri).description);
%   fprintf( fid, '%s%14.4e\n', '    cost function: ', res.results( ri).cost_function);
%   fprintf( fid, '%s\n'      , '  </result>');
% end
% 
% fprintf( fid, '%s\n', '</test>');

% Close the scoreboard file
fclose( fid_scoreboardfile);

  function remove_earlier_entry_for_this_commit( res, foldername_scoreboard)
    % Remove earlier entry for the current commit if it exists

    filename = [foldername_scoreboard '/scoreboard_' res.name '.txt'];

    if ~exist( filename,'file')
      % No scoreboard file exists yet)
      return
    end

    % Open the scoreboard file
    fid = fopen( filename,'r');
    complete_text = textscan( fid, '%s', 'delimiter', '\n', 'whitespace', ''); complete_text = complete_text{1};
    fclose( fid);

    % Find the end of the last test in the file
    i_end = length( complete_text);
    % Safety
    if ~contains( complete_text{ i_end},'</test>')
      error('whaa!')
    end

    % Find the start of the last test in the file
    foundit = false;
    i_start = -1;
    for i = i_end-1:-1:1
      if contains( complete_text{i},'<test>')
        foundit = true;
        i_start = i;
        break
      end
    end
    % Safety
    if ~foundit
      error('whaa!')
    end

    % Determine the git commit of the last test in the scoreboard file
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

    if strcmp( git_hash_string_file, res.git_hash_string)
      % The test in the file is from the current commit; remove it

      if i_start == 1
        % The file only contains this single test; simply remove the entire file
        delete( filename)
        return
      end

      % Remove this particular test from the text
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