function timeframes = get_UFEMISM_filelist( foldername, region)
  % Return a list of all the timeframes found in this simulatio,
  % and the files they are contained in

%     time                 = timeframes( tfi).time;
%     ti                   = timeframes( tfi).ti;
%     filename_restart     = timeframes( tfi).filename_restart;
%     filename_help_fields = timeframes( tfi).filename_help_fields;

filename_restart     = [foldername '/restart_'     region '_00001.nc'];
filename_help_fields = [foldername '/help_fields_' region '_00001.nc'];

if ~exist( filename_restart,'file')
  error(['Couldnt find any output in folder ' foldername '!'])
end

% Initialise list
timeframes = [];

nf = 1;
while exist( filename_restart,'file')
  
  % Find all timeframes contained in this file
  time = ncread( filename_restart,'time');
  ti   = 1:length(time);
  
  % Add them to the list
  for tfi = 1:length(time)
    timeframes(end+1).time                 = time( tfi);
    timeframes(end  ).ti                   = tfi;
    timeframes(end  ).filename_restart     = filename_restart;
    timeframes(end  ).filename_help_fields = filename_help_fields;
  end
  
  % Move to the next file
  nf = nf + 1;
  filename_restart     = [foldername '/restart_'     region '_' n2str( nf) '.nc'];
  filename_help_fields = [foldername '/help_fields_' region '_' n2str( nf) '.nc'];
  
end

  function n_str = n2str( n)
    if (n < 10)
      n_str = ['0000' num2str(n)];
    elseif (n < 100)
      n_str = ['000' num2str(n)];
    elseif (n < 1000)
      n_str = ['00' num2str(n)];
    elseif (n < 10000)
      n_str = ['0' num2str(n)];
    else
      n_str = num2str(n);
    end
  end
    
end