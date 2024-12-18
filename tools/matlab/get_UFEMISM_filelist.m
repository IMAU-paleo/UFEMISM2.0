function timeframes = get_UFEMISM_filelist( foldername, region)
% Return a list of all the timeframes found in this simulation,
% and the files they are contained in

filename = [foldername '/main_output_' region '_00001.nc'];

if ~exist( filename,'file')
  error(['Couldnt find any output in folder ' foldername '!'])
end

% Initialise list
timeframes = [];

nf = 1;
while exist( filename,'file')
  
  % Find all timeframes contained in this file
  time = ncread( filename,'time');
  ti   = 1:length(time);
  
  % Add them to the list
  for tfi = 1:length(time)
    timeframes(end+1).time     = time( tfi);
    timeframes(end  ).ti       = tfi;
    timeframes(end  ).filename = filename;
  end
  
  % Move to the next file
  nf = nf + 1;
  filename = [foldername '/main_output_' region '_' n2str( nf) '.nc'];
  
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