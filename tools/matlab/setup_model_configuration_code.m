clc
clear all
close all

filename_in  = '../../src/basic/model_configuration.f90';
filename_out = filename_in;

%% Read model_configuration.f90 source code

fid = fopen( filename_in,'r');
config = textscan( fid,'%s','delimiter','\n','whitespace','');
config = config{1};
fclose( fid);

%% Find start and end of the four blocks

parameters_i1 = 0;
parameters_i2 = 0;

type_i1       = 0;
type_i2       = 0;

namelist_i1   = 0;
namelist_i2   = 0;

valuecopy_i1  = 0;
valuecopy_i2  = 0;

for i = 1: length( config)-4
  
  % Section where the config parameters are defined
  if     ~isempty( strfind( config{ i  },'! ===== Configuration variables =====')) && ...
         ~isempty( strfind( config{ i+1},'! ===================================')) && ...
         ~isempty( strfind( config{ i+3},'  ! The "_config" variables, which will be collected into a NAMELIST, and replaced')) && ...
         ~isempty( strfind( config{ i+4},'  ! by the values in the external config file. Remember the "_config" extension!'))
    parameters_i1 = i+6;
  elseif ~isempty( strfind( config{ i  },'! ===== Configuration variables - end =====')) && ...
         ~isempty( strfind( config{ i+1},'! ========================================='))
    parameters_i2 = i-2;
  end
  
  % Section where the config type is defined
  if     ~isempty( strfind( config{ i  },'  TYPE type_config')) && ...
         ~isempty( strfind( config{ i+1},'    ! The different parameters that control a UFEMISM simulation'))
    type_i1 = i+3;
  elseif ~isempty( strfind( config{ i  },'  ! == Non-configurable variables')) && ...
         ~isempty( strfind( config{ i+1},'  ! ============================='))
    type_i2 = i-2;
  end
  
  % Section where the namelist is defined
  if     ~isempty( strfind( config{ i  },'    NAMELIST /CONFIG/&'))
    namelist_i1 = i+1;
  elseif ~isempty( strfind( config{ i  },'    ! End of the config NAMELIST'))
    namelist_i2 = i-1;
  end
  
  % Section where the config variables are copied to the config structure
  if     ~isempty( strfind( config{ i  },'    ! Copy the values of the _config variables to the C structure'))
    valuecopy_i1 = i+2;
  elseif ~isempty( strfind( config{ i  },'    ! Finished copying the values of the _config variables to the C structure'))
    valuecopy_i2 = i-2;
  end
  
end

%% Copy the parameters block
param_block = {};
for i = parameters_i1: parameters_i2
  param_block{ end+1} = config{ i};
end

%% Create the type block

% For lines that define parameters, remove the "_config" extension and 
% everything after that.

type_block = {};

for i = 1: length( param_block)
  
  % Read line from the parameters block
  str = param_block{ i};
  
  % Check if a parameter is defined in this line
  ii = strfind( str,'_config');
  if ~isempty( ii)
    % Remove the "_config" extension and everything after that.
    str = str( 1:ii-1);
  end
 
  % Skip the lines with the values for the SELEN moving time window
  if contains( str,'(/20._dp, 20._dp, 20._dp, 5._dp, 5._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, &') || ...
     contains( str,'0._dp,  0._dp,  0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, &') || ...
     contains( str,'0._dp,  0._dp,  0._dp, 0._dp, 0._dp  /)')
    continue
  end
  
  % Add to the type block
  type_block{ end+1} = str;
  
end

%% Create the namelist block

% Only keep the lines where the parameters are defined

namelist_block = {};

for i = 1: length( param_block)
  
  % Read line from the parameters block
  str = param_block{ i};
  
  % Check if a parameter is defined in this line
  ii = strfind( str,'_config');
  if ~isempty( ii)
    
    % Remove the value assignment and everything after that.
    jj = strfind( str,'::');
    str = str( jj+3:end);
    jj = strfind( str,'=');
    str = str( 1:jj-1);
    
    % Append whitespaces until all strings are of equal length
    while length( str) < 60
      str = [str ' '];
    end
    
    % Append '      ' at the start
    str = ['      ' str];
    
    % Append ', &' at the end
    str = [str ', &'];
 
    % Skip the lines with the values for the SELEN moving time window
    if contains( str,'(/20._dp, 20._dp, 20._dp, 5._dp, 5._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, &') || ...
       contains( str,'0._dp,  0._dp,  0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, &') || ...
       contains( str,'0._dp,  0._dp,  0._dp, 0._dp, 0._dp  /)')
      continue
    end
  
    % Add to the type block
    namelist_block{ end+1} = str;
    
  end
  
end

% Remove ', &' from the last entry
str = namelist_block{ end};
str = str( 1:end-3);
namelist_block{ end} = str;

%% Create the copy-values-to-structure block

valuecopy_block = {};

for i = 1: length( param_block)
  
  % Read line from the parameters block
  str = param_block{ i};
  
  % Check if a parameter is defined in this line
  ii = strfind( str,'_config');
  if ~isempty( ii)
    
    % Distill the parameter name
    jj = strfind( str,'::');
    str = str( jj+3:end);
    jj = strfind( str,'_config');
    str = str( 1:jj-1);
    param_name = str;

    % Append '    C%' to the start
    str = ['    C%' param_name];
    
    % Append whitespaces until all strings are of equal length
    while length( str) < 60
      str = [str ' '];
    end
    
    % Append value assignment
    str = [str ' = ' param_name '_config'];
    
  end
 
  % Skip the lines with the values for the SELEN moving time window
  if contains( str,'(/20._dp, 20._dp, 20._dp, 5._dp, 5._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, &') || ...
     contains( str,'0._dp,  0._dp,  0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, &') || ...
     contains( str,'0._dp,  0._dp,  0._dp, 0._dp, 0._dp  /)')
    continue
  end
  
  % Add to the copy-values-to-structure block
  valuecopy_block{ end+1} = str;
  
end

%% Substitute new blocks into source code

config_new = {};

% Part before parameters block
for i = 1: parameters_i1-1
  config_new{ end+1} = config{ i};
end

% Parameters block
for i = 1: length( param_block)
  config_new{ end+1} = param_block{ i};
end

% Part between parameters block and type block
for i = parameters_i2+1 : type_i1-1
  config_new{ end+1} = config{ i};
end

% Type block
for i = 1: length( type_block)
  config_new{ end+1} = type_block{ i};
end

% Part between type block and namelist block
for i = type_i2+1 : namelist_i1-1
  config_new{ end+1} = config{ i};
end

% Namelist block
for i = 1: length( namelist_block)
  config_new{ end+1} = namelist_block{ i};
end

% Part between namelist block and valuecopy block
for i = namelist_i2+1 : valuecopy_i1-1
  config_new{ end+1} = config{ i};
end

% Valuecopy block
for i = 1: length( valuecopy_block)
  config_new{ end+1} = valuecopy_block{ i};
end

% Part after valuecopy block
for i = valuecopy_i2+1 : length( config)
  config_new{ end+1} = config{ i};
end

% Remove trailing whitespaces
for i = 1: length( config_new)
  str = config_new{ i};
  if ~isempty( str)
    while strcmpi( str( end),' ')
      str = str( 1:end-1);
    end
  end
  config_new{ i} = str;
end

%% Write to new file
fid = fopen( filename_out,'w');
for i = 1: length( config_new)
  str = config_new{ i};
  fprintf( fid, '%s\n', str);
end
fclose( fid);