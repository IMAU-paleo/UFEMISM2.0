clc
clear all
close all

filename_in  = '../../src/basic/model_configuration.f90';
filename_out = '../../config-files/unit_tests.cfg';

%% Read model_configuration.f90 source code

fid = fopen( filename_in,'r');
config = textscan( fid,'%s','delimiter','\n','whitespace','');
config = config{1};
fclose( fid);

%% Find start and end of the parameter definition block

parameters_i1 = 0;
parameters_i2 = 0;

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
  
end

%% Copy the parameters block
param_block = {};
for i = parameters_i1: parameters_i2
  param_block{ end+1} = config{ i};
end

%% Convert to config file format

config_file = {};

for i = 1: length( param_block)
  
  % Read line from the parameters block
  str = param_block{ i};
  
  % Check if a parameter is defined in this line
  ii = strfind( str,'_config');
  if ~isempty( ii)
    % Remove the type declaration
    jj = strfind( str,'::');
    str = str( jj+3:end);
    % Append four leading spaces
    str = ['    ' str];
  end
 
  % Skip the lines with the values for the SELEN moving time window
  if contains( str,'(/20._dp, 20._dp, 20._dp, 5._dp, 5._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, &') || ...
     contains( str,'0._dp,  0._dp,  0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, &') || ...
     contains( str,'0._dp,  0._dp,  0._dp, 0._dp, 0._dp  /)')
    continue
  end
  
  % Replace '._dp' by '.0  '
  str = strrep( str, '._dp', '.0  ');
  
  % Replace '_dp' by '   '
  str = strrep( str, '_dp', '   ');
  
  % Fix array declaration brackets in SELEN stuff
  if contains( str,'SELEN_irreg_time_window_config')
    str = ['    SELEN_irreg_time_window_config          = 20.0  , 20.0  , 20.0  , 5.0  , 5.0  , 1.0  , 1.0  , '...
      '1.0  , 1.0  , 1.0  , 1.0  , 1.0  , 1.0  , 1.0  , 1.0  , 0.0  ,  0.0  ,  0.0  , 0.0  , 0.0  , 0.0  , 0.0  , '...
      '0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,  0.0  ,  0.0  ,  0.0  , 0.0  , 0.0  , 0.0  , 0.0  , '...
      '0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,  0.0  ,  0.0  , 0.0  , 0.0'];
  end
  if contains( str,'SELEN_visc_prof_config')
    str = ['    SELEN_visc_prof_config                       = 3.0  , 0.6   , 0.3              '...
      '! Viscosities of viscous asthenosphere layers [?]'];
  end
  
  % Make sure this config actually performs the unit tests
  if contains( str,'create_procedural_output_dir_config')
    str = ['    create_procedural_output_dir_config          = .TRUE.                           '...
      '! Automatically create an output directory with a procedural name (e.g. results_20210720_001/)'];
  end
  if contains( str,'fixed_output_dir_config')
    str = ['    fixed_output_dir_config                      = ''results_unit_tests''             '...
      '! If not, create a directory with this name instead (stops the program if this directory already exists)'];
  end
  if contains( str,'do_unit_tests_config')
    str = ['    do_unit_tests_config                         = .TRUE.                           '...
      '! Whether or not to (only) perform the unit tests in the main_validation module'];
  end
  
  % Add to the type block
  config_file{ end+1} = str;
  
end

%% Write to new file
fid = fopen( filename_out,'w');
fprintf( fid,'%s\n\n','&CONFIG');
for i = 1: length( config_file)
  str = config_file{ i};
  fprintf( fid, '%s\n', str);
end
fprintf( fid,'%s\n','');
fprintf( fid,'%s\n','/');
fclose( fid);