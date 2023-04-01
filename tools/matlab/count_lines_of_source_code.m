clc
clear all
close all

main_src_path = '../../src';

henk = dir( main_src_path);

n_tot = 0;

for i = 1: length( henk)
  if henk( i).isdir
    if strcmpi( henk( i).name,'.') || strcmpi( henk( i).name,'..'); continue; end
    f90_files = find_all_lower_f90_files( [main_src_path '/' henk( i).name], []);
    if ~isempty( f90_files)
      n = count_lines( f90_files);
      n_tot = n_tot + n;
      disp(['UFEMISM v2.0 module "' henk( i).name '" contains ' num2str( n) ' lines of code'])
    end
  end
end

disp(['UFEMISM v2.0 in total contains ' num2str( n_tot) ' lines of code'])

function f90_files = find_all_lower_f90_files( henk, f90_files)

piet = dir( henk);

for i = 1: length( piet)
  if strcmpi( piet( i).name,'.') || strcmpi( piet( i).name,'..'); continue; end
  if contains( piet( i).name,'.f90') || contains( piet( i).name,'.F90')
    f90_files{ end+1} = [henk '/' piet( i).name];
  elseif piet( i).isdir
    f90_files = find_all_lower_f90_files( [henk '/' piet( i).name], f90_files);
  end
end

end
function n = count_lines( f90_files)
n = 0;
for i = 1: length( f90_files)
  fid = fopen( f90_files{i});
  temp = textscan( fid,'%s','delimiter','\n');
  fclose( fid);
  n = n + length( temp{1});
end
end