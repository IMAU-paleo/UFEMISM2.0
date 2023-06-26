clc
clear all
close all

main_src_path = '../../src';

henk = dir( main_src_path);

n_tot = 0;

R.n     = [];
R.names = {};

for i = 1: length( henk)
  if henk( i).isdir
    if strcmpi( henk( i).name,'.') || strcmpi( henk( i).name,'..'); continue; end
    f90_files = find_all_f90_files( [main_src_path '/' henk( i).name]);
    if ~isempty( f90_files)
      n = count_lines( f90_files);
      n_tot = n_tot + n;
      R.n( end+1,1) = n;
      R.names{ end+1} = strrep( strrep( strrep( [henk( i).name ': ' num2str( n)], ...
        'surface_mass_balance','SMB'),...
        'basal_mass_balance','BMB'),...
        '_','\_');
    end
  end
end

% Sort by number of lines
[n_sorted,ind] = sortrows( R.n);

ind2 = [];
i = 0;
while ~isempty( ind)
  i = 1-i;
  if i==1
    ind2( end+1) = ind( 1);
    ind( 1) = [];
  else
    ind2( end+1) = ind( end);
    ind( end) = [];
  end
end
ind = ind2;

R.n     = R.n( ind);
names_new = {};
for i = 1: length( R.names)
  names_new{ end+1} = R.names{ ind( i)};
end
R.names = names_new;

disp(['UFEMISM v2.0 in total contains ' num2str( n_tot) ' lines of code'])

%% plot
p = pie( R.n, ones( size( R.n)), R.names);
for i = 1: length( p)
  if strcmpi( class( p( i)), 'matlab.graphics.primitive.Text')
    set( p( i), 'fontsize', 24)
  end
end
set( gcf,'position',[744   346   827   704],'color','w');
set( gca,'position',[0.15, 0.1, 0.55, 0.8],'fontsize',24,'ylim',[-1.2,1.4]);
title( gca,['Total: ' num2str( n_tot)])

function f90_files = find_all_f90_files( henk)

f90_files = {};

piet = dir( henk);

for i = 1: length( piet)
  if strcmpi( piet( i).name,'.') || strcmpi( piet( i).name,'..'); continue; end
  if contains( piet( i).name,'.f90') || contains( piet( i).name,'.F90')
    f90_files{ end+1} = [henk '/' piet( i).name];
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