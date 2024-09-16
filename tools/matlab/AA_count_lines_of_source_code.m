function AA_count_lines_of_source_code

n_tot = 0;

R.n     = [];
R.names = {};

%% src

main_path = '../../src';
henk = dir( main_path);

for i = 1: length( henk)
  if henk( i).isdir
    if strcmpi( henk( i).name,'.') || strcmpi( henk( i).name,'..'); continue; end
    code_files = find_all_code_files( [main_path '/' henk( i).name], '.f90');
    if ~isempty( code_files)
      n = count_lines( code_files);
      n_tot = n_tot + n;
      R.n( end+1,1) = n;
      R.names{ end+1} = henk( i).name;
    end
  end
end

%% automated_testing

main_path = '../../automated_testing';
henk = dir( main_path);

for i = 1: length( henk)
  if henk( i).isdir
    if strcmpi( henk( i).name,'.') || strcmpi( henk( i).name,'..'); continue; end
    code_files = find_all_code_files( [main_path '/' henk( i).name], '.m');
    if ~isempty( code_files)
      n = count_lines( code_files);
      n_tot = n_tot + n;
      R.n( end+1,1) = n;
      R.names{ end+1} = henk( i).name;
    end
  end
end

%% plot

disp(['UFEMISM v2.0 contains ' num2str( n_tot) ' lines of code'])

% Sort by number of lines
[~,ind] = sortrows( R.n);

wa = 800;
ha = 600;

H.Fig = figure('position',[200,200,wa,ha],'color','w');
H.Ax = axes('parent',H.Fig,'position',[0,0,1,1],'xtick',[],'ytick',[],'xlim',[0,1],'ylim',[-0.1,1.2]);

text( 0.15,1.1, 'Total:'       ,'fontsize',24,'fontweight','bold')
text( 0.25,1.1, num2str( n_tot),'fontsize',24,'fontweight','bold')

% draw patches
n_rel_sum = 0;
colors = lines( length( R.n));
for i = 1: length( R.n)
  
  ii = ind( i);
  
  n_rel = R.n( ii) / n_tot;
  
  % Patch
  xmin = 0.05;
  xmax = 0.15;
  ymin = n_rel_sum;
  ymax = n_rel_sum + n_rel;
  n_rel_sum = n_rel_sum + n_rel;
  
  patch('xdata',[xmin,xmax,xmax,xmin],'ydata',[ymin,ymin,ymax,ymax],'facecolor',colors( i,:));
  
  % Text
  x = 0.24;
  y = (i-1)/(length( R.n)-1);
  line( 'xdata',[xmax+0.01,x],'ydata',[(ymin+ymax)/2,y],'color',colors( i,:),'linewidth',2)
  
  x = 0.25;
  str = num2str( R.n( ii));
  text( x,y,str,'fontsize',24,'color',colors( i,:));
  
  x = 0.35;
  f = round( 1000 * R.n( ii) / n_tot) / 10;
  str = ['(' num2str( f) '%)'];
  text( x,y,str,'fontsize',24,'color',colors( i,:));
  
  x = 0.50;
  str = strrep( R.names{ ii},'_','\_');
  text( x,y,str,'fontsize',24,'color',colors( i,:));
  
end

function code_files = find_all_code_files( jan, extension)

code_files = {};
piet = dir( jan);

for p = 1: length( piet)
  if strcmpi( piet( p).name,'.') || strcmpi( piet( p).name,'..')
    continue
  end

  if piet( p).isdir
    code_files_more = find_all_code_files( [jan '/' piet( p).name], extension);
    for pp = 1: length( code_files_more)
      code_files{ end+1} = code_files_more{ pp};
    end
  end
  
  if contains( piet( p).name, extension, 'Ignorecase',true)
    code_files{ end+1} = [jan '/' piet( p).name];
  end
end

end
function n = count_lines( f90_files)
n = 0;
for fi = 1: length( f90_files)
  fid = fopen( f90_files{ fi});
  temp = textscan( fid,'%s','delimiter','\n');
  fclose( fid);
  n = n + length( temp{1});
end
end

end