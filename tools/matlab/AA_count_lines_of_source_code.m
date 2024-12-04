function AA_count_lines_of_source_code

n_lines_UFE = 0;

R.names      = {};
R.n_modules  = [];
R.n_lines    = [];
R.n_comments = [];
R.r_comments = [];

%% src

main_path = '../../src';
henk = dir( main_path);

for i = 1: length( henk)
  if henk( i).isdir
    if strcmpi( henk( i).name,'.') || strcmpi( henk( i).name,'..'); continue; end
    code_files = find_all_code_files( [main_path '/' henk( i).name], '.f90');
    if ~isempty( code_files)
      R.names{ end+1} = henk( i).name;

      [n_modules, n_lines, n_comments] = count_lines( code_files);

      n_lines_UFE = n_lines_UFE + n_lines;

      R.n_modules(  end+1,1) = n_modules;
      R.n_lines(    end+1,1) = n_lines;
      R.n_comments( end+1,1) = n_comments;
      R.r_comments( end+1,1) = n_comments / (n_lines - n_comments);
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
      R.names{ end+1} = henk( i).name;

      [n_modules, n_lines, n_comments] = count_lines( code_files);

      n_lines_UFE = n_lines_UFE + n_lines;

      R.n_modules(  end+1,1) = n_modules;
      R.n_lines(    end+1,1) = n_lines;
      R.n_comments( end+1,1) = n_comments;
      R.r_comments( end+1,1) = n_comments / (n_lines - n_comments);
    end
  end
end

%% plot

disp(['UFEMISM v2.0 contains ' num2str( n_lines_UFE) ' lines of code, spread over ' ...
  num2str( sum( R.n_modules)) ' module files.'])

% Sort by number of lines
[~,ind] = sortrows( R.n_lines);

wa = 1000;
ha = 600;

H.Fig = figure('position',[50,200,wa,ha],'color','w');
H.Ax = axes('parent',H.Fig,'position',[0,0,1,1],'xtick',[],'ytick',[],'xlim',[0,1.25],'ylim',[-0.05,1.2]);

text( 0.05, 1.15, ['Total number of lines: ' num2str( n_lines_UFE)],'fontsize',24,'fontweight','bold')
text( 0.25, 1.07, 'Component'  ,'fontsize',20,'fontweight','bold')
text( 0.65, 1.07, 'Lines'      ,'fontsize',20,'fontweight','bold')
text( 0.75, 1.07, '(fraction)' ,'fontsize',20,'fontweight','bold')
text( 0.87, 1.07, 'Mods'       ,'fontsize',20,'fontweight','bold')
text( 0.96, 1.07, 'Lin/mod'    ,'fontsize',20,'fontweight','bold')
text( 1.07, 1.07, 'Comts/lin'  ,'fontsize',20,'fontweight','bold')

% draw patches
n_rel_sum = 0;
colors = lines( length( R.n_lines));
for i = 1: length( R.n_lines)
  
  ii = ind( i);
  
  n_rel = R.n_lines( ii) / n_lines_UFE;
  
  % Patch
  xmin = 0.05;
  xmax = 0.15;
  ymin = n_rel_sum;
  ymax = n_rel_sum + n_rel;
  n_rel_sum = n_rel_sum + n_rel;
  patch('xdata',[xmin,xmax,xmax,xmin],'ydata',[ymin,ymin,ymax,ymax],'facecolor',colors( i,:));
  
  % Lines from patches to text
  x = 0.24;
  y = (i-1)/(length( R.n_lines)-1);
  line( 'xdata',[xmax+0.01,x],'ydata',[(ymin+ymax)/2,y],'color',colors( i,:),'linewidth',2)
  
  % Component name
  x = 0.25;
  str = strrep( R.names{ ii},'_','\_');
  text( x,y,str,'fontsize',24,'color',colors( i,:));
  
  % Lines (total)
  x = 0.65;
  str = num2str( R.n_lines( ii));
  text( x,y,str,'fontsize',24,'color',colors( i,:));
  
  % Lines (fraction)
  x = 0.75;
  f = round( 1000 * R.n_lines( ii) / n_lines_UFE) / 10;
  str = ['(' num2str( f) '%)'];
  text( x,y,str,'fontsize',24,'color',colors( i,:));
  
  % Number of modules
  x = 0.87;
  str = num2str( R.n_modules( ii));
  text( x,y,str,'fontsize',24,'color',colors( i,:));
  
  % Lines per module
  x = 0.96;
  str = num2str( floor( R.n_lines( ii) / R.n_modules( ii)));
  text( x,y,str,'fontsize',24,'color',colors( i,:));
  
  % Comment lines per functional line
  x = 1.07;
  str = num2str( round( 1000 * R.r_comments( ii)) / 1000);
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

  function [n_modules, n_lines, n_comments] = count_lines( f90_files)

    n_modules  = length( f90_files);
    n_lines    = 0;
    n_comments = 0;
    
    for fi = 1: length( f90_files)
      fid = fopen( f90_files{ fi});
      temp = textscan( fid,'%s','delimiter','\n');
      fclose( fid);
  
      % Total number of lines in this single f90 file
      n_lines_file =  length( temp{1});
      n_lines = n_lines + n_lines_file;
    
      % Count number of comment lines in this single f90 file
      n_comments_file = 0;
      for i = 1: length( temp{1})
        single_line = strjust( temp{1}{i},'left');
        if isempty(single_line); continue; end
        if strcmpi(single_line(1),'!') || strcmpi(single_line(1),'%')
          n_comments_file = n_comments_file + 1;
        end
      end
      n_comments = n_comments + n_comments_file;
  
    end
  
  end

end