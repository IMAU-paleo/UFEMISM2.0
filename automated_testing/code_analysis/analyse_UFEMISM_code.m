function analyse_UFEMISM_code
% Perform some basic analysis of the UFEMISM source code, e.g.:
% - number of lines per subroutine for each module,
% - number of comment lines per functional line for each module,
% - whether or not a modile has code privacy statements

src_dir = '../../src';

code_tree = initialise_empty_code_tree_branch;
code_tree = add_branch_to_code_tree( code_tree, src_dir);

  function code_tree = add_branch_to_code_tree( code_tree, src_dir)
    % Recursively add (sub-)directories or code files to the code tree

    code_tree.path = src_dir;

    list_of_files = dir( src_dir);

    for i = 1: length( list_of_files)
      if strcmpi( list_of_files(i).name,'.') || ...
          strcmpi( list_of_files(i).name,'..')
        continue
      end

      if list_of_files(i).isdir
        code_tree.sub_names{ end+1} = list_of_files(i).name;
        code_tree.subs{ end+1} = initialise_empty_code_tree_branch;
        code_tree.subs{ end} = add_branch_to_code_tree( code_tree.subs{ end}, ...
          [src_dir '/' list_of_files(i).name]);
      elseif strcmpi(list_of_files(i).name(end-3:end),'.f90')
        code_tree.sub_names{ end+1} = list_of_files(i).name;
        code_tree.subs{ end+1} = initialise_empty_code_tree_branch;
        code_tree.subs{ end} = add_leaf_to_code_tree( code_tree.subs{ end}, ...
          [src_dir '/' list_of_files(i).name]);
      end
    end

    % Calculate integrated values
    for si = 1: length( code_tree.subs)
      code_tree.n_lines = code_tree.n_lines + ...
        code_tree.subs{ si}.n_lines;
      code_tree.n_functional_lines = code_tree.n_functional_lines + ...
        code_tree.subs{ si}.n_functional_lines;
      code_tree.n_comment_lines = code_tree.n_comment_lines + ...
        code_tree.subs{ si}.n_comment_lines;
      code_tree.n_subroutines = code_tree.n_subroutines + ...
        code_tree.subs{ si}.n_subroutines;
      code_tree.n_modules = code_tree.n_modules + ...
        code_tree.subs{ si}.n_modules;
      code_tree.has_privacy_statement = code_tree.has_privacy_statement + ...
        code_tree.subs{ si}.has_privacy_statement;
    end

  end
  function code_tree = add_leaf_to_code_tree( code_tree, filename)

    code_tree.path = filename;

    code_tree = analyse_source_code_file( code_tree, filename);

  end

  function branch = initialise_empty_code_tree_branch
    branch.path      = '';
    branch.subs      = {};
    branch.sub_names = {};

    branch.n_lines            = 0;
    branch.n_functional_lines = 0;
    branch.n_comment_lines    = 0;
    branch.n_subroutines      = 0;
    branch.n_modules          = 0;
    branch.has_privacy_statement            = 0;
  end

  function code_tree = analyse_source_code_file( code_tree, filename)

    % Open, read, and close source code file
    fid = fopen( filename,'r');
    temp = textscan( fid,'%s','delimiter','\n'); lines = temp{1};
    fclose( fid);

    code_tree.n_lines               = n_lines( lines);
    code_tree.n_functional_lines    = n_functional_lines( lines);
    code_tree.n_comment_lines       = n_comment_lines( lines);
    code_tree.n_subroutines         = n_subroutines( lines);
    code_tree.n_modules             = n_modules( lines);
    code_tree.has_privacy_statement = has_privacy_statement( lines);

  end
  function n = n_lines( lines)
    n = length( lines);
  end
  function n = n_functional_lines( lines)
    n = 0;
    for li = 1: length( lines)
      line = strjust( lines{li}, 'left');
      if ~isempty( line) && ~startsWith( line,{'!','#'})
        n = n+1;
      end
    end
  end
  function n = n_comment_lines( lines)
    n = 0;
    for li = 1: length( lines)
      line = strjust( lines{li}, 'left');
      if startsWith( line,{'!','#'})
        n = n+1;
      end
    end
  end
  function n = n_subroutines( lines)
    n = 0;
    for li = 1: length( lines)
      line = strjust( lines{li}, 'left');
      if startsWith( line, 'subroutine', 'IgnoreCase',true)
        n = n+1;
      end
    end
  end
  function n = n_modules( lines)
    n = 0;
    for li = 1: length( lines)
      line = strjust( lines{li}, 'left');
      if startsWith( line, 'module', 'IgnoreCase',true)
        n = n+1;
      end
    end
  end
  function n = has_privacy_statement( lines)
    n = 0;
    for li = 1: length( lines)
      line = strjust( lines{li}, 'left');
      if startsWith( line, {'private','public'}, 'IgnoreCase',true)
        n = 1;
      end
    end
  end

end