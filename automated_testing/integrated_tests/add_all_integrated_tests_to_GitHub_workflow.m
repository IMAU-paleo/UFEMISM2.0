function add_all_integrated_tests_to_GitHub_workflow
% Add all the integrated tests to the GitHub Workflow file so they will be
% run automatically
%
% NOTE: meant to be run manually by the user after adding a new integrated
% test
%
% NOTE: this script must be run from inside automated_testing/integrated_tests!

foldername_workflows = '../.github/workflows';
foldername_integrated_tests = 'integrated_tests';

list_of_tests = list_all_integrated_tests( foldername_integrated_tests);

create_run_integrated_tests_workflow( list_of_tests)

for i = 1: length( list_of_tests)
  create_run_single_test_workflow( list_of_tests{ i});
end

create_analyse_integrated_tests_workflow( list_of_tests)

  function create_run_integrated_tests_workflow( list_of_tests)

    filename_run_integrated_tests_workflow  = [foldername_workflows ...
      '/UFE_test_suite_run_integrated_tests.yml'];

    fid = fopen( filename_run_integrated_tests_workflow,'w');

    fprintf( fid, '%s\n', 'name: UFEMISM Test Suite - run integrated tests');
    fprintf( fid, '%s\n', 'run-name: ${{ github.actor }} - UFEMISM Test Suite - run integrated tests');
    fprintf( fid, '%s\n', 'on:');
    fprintf( fid, '%s\n', '  workflow_call:');
    fprintf( fid, '%s\n', '');
    fprintf( fid, '%s\n', 'jobs:');

    for i = 1: length( list_of_tests)

      test_path = list_of_tests{ i};
      test_name = strrep( test_path,'/','_');

      fprintf( fid, '%s\n', '');
      fprintf( fid, '%s\n', ['  ' test_name ':']);
      fprintf( fid, '%s\n', ['    uses: ./.github/workflows/zz_' test_name '.yml']);

    end

    fclose( fid);

  end

  function create_run_single_test_workflow( test_path)

    ii = strfind( test_path,'/'); ii = ii( end);
    test_path_firstname = test_path( ii+1:end);

    test_name           = strrep( test_path          ,'/','_');
    test_name_firstname = strrep( test_path_firstname,'/','_');

    % Read dummy workflow file
    filename_workflow_dummy = 'integrated_test_single_workflow_dummy.txt';
    fid = fopen( filename_workflow_dummy,'r');
    temp = textscan( fid,'%s','delimiter','\n','whitespace',''); temp = temp{1};
    fclose( fid);

    % Place test name in file
    for ii = 1: length( temp)
      temp{ii} = strrep( temp{ii}, '!!test_path!!'          , test_path);
      temp{ii} = strrep( temp{ii}, '!!test_name!!'          , test_name);
      temp{ii} = strrep( temp{ii}, '!!test_name_firstname!!', test_name_firstname);
    end

    % Write to single test workflow file
    filename_workflow = [foldername_workflows '/zz_' test_name '.yml'];
    fid = fopen( filename_workflow,'w');
    for ii = 1: length( temp)
      fprintf( fid,'%s\n',temp{ii});
    end
    fclose( fid);

  end

  function create_analyse_integrated_tests_workflow( list_of_tests)

    filename_analyse_integrated_tests_workflow  = [foldername_workflows ...
      '/UFE_test_suite_analyse_integrated_tests.yml'];

    fid = fopen( filename_analyse_integrated_tests_workflow,'w');

    fprintf( fid, '%s\n', 'name: UFEMISM Test Suite - analyse integrated tests');
    fprintf( fid, '%s\n', 'run-name: ${{ github.actor }} - UFEMISM Test Suite - analyse integrated tests');
    fprintf( fid, '%s\n', 'on:');
    fprintf( fid, '%s\n', '  workflow_call:');
    fprintf( fid, '%s\n', '');
    fprintf( fid, '%s\n', 'jobs:');
    fprintf( fid, '%s\n', '  analyse_integrated_tests:');
    fprintf( fid, '%s\n', '    runs-on: macos-latest');
    fprintf( fid, '%s\n', '    steps:');
    fprintf( fid, '%s\n', '');
    fprintf( fid, '%s\n', '      - name: Checkout UFEMISM repository (from pull request)');
    fprintf( fid, '%s\n', "        if: ${{ github.event_name == 'pull_request' }}");
    fprintf( fid, '%s\n', '        uses: actions/checkout@v4');
    fprintf( fid, '%s\n', '        with:');
    fprintf( fid, '%s\n', '         repository: ${{ github.event.pull_request.head.repo.full_name }}');
    fprintf( fid, '%s\n', '         ref: ${{ github.event.pull_request.head.ref }}');
    fprintf( fid, '%s\n', '');
    fprintf( fid, '%s\n', '      - name: Checkout UFEMISM repository (from manual run)');
    fprintf( fid, '%s\n', "        if: ${{ github.event_name != 'pull_request' }}");
    fprintf( fid, '%s\n', '        uses: actions/checkout@v4');
    fprintf( fid, '%s\n', '');
    fprintf( fid, '%s\n', '      - name: Install MATLAB');
    fprintf( fid, '%s\n', '        uses: matlab-actions/setup-matlab@v2.2.0');
    fprintf( fid, '%s\n', '        with:');
    fprintf( fid, '%s\n', '          cache: true');

    for i = 1: length( list_of_tests)

      test_path = list_of_tests{ i};
      test_name = strrep( test_path,'/','_');

      fprintf( fid, '%s\n', '');
      fprintf( fid, '%s\n', ['# ' test_path]);
      fprintf( fid, '%s\n', '');
      fprintf( fid, '%s\n',['      - name: Download temporary scoreboard file for ' test_path]);
      fprintf( fid, '%s\n', '        uses: actions/download-artifact@v4');
      fprintf( fid, '%s\n', '        with:');
      fprintf( fid, '%s\n',['          name: temporary_scoreboard_' test_name]);
      fprintf( fid, '%s\n',['          path: automated_testing/' test_path]);
    end

    fclose( fid);

  end

end