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

create_workflow_file_run_and_analyse_integrated_tests( list_of_tests)

for i = 1: length( list_of_tests)
  create_workflow_file_run_single_test( list_of_tests{ i});
end

create_workflow_file_finalise_scoreboard( list_of_tests)

  function create_workflow_file_run_and_analyse_integrated_tests( list_of_tests)

    filename_run_and_analyse_integrated_tests_workflow  = [foldername_workflows ...
      '/UFE_test_suite_run_and_analyse_integrated_tests.yml'];

    fid = fopen( filename_run_and_analyse_integrated_tests_workflow,'w');

    fprintf( fid, '%s\n', '# NOTE: this script is created automatically by running');
    fprintf( fid, '%s\n', '# ''automated_testing/integrated_tests/add_all_integrated_tests_to_Github_workflow.m''');
    fprintf( fid, '%s\n', '');
    fprintf( fid, '%s\n', 'name: UFEMISM Test Suite - run and analyse integrated tests');
    fprintf( fid, '%s\n', 'run-name: ${{ github.actor }} - UFEMISM Test Suite - run and analyse integrated tests');
    fprintf( fid, '%s\n', 'on:');
    fprintf( fid, '%s\n', '  workflow_call:');
    fprintf( fid, '%s\n', '');
    fprintf( fid, '%s\n', 'jobs:');

    for ti = 1: length( list_of_tests)

      test_path = list_of_tests{ ti};
      test_name = strrep( test_path,'/','_');

      fprintf( fid, '%s\n', '');
      fprintf( fid, '%s\n', ['  ' test_name ':']);
      fprintf( fid, '%s\n', ['    uses: ./.github/workflows/zz_' test_name '.yml']);

    end

    fclose( fid);

  end

  function create_workflow_file_run_single_test( test_path)

    ii = strfind( test_path,'/'); ii = ii( end);
    test_path_firstname = test_path( ii+1:end);

    test_name           = strrep( test_path          ,'/','_');
    test_name_firstname = strrep( test_path_firstname,'/','_');

    % Read dummy workflow file
    if contains( test_name,'idealised')
      filename_workflow_dummy = 'integrated_test_idealised_single_workflow_dummy.txt';
    elseif contains( test_name,'realistic')
      filename_workflow_dummy = 'integrated_test_realistic_single_workflow_dummy.txt';
    else
      error('dont know if this test is idealised or realistic!')
    end
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

  function create_workflow_file_finalise_scoreboard( list_of_tests)

    filename_workflow_append_scoreboard_files  = [foldername_workflows ...
      '/UFE_test_suite_finalise_scoreboard.yml'];

    fid = fopen( filename_workflow_append_scoreboard_files,'w');

    fprintf( fid, '%s\n', '# NOTE: this script is created automatically by running');
    fprintf( fid, '%s\n', '# ''automated_testing/integrated_tests/add_all_integrated_tests_to_Github_workflow.m''');
    fprintf( fid, '%s\n', '');
    fprintf( fid, '%s\n', 'name: UFEMISM Test Suite - finalise scoreboard');
    fprintf( fid, '%s\n', 'run-name: ${{ github.actor }} - UFEMISM Test Suite - finalise scoreboard');
    fprintf( fid, '%s\n', 'on:');
    fprintf( fid, '%s\n', '  workflow_call:');
    fprintf( fid, '%s\n', '');
    fprintf( fid, '%s\n', 'jobs:');
    fprintf( fid, '%s\n', '  finalise_scoreboard:');
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
    fprintf( fid, '%s\n', '');
    fprintf( fid, '%s\n', '# =========================================================');
    fprintf( fid, '%s\n', '# ===== Download scoreboard files for component tests =====');
    fprintf( fid, '%s\n', '# =========================================================');
    fprintf( fid, '%s\n', '');
    fprintf( fid, '%s\n', '      - name: Download scoreboard files for component tests');
    fprintf( fid, '%s\n', '        uses: actions/download-artifact@v4');
    fprintf( fid, '%s\n', '        with:');
    fprintf( fid, '%s\n', '          name: scoreboard_files_component_tests');
    fprintf( fid, '%s\n', '          path: automated_testing/scoreboard/scoreboard_files');
    fprintf( fid, '%s\n', '');
    fprintf( fid, '%s\n', '# ==========================================================');
    fprintf( fid, '%s\n', '# ===== Download scoreboard files for integrated tests =====');
    fprintf( fid, '%s\n', '# ==========================================================');
    fprintf( fid, '%s\n', '#');
    fprintf( fid, '%s\n', '# NOTE: list created automatically; if you want to add new integrated tests,');
    fprintf( fid, '%s\n', '# just run ''automated_testing/integrated_tests/add_all_integrated_tests_to_Github_workflow.m'' again');

    for ti = 1: length( list_of_tests)

      test_path = list_of_tests{ ti};
      test_name = strrep( test_path,'/','_');

      fprintf( fid, '%s\n', '');
      fprintf( fid, '%s\n',['      - name: Download scoreboard file for ' test_path]);
      fprintf( fid, '%s\n', '        uses: actions/download-artifact@v4');
      fprintf( fid, '%s\n', '        with:');
      fprintf( fid, '%s\n',['          name: scoreboard_file_' test_name]);
      fprintf( fid, '%s\n', '          path: automated_testing/scoreboard/scoreboard_files');
    end

    fprintf( fid, '%s\n', '');
    fprintf( fid, '%s\n', '# ===============================');
    fprintf( fid, '%s\n', '# ===== Finalise scoreboard =====');
    fprintf( fid, '%s\n', '# ===============================');
    fprintf( fid, '%s\n', '');
    fprintf( fid, '%s\n', '      - name: Commit scoreboard files');
    fprintf( fid, '%s\n', '        # See https://github.com/marketplace/actions/add-commit');
    fprintf( fid, '%s\n', '        if: ${{ github.event_name == ''pull_request'' }} # Only do this for pull requests');
    fprintf( fid, '%s\n', '        uses: EndBug/add-and-commit@v9 # You can change this to use a specific version.');
    fprintf( fid, '%s\n', '        with:');
    fprintf( fid, '%s\n', '          add: automated_testing/scoreboard/scoreboard_files/*.xml');
    fprintf( fid, '%s\n', '          author_name: ${{ github.actor }} (from UFEMISM test suite workflow)');
    fprintf( fid, '%s\n', '          message: ''Update scoreboard files (from UFEMISM test suite workflow by ${{ github.actor }})''');
    fprintf( fid, '%s\n', '          push: true');
    fprintf( fid, '%s\n', '');
    fprintf( fid, '%s\n', '      - name: Create scoreboard visualisation');
    fprintf( fid, '%s\n', '        uses: matlab-actions/run-command@v2');
    fprintf( fid, '%s\n', '        with:');
    fprintf( fid, '%s\n', '          command: |');
    fprintf( fid, '%s\n', '            addpath(''automated_testing/scoreboard/scripts'')');
    fprintf( fid, '%s\n', '            create_scoreboard_html(''${{github.workspace}}/automated_testing'')');
    fprintf( fid, '%s\n', '');
    fprintf( fid, '%s\n', '      - name: Upload scoreboard visualisation as artifact');
    fprintf( fid, '%s\n', '        uses: actions/upload-artifact@v4.3.4');
    fprintf( fid, '%s\n', '        with:');
    fprintf( fid, '%s\n', '          name: scoreboard');
    fprintf( fid, '%s\n', '          path: automated_testing/test_reports/scoreboard.html');

    fclose( fid);

  end

end