import sys

print('Analysing integrated test idealised/MISMIPplus/MISOMIP ...')
print('')

# Get foldername automated testing

Nargs = len(sys.argv)
assert Nargs < 3, 'Too many arguments'
if Nargs == 1:
  # Assume this is a local run, not implemented
  sys.exit('Need to specify foldername_automated_testing')
else:
  # Assume this is a GitHub Workflow run
 
  foldername_automated_testing = sys.argv[1] 

# Import python script

sys.path.append(f'{foldername_automated_testing}/scoreboard/scripts')
from single_test_run import *

test_name = 'MISOMIP'
test_path = f'integrated_tests/idealised/MISMIPplus'

single_run = Single_run(test_name, test_path) 

single_run.write_scoreboard_file(foldername_automated_testing)
