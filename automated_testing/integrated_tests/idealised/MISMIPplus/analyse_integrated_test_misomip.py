import sys
import xarray as xr

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

# Basic info of run
test_name = 'MISOMIP'
test_path = f'integrated_tests/idealised/MISMIPplus'

# Initialise run
single_run = Single_run(test_name, test_path) 

# Add cost functions
ds = xr.open_dataset('results_5km_iceocean1r/transect_westeast.nc')
x_GL = ds.grounding_line_distance_from_start.values
ds.close()

err_x_GL_final_lo = abs( min( 0, x_GL[-1] - 430e3));
err_x_GL_final_hi = abs( max( 0, x_GL[-1] - 450e3));

single_run.add_cost_function('err_x_GL_final_lo', 'abs( min( 0, x_GL[-1] - 430e3))', err_x_GL_final_lo)
single_run.add_cost_function('err_x_GL_final_hi', 'abs( max( 0, x_GL[-1] - 450e3))', err_x_GL_final_hi)

# Write scoreboard file
single_run.write_scoreboard_file(foldername_automated_testing)
