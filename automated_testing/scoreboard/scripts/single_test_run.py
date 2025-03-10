import numpy as np
import datetime as dt
import git

class Single_run:
  def __init__(self, test_name, test_category):
    self.name            = test_name
    self.category        = test_category
    self.date_and_time   = dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    repo = git.Repo(search_parent_directories=True)
    self.git_hash_string = repo.head.object.hexsha

    assert isinstance(self.name, str), "Single run test_name should be string"
    assert isinstance(self.category, str), "Single run test_category should be string"

    self.cost_functions = {}

    cost_function_dummy = Cost_function('','',0.)
    self.cost_functions[0] = cost_function_dummy 

    print(self.name)
    print(self.category)
    print(self.date_and_time)
    print(self.git_hash_string)

  def write_scoreboard_file(self, foldername_automated_testing):

    string_replacements_test_category = [ 
      ['/',                        '_'],
      ['component_tests',          'ct'],
      ['integrated_tests',         'it'],
      ['discretisation',           'disc'],
      ['mapping_and_derivatives',  'map_deriv'],
      ['remapping',                'remap'],
      ['mesh_to_grid',             'm2g'],
      ['grid_to_mesh',             'g2m'],
      ['mesh_to_mesh',             'm2m'],
      ['idealised',                'ideal'],
      ['Halfar',                   'Hlf']]

    replcategory = self.category
    for i,repl in enumerate(string_replacements_test_category):
      replcategory = replcategory.replace(repl[0],repl[1])
    print(self.category,replcategory)

    filename = f"{foldername_automated_testing}/scoreboard/temporary_scoreboard_files/{replcategory}_{self.name}_{self.git_hash_string}.xml"

    print(filename)

class Cost_function:
  def __init__(self, name, definition, value):
    self.name       = name
    self.definition = definition
    self.value      = value

    assert isinstance(self.name, str), "Cost function name should be string"
    assert isinstance(self.definition, str), "Cost function definition should be string"
    assert isinstance(self.value, float), "Cost function value should be float"
