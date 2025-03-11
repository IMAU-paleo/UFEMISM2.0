import numpy as np
import datetime as dt
import git
import xml.etree.cElementTree as ET

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

  def add_cost_function(self, name, definition, value):

    # Get index of new cost function
    Navail = len(self.cost_functions)
    if (Navail == 1 and self.cost_functions[0].name == ''):
      nnext = 0
    else:
      nnext = Navail

    # Add new cost function
    self.cost_functions[nnext] = Cost_function(name, definition, value)

  def get_filename(self, foldername_automated_testing):

    # Some replacements to shorten filenames
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

    # Create filename
    replcategory = self.category
    for i,repl in enumerate(string_replacements_test_category):
      replcategory = replcategory.replace(repl[0],repl[1])

    self.filename = f"{foldername_automated_testing}/scoreboard/temporary_scoreboard_files/{replcategory}_{self.name}_{self.git_hash_string}.xml"

  def write_scoreboard_file(self, foldername_automated_testing):

    # Get filename of scoreboard file
    self.get_filename(foldername_automated_testing)

    # Set up header
    root = ET.Element('single_run')
    ET.SubElement(root, "name").text = self.name
    ET.SubElement(root, "category").text = self.category
    ET.SubElement(root, "date_and_time").text = self.date_and_time
    ET.SubElement(root, "git_hash_string").text = self.git_hash_string
    
    # Add cost functions
    for c in range(len(self.cost_functions)):
      costfunction = self.cost_functions[c]
      cf = ET.SubElement(root,'cost_functions')
    
      ET.SubElement(cf, "name").text = costfunction.name
      ET.SubElement(cf, "definition").text = costfunction.definition
      if isinstance(costfunction.value, int):
        ET.SubElement(cf, "value").text = f"{costfunction.value:.0f}"
      elif isinstance(costfunction.value, float):
        ET.SubElement(cf, "value").text = f"{costfunction.value:.4f}"
      else:
        sys.exit('Invalid value type in cost function, must be int or float')
    
    # Write xml file
    tree = ET.ElementTree(root)
    ET.indent(tree, '    ')
    tree.write(self.filename,encoding='UTF-8',xml_declaration=True)

    # Print success
    print(f'Successfully wrote scoreboard file {self.filename}')

class Cost_function:
  def __init__(self, name, definition, value):
    self.name       = name
    self.definition = definition
    self.value      = value

    assert isinstance(self.name, str), "Cost function name should be string"
    assert isinstance(self.definition, str), "Cost function definition should be string"
