import os
import glob
from main_mesh import Mesh

class Run(object):
    """ Properties and functions from a Ufemism run """

    def __init__(self, directory):
        """ Gather basic info from run """

        self.directory  = directory
        self.runname = os.path.basename(self.directory)

        self.Nmeshes = int(sorted(glob.glob(f'{self.directory}/main_output_ANT_0*.nc'))[-1][-8:-3])
    
    def get_mesh(self,mesh):
        """ Gather info on a given mesh """

        assert mesh <= self.Nmeshes, 'Mesh number too high, not available in output'

        mesh = Mesh(self.directory,mesh)
        
        return mesh