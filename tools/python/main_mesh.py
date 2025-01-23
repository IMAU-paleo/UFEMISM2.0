import os
import glob
import xarray as xr
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

from main_run import Run
from colormaps import *

class Mesh(Run):
    """ Properties and functions of a single mesh """

    def __init__(self, directory, mesh):
        """ Gather basic info from run """
        super().__init__(directory)

        self.mesh = mesh
        
        assert self.mesh <= self.Nmeshes, 'Mesh number too high, not available in output'

    def open(self):
        self.ds = xr.open_dataset(f'{self.directory}/main_output_ANT_{self.mesh:05d}.nc')
        self.Ntimes = len(self.ds.time)
    
    def close(self):
        self.ds.close()

    def get_voronoi(self):
        """ Extract Voronoi cells as patches """
        self.voronois = []
        for vi in range(0,len(self.ds.vi)):
            nVVor = self.ds.nVVor[vi].values
            VVor = self.ds.VVor[:nVVor,vi].values
            Vor = self.ds.Vor[:,VVor-1].values
            self.voronois.append(Polygon(Vor.T))

    def get_triangles(self):
        """ Extract triangles as patches """
        self.triangles = []
        #for ti in range(0,len(self.ds.ti)):
            #nVVor = self.ds.nVVor[vi].values
            #VVor = self.ds.VVor[:nVVor,vi].values
            #Vor = self.ds.Vor[:,VVor-1].values
            #self.triangles.append(Polygon(Tri.T))

    def get_pcoll(self,varname,t):
        """ Get patch collection """

        #Get data
        var = self.get_data(varname,t)

        #Get colormap info
        cmap,norm = get_cmap(varname)

        #Check type (voronoi / triangle)
        if 'vi' in var.dims:
            pcoll = PatchCollection(self.voronois,cmap=cmap,norm=norm)
        elif 'ti' in var.dims:
            pcoll = PatchCollection(self.triangles,cmap=cmap,norm=norm)
        else:
            print(f'ERROR: variable {var} is not on vertices or triangles')

        # Fill array
        if varname == 'BMB':
            #Reverse values for BMB
            pcoll.set_array(-var.values)
        else:
            pcoll.set_array(var.values)

        return pcoll

    def get_data(self,varname,t):
        """ Get data array of variable for time slice t"""

        # Read variable 
        try: 
            var = self.ds[varname]
        except:
            print(f'ERROR: could not read variable {varname}, make sure dataset is open and variable exists')
            return

        # Extract time slice
        try:
            var = var.isel(time=t)
        except:
            print(f'ERROR: requested time slice {t} is not available')
            return
    
        return var