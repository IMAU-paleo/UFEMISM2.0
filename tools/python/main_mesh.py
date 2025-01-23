import os
import glob
import xarray as xr
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

from colormaps import *
from utils import *

class Mesh(object):
    """ Properties and functions of a single mesh """

    def __init__(self, directory, mesh):
        """ Gather basic info from run """

        self.directory = directory
        self.mesh = mesh

        self.got_voronois = False
        self.got_triangles = False

        self.open()
        self.close()

    def make_plot(self,variables,t,f=None,ncols=1,dpi=1200,format='png',purpose='single'):
        """ Make a figure of requested variables at time slice t """

        assert purpose in ['single','movie'], 'ERROR: invalid purpose in make_plot'

        if purpose == 'movie':
            assert f is not None, "ERROR: provide frame number f in make_plot or change purpose to 'single' "
            assert format == 'png', "ERROR: format for movie must be 'png'"

        # Open data set
        self.open()

        #Prepare figure
        fig,ax = plt.subplots(len(variables),ncols,sharex=True,sharey=True)

        #Add data (voronois / triangles)
        for v,varname in enumerate(variables):
            pcoll = self.get_pcoll(varname,t)

            ax[v].add_collection(pcoll)

            #Make up subplot
            ax[v].set_xlim([self.ds.xmin,self.ds.xmax])
            ax[v].set_ylim([self.ds.ymin,self.ds.ymax])
            ax[v].set_aspect(1)

        #Add time if movie
        if purpose=='movie':
            ax[0].set_title(f'Year {self.ds.time[t].values:.0f}')

        #Save figure
        if purpose=='single':
            check_create_dir(f'{self.directory}/figures')
            figname = f'{self.directory}/figures/'
            for varname in variables:
                figname += f'{varname}_'
            figname += f'{int(self.ds.time[t].values):06d}.{format}'
        elif purpose=='movie':
            check_create_dir(f'{self.directory}/movie')
            figname = f'{self.directory}/movie/frame_{f:03d}.png'

        plt.savefig(figname,dpi=dpi)
        print(f'Created {figname}')

        #Clean up
        plt.close()
        del fig, ax, pcoll

        # Close data set
        self.close()

    def open(self):
        self.ds = xr.open_dataset(f'{self.directory}/main_output_ANT_{self.mesh:05d}.nc')
        self.Ntimes = len(self.ds.time)
    
    def close(self):
        self.ds.close()

    def get_voronois(self):
        """ Extract Voronoi cells as patches """
        self.voronois = []
        for vi in range(0,len(self.ds.vi)):
            nVVor = self.ds.nVVor[vi].values
            VVor = self.ds.VVor[:nVVor,vi].values
            Vor = self.ds.Vor[:,VVor-1].values
            self.voronois.append(Polygon(Vor.T))
        
        self.got_voronois = True

    def get_triangles(self):
        """ Extract triangles as patches """
        self.triangles = []
        #for ti in range(0,len(self.ds.ti)):
            #nVVor = self.ds.nVVor[vi].values
            #VVor = self.ds.VVor[:nVVor,vi].values
            #Vor = self.ds.Vor[:,VVor-1].values
            #self.triangles.append(Polygon(Tri.T))
        
        self.got_triangles = True

    def get_pcoll(self,varname,t):
        """ Get patch collection """

        #Get data
        var = self.get_data(varname,t)

        #Get colormap info
        cmap,norm = get_cmap(varname)

        #Check type (voronoi / triangle)
        if 'vi' in var.dims:
            if not self.got_voronois:
                self.get_voronois()
            pcoll = PatchCollection(self.voronois,cmap=cmap,norm=norm)
        elif 'ti' in var.dims:
            if not self.got_triangles:
                self.get_triangles()
            pcoll = PatchCollection(self.triangles,cmap=cmap,norm=norm)
        else:
            print(f'ERROR: variable {varname} is not on vertices or triangles')

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