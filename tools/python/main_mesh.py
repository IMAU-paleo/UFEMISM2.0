import os
import glob
import xarray as xr
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

from colormaps import *
from utils import *

class Mesh(object):
    """ Properties and functions of a single mesh """

    def __init__(self, directory, mesh_number, file='main_output_ANT'):
        """ Gather basic info from run """

        self.directory = directory
        self.mesh_number = mesh_number
        self.file = file

        self.got_voronois = False
        self.got_triangles = False

        self.open()
        self.close()

    def __repr__(self):
        return f"Mesh('{self.directory}',{self.mesh_number},'{self.file}')"

    def __str__(self):
        return f"Mesh number {self.mesh_number} of Run '{self.directory}'"

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

            #Add data
            pcoll = self.get_pcoll(varname,t)
            im = ax[v].add_collection(pcoll)
            im.cmap.set_bad('magenta')

            #Add grounding line
            gl = self.get_GL(t)
            ax[v].plot(gl[0,:],gl[1,:],c='yellow',lw=.25)

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
        self.ds = xr.open_dataset(f'{self.directory}/{self.file}_{self.mesh_number:05d}.nc')
        self.Ntimes = len(self.ds.time)
    
    def close(self):
        self.ds.close()

    def get_voronois(self):
        """ Extract Voronoi cells as patches """
        self.voronois = []

        print(f"Computing {len(self.ds.vi)} voronoi polygons ...")

        for vi in range(0,len(self.ds.vi)):
            nVVor = self.ds.nVVor[vi].values
            VVor = self.ds.VVor[:nVVor,vi].values
            Vor = self.ds.Vor[:,VVor-1].values
            self.voronois.append(Polygon(Vor.T))
        
        self.got_voronois = True

        return f"Finished computing voronoi polygons"

    def get_triangles(self):
        """ Extract triangles as patches """
        self.triangles = []

        print(f"Computing {len(self.ds.ti)} triangle polygons ...")

        for ti in range(0,len(self.ds.ti)):
            Tri = self.ds.Tri[:,ti].values
            V = self.ds.V[:,Tri-1].values
            self.triangles.append(Polygon(V.T))
        
        self.got_triangles = True

        return f"Finished computing triangle polygons"

class Timeframe(object):
    """ Single timeframe of a mesh """

    def __init__(self, Mesh, t):
        """ Gather basic info from run """

        self.Mesh = Mesh
        self.t = t

        self.ds = self.Mesh.ds.isel(time=t)

        self.got_gl = False
        self.got_mask = False

    def __repr__(self):
        return f"Timeframe({repr(self.Mesh)},{self.t})"

    def __str__(self):
        return f"Timeframe {self.t} of Mesh number {self.mesh_number} of Run '{self.directory}'"

    def get_gl(self):
        """ Extract grounding line """

        # Read variable

        try:
            var = self.ds['grounding_line'].values
        except KeyError:
            print(f"ERROR: 'grounding_line' not in output files")
            return

        self.gl = var
        self.got_gl = True

        return

    def get_mask(self):
        """ Get a reduced mask separating grounded ice, ice shelf and ocean """

        mask = self.ds['mask']
        mask = xr.where(mask==1,3,mask)
        mask = xr.where(mask==5,3,mask)
        mask = xr.where(mask==7,3,mask)
        mask = xr.where(mask==9,3,mask)
        mask = xr.where(mask==10,3,mask)
        mask = xr.where(mask==6,4,mask)
        mask = xr.where(mask==8,2,mask)

        self.mask = mask.values
        self.got_mask = True

        return

    def get_pcoll(self,varname):
        """ Get patch collection """

        #Get data
        if varname == 'Uabs_lad':
            var1 = self.get_data('U_lad')
            var2 = self.get_data('V_lad')
            var = (var1**2+var2**2)**.5
        elif varname[:3] == 'BMB':
            var = self.get_data('BMB')
        else:
            var = self.get_data(varname)

        #Get colormap info
        cmap,norm = get_cmap(varname)

        #Check type (voronoi / triangle)
        if 'vi' in var.dims:
            if not self.Mesh.got_voronois:
                self.Mesh.get_voronois()
            pcoll = PatchCollection(self.Mesh.voronois,cmap=cmap,norm=norm)
        elif 'ti' in var.dims:
            if not self.Mesh.got_triangles:
                self.Mesh.get_triangles()
            pcoll = PatchCollection(self.Mesh.triangles,cmap=cmap,norm=norm)
        else:
            print(f'ERROR: variable {varname} is not on vertices or triangles')

        # Fill array
        if varname[:3] == 'BMB':
            #Reverse values for BMB
            pcoll.set_array(-var.values)
        else:
            pcoll.set_array(var.values)

        return pcoll

    def get_data(self,varname):
        """ Get data array of variable """

        # Read variable
        try: 
            var = self.ds[varname]
        except:
            print(f'ERROR: could not read variable {varname}, make sure dataset is open and variable exists')
            return

        return var
