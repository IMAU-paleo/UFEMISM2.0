import os
import glob
import xarray as xr
import matplotlib.pyplot as plt
from main_mesh import Mesh

from utils import *

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
    
    def plot_scalars(self):
        """ Make plot of main scalars """

        #Make outputfolder if necessary
        check_create_dir(f'{self.directory}/figures')

        #Read scalar data
        ds = xr.open_dataset(f'{self.directory}/scalar_output_ANT_00001.nc')

        #Prepare figure
        fig = plt.figure()

        ax = fig.add_subplot(121)
        ax.plot(ds.time,ds.ice_volume,label='ice volume [msle]')
        ax.plot(ds.time,ds.ice_volume_af,label='ice volume af [msle]')
        ax.legend()

        ax = fig.add_subplot(122)
        ax.plot(ds.time,ds.SMB_total,label='SMB [Gt/yr]')
        ax.plot(ds.time,ds.BMB_total,label='BMB [Gt/yr]')
        ax.plot(ds.time,ds.LMB_total,label='LMB [Gt/yr]')
        ax.plot(ds.time,ds.margin_ocean_flux,label='Margin [Gt/yr]')
        ax.plot(ds.time,ds.SMB_total+ds.BMB_total+ds.LMB_total+ds.margin_ocean_flux,c='k',label='Total [Gt/yr]')
        ax.legend()

        #Save figure
        figname = f'{self.directory}/figures/scalars.pdf'
        plt.savefig(figname)
        print(f'Created {figname}')

        #Clean up
        plt.close()
        del fig, ax

        # Close data set
        ds.close()