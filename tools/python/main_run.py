import os
import glob
import xarray as xr
import matplotlib.pyplot as plt
from main_mesh import Mesh

from utils import *
from analyse_resources import *

class Run(object):
    """ Properties and functions from a Ufemism run """

    def __init__(self, directory):
        """ Gather basic info from run """

        self.directory  = directory
        self.runname = os.path.basename(self.directory)

        try:
           self.Nmeshes = int(sorted(glob.glob(f'{self.directory}/main_output_ANT_0*.nc'))[-1][-8:-3])
        except:
           self.Nmeshes = 0

    def __repr__(self):
        return f"Run('{self.directory}')"

    def __str__(self):
        """ Spit out info on this run """
        return f"Run '{self.runname}' containing {self.Nmeshes} meshes"
    
    def get_mesh(self,mesh_number,file='main_output_ANT'):
        """ Gather info on a given mesh """

        assert mesh_number <= self.Nmeshes, 'Mesh number too high, not available in output'
        assert mesh_number >= 1, 'Mesh number too low, should be at least 1'

        mesh = Mesh(self,mesh_number,file=file)

        return mesh
    
    def get_laddie_mesh(self,mesh_number,file='laddie_output_fields_ANT'):
        """ Gather info on a given mesh """

        mesh = Mesh(self,mesh_number,file=file)

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

    def make_movie(self,variables,framerate=10):
        """ Make a movie of variables """

        f = 1 #Frame counter

        #Loop over available meshes
        for m in range(1,self.Nmeshes+1):
            mesh = self.get_mesh(m)

            #Loop over time slices
            for t in range(0,mesh.Ntimes):
                mesh.make_plot(variables,t,f,purpose='movie')
                f += 1
        
        #Make video
        moviename = f'{self.directory}/movie/'
        for varname in variables:
                moviename += f'{varname}_'
        moviename = moviename[:-1] #Remove last underscore

        os.system(f'ffmpeg -r {framerate} -f image2 -i {self.directory}/movie/frame_%03d.png -pix_fmt yuv420p -vcodec libx264 -crf 24 {moviename}.mp4')

        #Remove frames
        os.system(f'rm {self.directory}/movie/frame*.png')

    def plot_comptime(self,Nlevels=9,toplevel=0,toproutine='UFEMISM_program'):
        """ Plot computation time of subroutines """

        #Make outputfolder if necessary
        check_create_dir(f'{self.directory}/figures')

        #Read scalar data
        ds = xr.open_dataset(f'{self.directory}/resource_tracking.nc')
        ds = ds.sum(dim='time')

        #Get names and computation times of subroutines down to required level
        names,tcomps = get_allnames(ds,maxlev=toplevel+Nlevels+2)

        #Output name
        savename = f'{self.directory}/figures/resources_{toproutine}_{Nlevels:.0f}.pdf'

        #Make plot
        plot_resources(names,tcomps,Nlevels,savename=savename,toplevel=toplevel,toproutine=toproutine)

        # Close data set
        ds.close()
        
