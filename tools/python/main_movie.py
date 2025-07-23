import os
import glob
import xarray as xr

from colormaps import *
from utils import *
from main_run import *
from main_mesh import *
from main_figure import *

class Movie(object):
    """ Movie of subplots of Ufe/Lad output """

    def __init__(self, Run, variables):

        self.Run = Run
        self.variables = variables

        self.directory = os.path.join(self.Run.directory,'movie')

        self.orientation = 'vertical'
        self.figsize = (7,5)

        self.dxmin = 3e5
        self.dxmax = 1.5e5
        self.dymin = 0
        self.dymax = 0

    def __repr__(self):
        return f"Movie()"

    def __str__(self):
        return f"Movie"

    def make(self, interval=1, length=30):

        ft = 1 # counter for available time frames
        self.fp = 1 # counter for plotted time frame

        if not os.path.isdir(self.directory):
            os.mkdir(f'{self.directory}')

        # Loop over available meshes
        for m in range(1,self.Run.Nmeshes+1):
            # Open mesh
            mesh = Mesh(self.Run, m)

            for t in range(0,mesh.Ntimes):
                if ft in np.arange(1,10001,interval):
                    # Open timeframe
                    tf = Timeframe(mesh, t)

                    # Open figure
                    fig = Figure(self.figsize, directory = self.directory)
                    fig.set_orientation(self.orientation)

                    # Add variable fields
                    for var in self.variables:
                        fig.add_field(tf, var, mask='shelf')

                    # Make and save figure
                    fig.make(f'frame_{self.fp:03d}', dxmin=self.dxmin, dxmax=self.dxmax,
                        dymin=self.dymin, dymax=self.dymax, add_time=True)
                    self.fp += 1
                ft += 1

        # Determine framerate
        framerate = self.fp / length

        # Moviename
        moviename = ''
        for var in self.variables:
            moviename += f'{var}_'
        moviename = moviename[:-1]

        fullmoviename = os.path.join(self.directory,f'{moviename}.mp4')
        # Make movie
        os.system(f'ffmpeg -r {framerate} -f image2 -i {self.directory}/frame_%03d.png -pix_fmt yuv420p -vcodec libx264 -crf 24 {fullmoviename}')

        # Remove frames
        os.system(f'rm {self.directory}/frame*.png')
