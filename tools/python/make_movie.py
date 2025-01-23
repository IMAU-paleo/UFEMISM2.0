import os
import glob
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from colormaps import *

from single_figure import *


def make_movie(vars,inputfolder,outputfolder,run,ncols=1,dpi=1200):
    #Get available meshes
    meshmax = int(sorted(glob.glob(f'{inputfolder}/{run}/main_output_ANT_0*.nc'))[-1][-8:-3])
    
    #Loop over available meshes
    for mesh in range(1,meshmax+1):
        #Get available times
        fname = f'{inputfolder}/{run}/main_output_ANT_{mesh:05d}.nc'
        dss = xr.open_dataset(fname)
        lent = len(dss.time)
    
        #Loop over available time slices:
        for t in range(0,lent):
            plot_figure(vars,inputfolder,outputfolder,run,mesh,t,ncols=ncols,dpi=dpi)
        dss.close()
        print(f'== Finished mesh {mesh} of {meshmax} ==')

    #Make video
    moviename = 'movie'
    framerate = 10

    command = f'ffmpeg -r {framerate} -f image2 -i {outputfolder}/frame_020%03d.png -pix_fmt yuv420p -vcodec libx264 -crf 24 {outputfolder}/{moviename}.mp4'
    os.system(command)

    return