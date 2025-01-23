import os
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from colormaps import *

def plot_figure(vars,inputfolder,outputfolder,run,mesh,t,ncols=1,dpi=1200):
    #Read data
    fname = f'{inputfolder}/{run}/main_output_ANT_{mesh:05d}.nc'
    ds = xr.open_dataset(fname)

    #Extract Voronoi cells
    patches = []
    for vi in range(0,len(ds.vi)):
        nVVor = ds.nVVor[vi].values
        VVor = ds.VVor[:nVVor,vi].values
        Vor = ds.Vor[:,VVor-1].values
        patches.append(Polygon(Vor.T))

    #Get colormap data
    cmap,norm = get_colormaps()

    #Prepare figure
    fig,ax = plt.subplots(len(vars),ncols)

    #Loop over requested variables
    for v,var in enumerate(vars):

        #Create patch with variable-dependent colormap data
        p = PatchCollection(patches,cmap=cmap[var],norm=norm[var])

        #Make subplot
        dax = ax[v]

        # Get data
        if var == 'BMB':
            #Reverse values for BMB
            cols = [-ds[var][t,vi].values for vi in range(0,len(ds.vi))]
        else:
            cols = [ds[var][t,vi].values for vi in range(0,len(ds.vi))]
        
        #Attach data to patch
        p.set_array(cols)

        #Plot cells
        dax.add_collection(p)

        #Make up grid of subplot
        dax.set_xlim([ds.xmin,ds.xmax])
        dax.set_ylim([ds.ymin,ds.ymax])
        dax.set_aspect(1)

    #Make outputfolder if necessary
    if not os.path.isdir(outputfolder):
        os.makedirs(outputfolder)

    #Save figure
    plt.savefig(f'{outputfolder}/frame_{int(ds.time[t].values):06d}',dpi=dpi)
    plt.close()
    print(f'Created frame_{int(ds.time[t].values):06d}')
