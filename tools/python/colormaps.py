import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import cmocean as cmo
from copy import copy

def get_colormaps():
    cmap = {}
    norm = {}

    #Create BMB colormap
    vmax = 100
    vmin = -10
    linthresh = .3
    linscale = .25
    fracpos = (np.log10(vmax/linthresh)+linscale)/(np.log10(vmax/linthresh)+np.log10(-(vmin/linthresh))+2*linscale)
    nneg = np.int_((1-fracpos)*256)
    colors1 = plt.get_cmap('cmo.dense_r')(np.linspace(0,1.,nneg+1))
    colors2 = plt.get_cmap('gist_heat_r')(np.linspace(0., 1, 256-nneg-1))
    colors = np.vstack((colors1, colors2))

    cmap['BMB'] = mpl.colors.LinearSegmentedColormap.from_list('my_colormap', colors)
    norm['BMB'] = mpl.colors.SymLogNorm(linthresh, vmin=vmin, vmax=vmax, linscale=linscale)

    cmap['Hi'] = copy(plt.get_cmap('cmo.ice'))
    norm['Hi'] = mpl.colors.Normalize(vmin=0,vmax=4000,clip=True)

    cmap['Hs'] = copy(plt.get_cmap('cmo.ice'))
    norm['Hs'] = mpl.colors.Normalize(vmin=0,vmax=1000,clip=True)

    cmap['uabs_surf'] = copy(plt.get_cmap('CMRmap_r'))
    norm['uabs_surf'] = mpl.colors.LogNorm(vmin=1.,vmax=3000,clip=True)

    return cmap,norm