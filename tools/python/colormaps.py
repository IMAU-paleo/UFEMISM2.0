import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import cmocean as cmo
from copy import copy

def get_cmap(varname):

    if varname == 'BMB':
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

        cmap = mpl.colors.LinearSegmentedColormap.from_list('my_colormap', colors)
        norm = mpl.colors.SymLogNorm(linthresh, vmin=vmin, vmax=vmax, linscale=linscale)

    elif varname == 'BMB_v2':
        #Create BMB colormap
        vmax = 100
        vmin = -10
        linthresh = .3
        linscale = .25
        fracpos = (np.log10(vmax/linthresh)+linscale)/(np.log10(vmax/linthresh)+np.log10(-(vmin/linthresh))+2*linscale)
        nneg = np.int_((1-fracpos)*256)
        colors1 = plt.get_cmap('cmo.ice')(np.linspace(0,1.,nneg+1))
        colors2 = plt.get_cmap('inferno')(np.linspace(0., 1, 256-nneg-1))
        colors = np.vstack((colors1, colors2))

        cmap = mpl.colors.LinearSegmentedColormap.from_list('my_colormap', colors)
        norm = mpl.colors.SymLogNorm(linthresh, vmin=vmin, vmax=vmax, linscale=linscale)
        
    elif varname == 'Hi':
        cmap = copy(plt.get_cmap('cmo.ice'))
        norm = mpl.colors.Normalize(vmin=0,vmax=4000,clip=True)

    elif varname == 'Hib':
        cmap = copy(plt.get_cmap('cmo.deep'))
        norm = mpl.colors.Normalize(vmin=-500,vmax=0,clip=True)
    
    elif varname == 'Hs':
        cmap = copy(plt.get_cmap('cmo.ice'))
        norm = mpl.colors.Normalize(vmin=0,vmax=1000,clip=True)

    elif varname == 'H_lad':
        cmap = copy(plt.get_cmap('gist_stern'))
        norm = mpl.colors.Normalize(vmin=1,vmax=300,clip=True)

    elif varname == 'T_lad':
        cmap = copy(plt.get_cmap('cmo.thermal'))
        norm = mpl.colors.Normalize(vmin=-2,vmax=1,clip=True)

    elif varname == 'S_lad':
        cmap = copy(plt.get_cmap('cmo.haline'))
        norm = mpl.colors.Normalize(vmin=32,vmax=34,clip=True)

    elif varname in ['U_lad', 'V_lad']:
        cmap = copy(plt.get_cmap('cmo.balance'))
        norm = mpl.colors.Normalize(vmin=-.2,vmax=.2,clip=True)

    elif varname in ['Uabs_lad']:
        cmap = copy(plt.get_cmap('Greens'))
        norm = mpl.colors.Normalize(vmin=0,vmax=.2,clip=True)

    elif varname in ['uabs_surf', 'uabs_vav']:
        cmap = copy(plt.get_cmap('turbo'))
        norm = mpl.colors.Normalize(vmin=0,vmax=2000,clip=True)
        #norm = mpl.colors.LogNorm(vmin=1.,vmax=3000,clip=True)
    
    else:
        print(f'ERROR: no colormap available yet for {varname}, add one to colormaps.py')
        return

    return cmap,norm