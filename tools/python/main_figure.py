import os
import glob
import xarray as xr
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.gridspec as gridspec

from colormaps import *
from utils import *

class Figure(object):
    """ Figure of subplots of Ufe/Lad output """

    def __init__(self, figsize, directory='../figures/'):

        self.directory = directory

        self.orientation = 'horizontal'

        self.cbar_ratio = 30
        self.figsize = figsize

        self.fields = {}

    def __repr__(self):
        return f"Figure()"

    def __str__(self):
        return f"Figure"

    def set_orientation(self,orientation):

        if orientation in ['horizontal','h']:
            self.orientation = 'horizontal'
        elif orientation in ['vertical','v']:
            self.orientation = 'vertical'
        else:
            raise ValueError('Invalid option for orientation')

        return f"Set figure orientation to {self.orientation}"

    def add_field(self,Field):

        self.fields[Field.name] = Field
        return 

    def delete_field(self,fieldname):

        del self.fields[fieldname]
        return

    def print_fields(self):

        print(self.fields.keys())
        return

    def make(self,figname,dxmin=0,dxmax=0,dymin=0,dymax=0, add_gl=True):

        fig = plt.figure(figsize=self.figsize,constrained_layout=True)

        if self.orientation=='horizontal':
            nrows = 2
            ncols = len(self.fields)
            spec = gridspec.GridSpec(nrows,ncols,figure=fig,height_ratios=[self.cbar_ratio,1],hspace=0.01,wspace=0.01)
            for k,key in enumerate(self.fields.keys()):
                self.fields[key].ax = fig.add_subplot(spec[0,k])
                self.fields[key].cax = fig.add_subplot(spec[1,k])
        else:
            nrows = len(self.fields)
            ncols = 2
            spec = gridspec.GridSpec(nrows,ncols,figure=fig,width_ratios=[self.cbar_ratio,1],hspace=0.01,wspace=0.01)
            for k,key in enumerate(self.fields.keys()):
                self.fields[key].ax = fig.add_subplot(spec[k,0])
                self.fields[key].cax = fig.add_subplot(spec[k,1])

        for k,key in enumerate(self.fields.keys()):
            field = self.fields[key]
            if field.mask is not None:
                field.get_mcoll()
                im = field.ax.add_collection(field.mcoll)
            field.get_pcoll()
            im = field.ax.add_collection(field.pcoll)
            cbar = plt.colorbar(im,cax=field.cax,orientation=self.orientation)
            cbar.set_label(field.varname)

            if add_gl:
                if not field.Timeframe.got_gl:
                    field.Timeframe.get_gl()
                field.ax.plot(field.Timeframe.gl[0,:],field.Timeframe.gl[1,:],c='k',lw=.5)

            field.ax.set_xlim([field.xmin+dxmin,field.xmax-dxmax])
            field.ax.set_ylim([field.ymin+dymin,field.ymax-dymax])
            field.ax.set_aspect(1)
            field.ax.set_xticks([])
            field.ax.set_yticks([])

        fig.suptitle(f'Year {field.Timeframe.time:.0f}')

        fullfigname = os.path.join(self.directory,f'{figname}.png')
        plt.savefig(fullfigname, bbox_inches = 'tight', pad_inches = 0,dpi=450)
        print(f'Created {figname}')
        plt.close()

        return


class Field(object):
    """ Field to include in subplot """

    def __init__(self, Timeframe, varname, mask=None):

        self.Timeframe = Timeframe
        self.t = self.Timeframe.t
        self.varname = varname
        self.name = f"{varname}_{self.t}"
        self.mask = mask

        self.xmin = self.Timeframe.ds.xmin
        self.xmax = self.Timeframe.ds.xmax
        self.ymin = self.Timeframe.ds.ymin
        self.ymax = self.Timeframe.ds.ymax

        self.get_cmap()

        self.get_data()

    def __repr__(self):
        return f"Field({repr(self.Timeframe)}, {self.name})"

    def __str__(self):
        return f"Field {self.name} at time {self.Timeframe.t}"

    def get_data(self):
        """ Extract data values """


        if self.varname == 'Uabs_lad':
            uvar = self.Timeframe.ds['U_lad']
            vvar = self.Timeframe.ds['V_lad']
            self.data = (uvar**2+vvar**2)**.5
        elif self.varname[:3] == 'BMB':
            self.data = -self.Timeframe.ds['BMB']
        else:
            try:
                self.data = self.Timeframe.ds[self.varname]
            except KeyError:
                print(f"ERROR: {self.varname} not in Timeframe")
                return

        if self.mask is not None:
            if 'vi' in self.data.dims:
                if not self.Timeframe.got_mask:
                    self.Timeframe.get_mask()

                if self.mask == 'shelf':
                    self.data = xr.where(self.Timeframe.mask == 4, self.data, np.nan)
                elif self.mask == 'ice':
                    self.data = xr.where([x in [3,4] for x in self.Timeframe.mask], self.data, np.nan)
                else:
                    raise ValueError('Invalid option for mask')
            elif 'ti' in self.data.dims:
                # Just mask out zero values, no clean option yet
                self.data = xr.where(self.data==0, np.nan, self.data)
            else:
                print(f'ERROR: variable {varname} is not on vertices or triangles')

            message = f"Extracted {self.varname} with mask {self.mask}"
        else:
            message = f"Extracted {self.varname}"

        return message

    def get_cmap(self):
        """ Get colormap info """
        cmap,norm = get_cmap(self.varname)

        self.cmap = cmap
        self.norm = norm

        return

    def get_pcoll(self):
        """ Get patch collection """

        #Check type (voronoi / triangle)
        if 'vi' in self.data.dims:
            if not self.Timeframe.Mesh.got_voronois:
                self.Timeframe.Mesh.get_voronois()
            self.pcoll = PatchCollection(self.Timeframe.Mesh.voronois,cmap=self.cmap,norm=self.norm)
        elif 'ti' in self.data.dims:
            if not self.Timeframe.Mesh.got_triangles:
                self.Timeframe.Mesh.get_triangles()
            self.pcoll = PatchCollection(self.Timeframe.Mesh.triangles,cmap=self.cmap,norm=self.norm)
        else:
            print(f'ERROR: variable {varname} is not on vertices or triangles')

        # Fill array
        self.pcoll.set_array(self.data.values)

        return

    def get_mcoll(self):
        """ Get mask collection """

        mcmap = 'ocean'
        mnorm = mpl.colors.Normalize(vmin=1.6, vmax=3.1, clip=True)

        if not self.Timeframe.Mesh.got_voronois:
            self.Timeframe.Mesh.get_voronois()

        self.mcoll = PatchCollection(self.Timeframe.Mesh.voronois,cmap=mcmap,norm=mnorm)

        # Fill array
        self.mcoll.set_array(self.Timeframe.mask)

        return
