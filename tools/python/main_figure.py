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

    def add_field(self, Timeframe, varname, mask=None):
        field = Field(Timeframe, varname, mask=mask)

        print(f'Added {varname} with hash {hash(field)}')

        self.fields[hash(field)] = field
        return 

    def add_diff(self, Timeframe1, varname1, Timeframe2, varname2, name='diff', mask=None, vmax=None, colmap='cmo.diff'):
        field1 = Field(Timeframe1, varname1)
        field2 = Field(Timeframe2, varname2)
        field = DiffField(field1, field2, name=name, mask=mask, vmax=vmax, colmap=colmap)

        print(f'Added difference field {name} with hash {hash(field)}')

        self.fields[hash(field)] = field
        return 

    def delete_field(self,fieldhash):

        del self.fields[fieldhash]
        return

    def print_fields(self):

        print(self.fields.keys())
        return

    def make(self,figname,dxmin=0,dxmax=0,dymin=0,dymax=0, add_gl=True, add_time=False):

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

        if add_time:
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
        self.Mesh = self.Timeframe.Mesh
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
            if not self.Mesh.got_voronois:
                self.Mesh.get_voronois()
            self.pcoll = PatchCollection(self.Mesh.voronois,cmap=self.cmap,norm=self.norm)
        elif 'ti' in self.data.dims:
            if not self.Mesh.got_triangles:
                self.Mesh.get_triangles()
            self.pcoll = PatchCollection(self.Mesh.triangles,cmap=self.cmap,norm=self.norm)
        else:
            print(f'ERROR: variable {varname} is not on vertices or triangles')

        # Fill array
        self.pcoll.set_array(self.data.values)

        return

    def get_mcoll(self):
        """ Get mask collection """

        mcmap = 'ocean'
        mnorm = mpl.colors.Normalize(vmin=1.6, vmax=3.1, clip=True)

        if not self.Mesh.got_voronois:
            self.Mesh.get_voronois()

        self.mcoll = PatchCollection(self.Mesh.voronois,cmap=mcmap,norm=mnorm)

        # Fill array
        self.mcoll.set_array(self.Timeframe.mask)

        return

class DiffField(object):
    """ Difference field between two Fields """

    def __init__(self, Field1, Field2, name, mask, vmax, colmap):

        self.Field1 = Field1
        self.Field2 = Field2
        self.Timeframe = self.Field1.Timeframe

        self.name = name
        self.varname = 'diff'
        self.mask = mask

        self.xmin = self.Field1.xmin
        self.xmax = self.Field1.xmax
        self.ymin = self.Field1.ymin
        self.ymax = self.Field1.ymax

        self.vmax = vmax
        self.colmap = colmap

        self.check_compatibility()

        self.get_cmap()

        self.get_data()

    def __repr__(self):
        return f"DiffField({repr(self.Field1)}, {repr(self.Field2)}, {self.name})"

    def __str__(self):
        return f"DiffField between {self.Field1.name} at time {self.Field1.Timeframe.t} and {self.Field1.name} at time {self.Field1.Timeframe.t}"

    def check_compatibility(self):
        """ Check whether difference field can be computed """

        assert self.xmin == self.Field2.xmin, 'ERROR: xmin is not equal between input Fields'
        assert self.xmax == self.Field2.xmax, 'ERROR: xmax is not equal between input Fields'
        assert self.ymin == self.Field2.ymin, 'ERROR: ymin is not equal between input Fields'
        assert self.ymax == self.Field2.ymax, 'ERROR: ymax is not equal between input Fields'

        if 'vi' in self.Field1.data.dims:
            assert 'vi' in self.Field2.data.dims, 'ERROR: cannot take difference between voronois and triangles'
            assert (self.Field1.Timeframe.ds['V'] == self.Field2.Timeframe.ds['V']).all(), 'ERROR: V is not equal between input Fields'
        elif 'ti' in self.Field1.data.dims:
            assert 'ti' in self.Field2.data.dims, 'ERROR: cannot take difference between voronois and triangles'
            assert (self.Field1.Timeframe.ds['Tri'] == self.Field2.Timeframe.ds['Tri']).all(), 'ERROR: ti is not equal between input Fields'

        print('Input Fields are compatible')
        self.Mesh = self.Field1.Timeframe.Mesh
        return

    def get_cmap(self):
        """ Define cmap and norm """
        self.cmap = plt.get_cmap(self.colmap)

        if self.vmax is None:
            self.dmax = max(abs(self.Field1.data-self.Field2.data).values)
            if self.dmax == 0:
                print('WARNING: both fields are equal')
                vmax = 1
            else:
                vmax = self.dmax
        else:
            vmax = self.vmax

        self.norm = mpl.colors.Normalize(vmin=-vmax,vmax=vmax,clip=True)
        return

    def get_data(self):
        """ Get data and add mask if necessary """

        self.data = self.Field1.data - self.Field2.data

        if self.mask is not None:
            if 'vi' in self.data.dims:
                if not self.Timeframe.got_mask:
                    self.Timeframe.get_mask()

                if self.mask == 'shelf':
                    self.data = xr.where(self.Timeframe.mask == 4, self.data, np.nan)
                elif self.mask == 'sheet':
                    self.data = xr.where(self.Timeframe.mask == 3, self.data, np.nan)
                elif self.mask == 'ice':
                    self.data = xr.where([x in [3,4] for x in self.Timeframe.mask], self.data, np.nan)
                else:
                    raise ValueError('Invalid option for mask')
            elif 'ti' in self.data.dims:
                # Just mask out zero values, no clean option yet
                self.data = xr.where(self.data==0, np.nan, self.data)
            else:
                print(f'ERROR: something went wrong with vertices and triangles')

            message = f"Extracted {self.name} with mask {self.mask}"
        else:
            message = f"Extracted {self.name}"

        return message

    def get_pcoll(self):
        """ Get patch collection """

        #Check type (voronoi / triangle)
        if 'vi' in self.data.dims:
            if not self.Mesh.got_voronois:
                self.Mesh.get_voronois()
            self.pcoll = PatchCollection(self.Mesh.voronois,cmap=self.cmap,norm=self.norm)
        elif 'ti' in self.data.dims:
            if not self.Mesh.got_triangles:
                self.Mesh.get_triangles()
            self.pcoll = PatchCollection(self.Mesh.triangles,cmap=self.cmap,norm=self.norm)
        else:
            print(f'ERROR: variable {varname} is not on vertices or triangles')

        # Fill array
        self.pcoll.set_array(self.data.values)

        return

    def get_mcoll(self):
        """ Get mask collection """

        mcmap = 'ocean'
        mnorm = mpl.colors.Normalize(vmin=1.6, vmax=3.1, clip=True)

        if not self.Mesh.got_voronois:
            self.Mesh.get_voronois()

        self.mcoll = PatchCollection(self.Mesh.voronois,cmap=mcmap,norm=mnorm)

        # Fill array
        self.mcoll.set_array(self.Timeframe.mask)

        return
