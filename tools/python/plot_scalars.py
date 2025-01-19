import os
import xarray as xr
import matplotlib.pyplot as plt


def plot_scalars(inputfolder,outputfolder,run):
    #Read data
    fname = f'{inputfolder}/{run}/scalar_output_ANT_00001.nc'
    ds = xr.open_dataset(fname)

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

    ds.close()

    #Make outputfolder if necessary
    if not os.path.isdir(outputfolder):
        os.makedirs(outputfolder)

    #Save figure
    plt.savefig(f'{outputfolder}/scalars.pdf')
