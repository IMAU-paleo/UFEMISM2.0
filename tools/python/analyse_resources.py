import xarray as xr
import numpy as np

def print_line(ds,n,lim=1):
    c = 0 #counter of /

    frac = 100*ds.tcomp[0,n].values/np.max(ds.tcomp[0,:].values)
    if frac<.01:
        return
    vals = ds.routine_names_encoded[:,n].values
    name = []
    for i in range(len(ds.name_length)):
        if vals[i] == -1 and c>lim:
            break
        elif vals[i]<10 and c>lim:
            name.append(str(vals[i]))
        elif vals[i]<37 and c>lim:
            name.append(chr(ord('`')+vals[i]-10))
        elif vals[i]<63 and c>lim:
            name.append(chr(ord('@')+vals[i]-36))
        elif vals[i]==63 and c>lim:
            name.append('_')
        elif vals[i]==64:
            if c>lim:
                name.append('/')
            c += 1
        elif vals[i]==65 and c>lim:
            name.append('(')
        elif vals[i]==66 and c>lim:
            name.append(')')
    strname = ''.join(name)
    #if strname[:13] != 'run_BMB_model':
    #    return
    if strname == 'subroutine_placeholder' or c>7:
        return
    print(f'{frac: 7.2f}% {strname}')

def print_resources(inputfolder,run):
    #Read data
    fname = f'{inputfolder}/{run}/resource_tracking.nc'
    ds = xr.open_dataset(fname)

    for n in range(len(ds.n_routines)):
        print_line(ds,n)

    ds.close()
