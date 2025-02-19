import xarray as xr
import numpy as np

def get_name(rname):
    lev = 0
    dname = []
    name = {}
    # Loop over individual values
    for i in range(len(rname)):
        if rname[i] == -1:
            # End of name reached
            break
        elif rname[i]< 10:
            dname.append(str(rname[i]))
        elif rname[i]<37:
            dname.append(chr(ord('`')+rname[i]-10))
        elif rname[i]<63:
            dname.append(chr(ord('@')+rname[i]-36))
        elif rname[i]==63:
            dname.append('_')
        elif rname[i]==64:
            # Found a '/', store this level
            name[lev] = ''.join(dname)
            # Go to next level
            lev += 1
            dname = []
            continue
        elif rname[i]==65:
            dname.append('(')
        elif rname[i]==66:
            dname.append(')')
    # Store last level
    name[lev] = ''.join(dname)
    return name

def get_allnames(ds,maxlev=7):
    names = []
    tcomps = []
    for r in range(len(ds.n_routines)):
        name = get_name(ds.routine_names_encoded[:,r].values)
        if len(name)==1:
            # Reached last non-placeholder name
            break
        if len(name)>maxlev:
            # Not looking this deep, skip
            continue
        names.append(name)
        tcomps.append(ds.tcomp[r].values)
    return names,tcomps

def plot_resources(names,tcomps,baselevel,baseroutine,toplevel=0,toproutine='UFEMISM_program'):

    plt.rcParams["axes.prop_cycle"] = plt.cycler("color", plt.get_cmap('Dark2')(np.linspace(0,1,8)))

    fig,ax = plt.subplots(1,1)

    for i,tcomp in enumerate(tcomps):
        if len(names[i]) < baselevel+1:
            continue
        if names[i][toplevel] != toproutine:
            continue
        if names[i][baselevel] != baseroutine:
            continue
        elif len(names[i])==baselevel+1:
            # Get total tcomp for x-scaling
            x1 = 0
            xscale = tcomp
        elif len(names[i])==baselevel+2:
            # Get width
            x0 = x1
            x1 += tcomp/xscale
            y1 = 0
            yscale = tcomp
            if tcomp/xscale > .01:
                # At least 1 percent of total
                ax.text((x0+x1)/2-.02,0,names[i][baselevel+1],rotation=-60,va='top',ha='left')
        elif len(names[i])==baselevel+3:
            # Get height
            y0 = y1
            y1 += tcomp/yscale

            ax.fill_betweenx([y0,y1],x0,x1)
            if tcomp/xscale > .01:
                if (y1-y0>x1-x0):
                    # Vertical
                    ax.text((x0+x1)/2,(y0+y1)/2,names[i][baselevel+2],rotation=-90,va='center',ha='center')
                else:
                    # Horizontal
                    ax.text((x0+x1)/2,(y0+y1)/2,names[i][baselevel+2],va='center',ha='center')
    ax.axis("off")
    ax.set_title(f'{baseroutine}')
    plt.show()
