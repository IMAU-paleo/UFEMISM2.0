import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

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
        if (len(name)==1):
            # Reached last non-placeholder name
            break
        if len(name)>maxlev:
            # Not looking this deep, skip
            continue
        names.append(name)
        tcomps.append(ds.tcomp[r].values)
    return names,tcomps

def plot_resources(names,tcomps,Nlevels,savename='test.png',toplevel=0,toproutine='UFEMISM_program'):

    plt.rcParams["axes.prop_cycle"] = plt.cycler("color", plt.get_cmap('Dark2')(np.linspace(0,1,8)))
    plt.rcParams["figure.subplot.left"] = .01
    plt.rcParams["figure.subplot.right"] = .99
    plt.rcParams["figure.subplot.bottom"] = .01
    plt.rcParams["figure.subplot.top"] = .99

    width = 14
    height = 8

    fig,ax = plt.subplots(1,1,figsize=(width,height))

    y0 = np.zeros((Nlevels,1))
    y1 = np.zeros((Nlevels,1))
    x = np.linspace(0,1,Nlevels+1)
    x0 = x[:-1]
    x1 = x[1:]
    lev = -1

    # Get ttot
    ttot = 0
    for i,tcomp in enumerate(tcomps):
        if len(names[i]) == toplevel+2:
            ttot += tcomp
    print(ttot)

    for i,tcomp in enumerate(tcomps):

        if len(names[i]) == toplevel+1:
            ttot = tcomp

        if len(names[i]) < toplevel+1:
            # Not enough levels
            continue
        elif names[i][toplevel] != toproutine:
            # Doesn't fall within toproutine, so skip
            continue
        elif len(names[i]) > toplevel+Nlevels+1:
            # Too deep, don't show
            continue
        else:
            
            levnew = len(names[i])-toplevel-2
            if levnew == lev:
                y0[levnew] = y1[levnew]
            elif levnew > lev:
                y0[levnew] = y0[levnew-1]
            else:
                y0[levnew] = y1[levnew]

            lev = levnew
            y1[lev] = y0[lev] + tcomp/ttot
            if 'gather_to_all' in names[i][lev+1]:
                ax.fill_between([x0[lev],x1[lev]],y0[lev],y1[lev],lw=0,color='.3')
            elif 'multiply_CSR_matrix' in names[i][lev+1]:
                ax.fill_between([x0[lev],x1[lev]],y0[lev],y1[lev],lw=0,color='tab:red')
            else:
                ax.fill_between([x0[lev],x1[lev]],y0[lev],y1[lev],lw=0)

            if tcomp/ttot > .01:
                # At least 1 percent of total
                rotation = 180/3.14* np.arctan((y1[lev][0]-y0[lev][0])/(x1[lev]-x0[lev])*height/width)
                if rotation < 10:
                    rotation = 0
                ax.text((x0[lev]+x1[lev])/2,(y0[lev]+y1[lev])/2,f'{names[i][lev+1]} ({100*tcomp/ttot:.1f}%)',rotation=rotation,va='center',ha='center',fontsize=6)

    ax.set_xlim([0,1])
    ax.set_ylim([0,1])
    ax.set_xticks([])
    ax.set_yticks([])
    #ax.set_title(f'{baseroutine}')
    plt.savefig(savename)
    #plt.close()
    print(f'Created {savename}')