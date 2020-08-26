#This python script plots the output files of the TB branch in the data_state_dist directory in the case of a rectangular lattice.
#It need to be in the data_state_dist directory with the occ_dfermi_*.dat and/or occ_fermi_*.dat files with the position.dat, dos_r_sc.dat  and input files in the directory above.

import matplotlib.pyplot as plt
import numpy as np
import os

def get_lattice(filename):
    f = open(filename, 'r')
    lines=f.readlines()
    Nlines=len(lines)
    lattice=np.zeros([3,3])
    nspin=1
    for i in range(0,Nlines):
        if "Nsize" in lines[i]:
            tmp = lines[i].split()
            shape=np.array([int(st) for st in tmp[1:4]])
        if lines[i][0:4]== 'alat':
            tmp = lines[i].split()
            alat=np.array([float(st) for st in tmp[1:4]])
        if lines[i][0:7]== 'lattice':
            lattice[:,0]=np.array([float(st) for st in lines[i+1].split()])
            lattice[:,1]=np.array([float(st) for st in lines[i+2].split()])
            lattice[:,2]=np.array([float(st) for st in lines[i+3].split()])
        if 'activate_mag_TB' in lines[i]:
            tmp = lines[i].split()
            if 't' in tmp[1].lower():
                nspin=2
    lattice_out=lattice*alat* shape
    return [lattice_out,shape,nspin]

def get_extema(intarr):
    ext_f=[np.inf,-np.inf]
    ext_df=[np.inf,-np.inf]
    for i in intarr:
        file_fermi='occ_fermi_{:03d}.dat'.format(i)
        file_dfermi='occ_dfermi_{:03d}.dat'.format(i)
        f=np.loadtxt(file_fermi,skiprows=1)
        ext_f[0]=min(ext_f[0],np.amin(f))
        ext_f[1]=max(ext_f[1],np.amax(f))
        df=np.loadtxt(file_dfermi,skiprows=1)
        ext_df[0]=min(ext_df[0],np.amin(df))
        ext_df[1]=max(ext_df[1],np.amax(df))
    return [ext_f,ext_df]

def plot_states(i_dat,pos,E,dos,extema,shape,plotdict):
    file_fermi='occ_fermi_{:03d}.dat'.format(i_dat)
    file_dfermi='occ_dfermi_{:03d}.dat'.format(i_dat)
    bgcol='darkgrey'
    with open(file_fermi) as f: 
        Ef = float(f.readline())

    #prepare plotting axes
    fig,axes=plt.subplots(2,3, gridspec_kw={'width_ratios': [1, 2,2]})
    fig.subplots_adjust(left=0.1,hspace=0.08,wspace=0.08,top=0.95)
    i_col=0
    
    #plot DOS
    gs=axes[0,i_col].get_gridspec()
    for ax in axes[0:, 0]:
        ax.remove()
    axbig=fig.add_subplot(gs[0:,i_col])
    ax=axbig
    ax.plot(dos,E,zorder=10,c='black')
    ax.scatter(dos,E,zorder=10,c='black')
    ax.set_xlim(left=0.0)
    ax.plot([0.0, np.amax(dos)*100], [Ef,Ef], 'red', lw=2,zorder=5)
    ax.set_ylim([np.amin(E),np.amax(E)])
    ax.set_xlabel('DOS',labelpad=-1)
    ax.set_ylabel('E[eV]',labelpad=-3)
    i_col=i_col+1

    text_label=['spin\nup','spin\ndown'] 
    bbox=dict(boxstyle="round",ec=(0.0, 0.0, 0.0),fc=(0.8, 0.8, 0.8),alpha=1.0 )
    text_dict={'x':0.05,'y':0.98,'bbox':bbox,'size': 8, 'ha':"left", 'va':"top",'zorder':100,'multialignment':'center'}
    #plot fermi distribution
    if True:
        occ=np.loadtxt(file_fermi,skiprows=1)
        occ.shape=shape
        axs=axes[:,i_col]
        for i in range(0,shape[2]):
          im=axs[i].hexbin(pos[0],pos[1],C=occ[:,:,i].flatten(),vmin=extrema[0][0],vmax=extrema[0][1],cmap='inferno',**plotdict)
        for i,ax in enumerate(axs):
         ax.set_aspect(1)
         ax.set_facecolor(bgcol)
         ax.set_xticks([])
         ax.set_yticks([])
         ax.set_xlim([np.amin(pos[0]),np.amax(pos[0])])
         ax.set_ylim([np.amin(pos[1]),np.amax(pos[1])])
         ax.text(s=text_label[i],transform=ax.transAxes,**text_dict)
        cbar=fig.colorbar(im,ax=axs[:])
        cbar.set_label(r'occupation fermi-dirac',labelpad=1)
        i_col=i_col+1

    #plot derivative fermi distribution
    if True:
        dE=np.loadtxt(file_dfermi,skiprows=1)
        dE.shape=shape
        axs=axes[:,i_col]
        contour_args={'origin':'lower','cmap':'inferno_r','interpolation':'nearest','vmin':extrema[1][0],'vmax':extrema[1][1]}
        for i in range(0,shape[2]):
         im = axs[i].hexbin(pos[0],pos[1],C=dE[:,:,i].flatten(),vmin=extrema[1][0],vmax=extrema[1][1],cmap='inferno_r',**plotdict)
        for i,ax in enumerate(axs):
         ax.set_facecolor(bgcol)
         ax.set_aspect(1)
         ax.set_xticks([])
         ax.set_yticks([])
         ax.set_xlim([np.amin(pos[0]),np.amax(pos[0])])
         ax.set_ylim([np.amin(pos[1]),np.amax(pos[1])])
         ax.text(s=text_label[i],transform=ax.transAxes,**text_dict)
        cbar=fig.colorbar(im,ax=axs)
        cbar.set_label(r'fermi-dirac derivative',labelpad=1)
        i_col=i_col+1

    #finalize plotting
    for ax in np.ravel(axes[:,i_col:]):
        ax.remove()
    fig.savefig('state_dist_{:03d}.png'.format(i_dat),dpi=150)
    plt.close(fig)


#prepare position and dos data
[_,shape,nspin]=get_lattice('../input')
Nsize=shape[:2]
shape=[*Nsize,nspin]
pos=np.loadtxt('../positions.dat')
pos.shape=[*Nsize,3]
x_step=np.abs(pos[1,0,0]-pos[0,0,0])
y_step=np.abs(pos[0,1,1]-pos[0,0,1])
gridshape=shape[0:2]
x=pos[:,:,0].flatten()
y=pos[:,:,1].flatten()
x=x-(np.amax(x)+np.amin(x))/2
y=y-(np.amax(y)+np.amin(y))/2
extent=np.array([np.amin(x),np.amax(x),np.amin(y),np.amax(y)])
gridsize=[gridshape[0]-1,gridshape[1]]
diff=np.array([0.0,0.0,-y_step,y_step])
extent=extent+diff
pos=[x,y]
plotdict={"gridsize":gridsize,"aa":False,"edgecolor":None,"extent":extent,"zorder":10}



dos_data=np.loadtxt('../dos_r_sc.dat')
E=dos_data[:,0]
dos=dos_data[:,1]

#get maximal used indice of data file
files=os.listdir('./')
n_ext=[1,1]
for file in sorted(files,reverse = True):
 file='./'+file
 if not os.path.isfile(file):
  continue
 if 'occ_dfermi_' in file:
  n_ext[1]=int(file[13:16])+1
  break
 elif 'occ_fermi_' in file:
  n_ext[1]=int(file[12:15])+1
  break

#get range of color values for plot
extrema=get_extema(list(range(n_ext[0], n_ext[1])))
#actual plotting
for i in range(n_ext[0],n_ext[1]):
 plot_states(i,pos,E,dos,extrema,shape,plotdict)
