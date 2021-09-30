import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
import os


def plot_spin(position,spin,shape,name):
    fontsmall=11
    fontlarge=12
    size_cm=[18.0,9.0]
    fig, axes = plt.subplots(ncols=3,figsize=(size_cm[0]/2.54,size_cm[1]/2.54))
    #fig.patch.set_alpha(0.0)
    cmap = 'seismic'
    #cmap = 'bwr'
    x_step=np.abs(position[1,0,0]-position[0,0,0])
    y_step=np.abs(position[0,1,1]-position[0,0,1])
    gridshape=shape[0:2]
    x=position[:,:,0].flatten()
    y=position[:,:,1].flatten()
    x=x-(np.amax(x)+np.amin(x))/2
    y=y-(np.amax(y)+np.amin(y))/2
    label=['$m_x$','$m_y$','$m_z$']
    for i in range(0,3):
      ax=axes[i]
      ax.set_aspect(1)
      ax.set_facecolor('black')
      extent=np.array([np.amin(x),np.amax(x),np.amin(y),np.amax(y)])
      gridsize=[gridshape[0]-1,gridshape[1]]
      diff=np.array([0.0,0.0,-y_step,y_step])
      extent=extent+diff
      im=ax.hexbin(x,y,C=spin[:,:,i].flatten(),gridsize=gridsize,cmap=cmap,vmin=-1,vmax=1,aa=False,edgecolor=None,extent=extent,zorder=10)
      ax.set_xticks([])
      ax.set_yticks([])
      [i.set_linewidth(3.0) for i in ax.spines.values()]
    
      ax.text(0.05, 0.95, label[i], size=12, rotation=0.,color='black',
           ha="left", va="top",transform=ax.transAxes,  multialignment='center',
           bbox=dict(boxstyle="round",pad=0.15,lw=1.0,
                     ec=(0.0, 0.0, 0.0),
                     fc=(1.0, 1.0, 1.0),
                     )
           )
    
    
    fig.subplots_adjust(left=0.05,right=0.95,hspace=0.03,wspace=0.03)
    cbar=fig.colorbar(im,ax=axes,ticks=[-1,0,1],pad=0.03)
    cbar.ax.tick_params(labelsize=fontsmall,direction='in',pad=2,length=2,width=0.5)
    cbar.set_label(r'magnetization',fontsize=fontlarge,labelpad=-8)
    fig.savefig(name,dpi=300)#,transparent=True)
    plt.close(fig)


with open('input', 'r') as f:
    for line in f:
        if 'Nsize' == line[0:5]:
            shape=list(map(int,(line[5:].split())))
shape=[shape[1],shape[0],3]

position=np.loadtxt('positions.dat')
position.shape=shape
spindats=[]
for file in os.listdir("./"):
    if 'magnetic_' in file and file.endswith(".dat"):
        spindats.append(file)

for f in spindats:
    spin=np.loadtxt(f)
    spin.shape=shape
    name=f[:-4]+'.png'
    plot_spin(position,spin,shape,name)
