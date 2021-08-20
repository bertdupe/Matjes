import numpy as np
import matplotlib.pyplot as plt

pos=np.loadtxt('positions.dat')
mag_init=np.loadtxt('magnetic_start.dat')

fig,axes=plt.subplots(1,3,sharex='all',sharey='all')
pltdict={'vmin':-1.0,'vmax':1.0,'cmap':'coolwarm'}
title=['spin-x','spin-y','spin-z']

text_args=dict(x=0.05, y=0.95, size=8, rotation=0.,color='black',
               ha="left", va="top", multialignment='center',
               bbox=dict(boxstyle="round",pad=0.15,lw=1.0,
                ec=(0.0, 0.0, 0.0),
                fc=(1.0, 1.0, 1.0),
                )
           )
shape=[60,60]
for i in range(0,3):
    ax=axes[i]
    #im=ax.scatter(pos[:,0],pos[:,1],c=mag_init[:,i],**pltdict)
    mag_init.shape=[*shape,3]
    im=ax.imshow(mag_init[:,:,i],origin='lower',**pltdict)
for i,ax in enumerate(axes):
    ax.text(s=title[i], transform=ax.transAxes,**text_args)
    ax.set_aspect(1)
#    ax.set_xlim(np.amin(pos[:,0]),np.amax(pos[:,0]))
#    ax.set_ylim(np.amin(pos[:,1]),np.amax(pos[:,1]))
    ax.set_xticks([])
    ax.set_yticks([])
fig.subplots_adjust(wspace=0.05,hspace=0.01,top=0.95,bottom=0.05,left=0.01,right=0.99)
plt.colorbar(im,ax=axes)
fig.set_size_inches(10,2.5)
fig.savefig('magnetic_texture_init.png',dpi=150)
