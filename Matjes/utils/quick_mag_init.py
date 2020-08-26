import numpy as np
import matplotlib.pyplot as plt

pos=np.loadtxt('positions.dat')
mag_init=np.loadtxt('SpinSTM_start.dat')

fig,axes=plt.subplots(1,3,sharex='all',sharey='all')
pltdict={'vmin':-1.0,'vmax':1.0,'cmap':'coolwarm'}
title=['spin-x','spin-y','spin-z']
for i in range(0,3):
    ax=axes[i]
    im=ax.scatter(pos[:,0],pos[:,1],c=mag_init[:,i],**pltdict)
for i,ax in enumerate(axes):
    ax.set_title(title[i])
    ax.set_aspect(1)
    ax.set_xlim(np.amin(pos[:,0]),np.amax(pos[:,0]))
    ax.set_ylim(np.amin(pos[:,1]),np.amax(pos[:,1]))
    ax.set_xticks([])
    ax.set_yticks([])
plt.colorbar(im,ax=axes)
plt.show()
