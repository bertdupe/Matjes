import matplotlib.pyplot as plt
import numpy as np
import os
import sys
sys.path.insert(1, '/home/pbuhl/scripts/matjes/python_scripts/tight-binding-fs')
from func import *


def plot_dos_arr(k,klat,dosarr,titlearr=None,plot_shape=[2,2],name='dos.png',title='',vext=None,cmap='viridis',lc='white',clabel=''):
  fig,axes=plt.subplots(*plot_shape)
  ax_rav=np.ravel(axes)
  for i in range(0,len(dosarr)):
    ax=ax_rav[i]
    dos=dosarr[i]
    if titlearr:
        ax.set_title(titlearr[i])
    if vext:
        im=ax.tricontourf(k[:,0],k[:,1],dos,20,vmin=vext[0],vmax=vext[1],cmap=cmap)
    else:
        im=ax.tricontourf(k[:,0],k[:,1],dos,20,cmap=cmap)

  for ax in ax_rav:
    ax.set_aspect(1)
    ax.set_xlim([-2/3*klat[0,0],2/3*klat[0,0]])
    ax.set_ylim([-klat[1,0],klat[1,0]])
    ax.set_xticks([])
    ax.set_yticks([])
  for ax in ax_rav[len(dosarr):]:
    ax.set_axis_off()
  fig.suptitle(title)
  cbar=plt.colorbar(im,ax=axes)
  cbar.set_label(clabel)

  #add boundary
  edges=np.array([[1,1,0],[2,-1,0],[1,-2,0],[-1,-1,0],[-2,1,0],[-1,2,0],[1,1,0]])/3
  edges= edges @ klat.T
  for ax in ax_rav:
      ax.plot(edges[:,0],edges[:,1],color=lc,lw=4)
  size_cm=[26.0,20.0]
  fig.set_size_inches(size_cm[0]/2.54,size_cm[1]/2.54)
  fig.savefig(name,dpi=150)
  plt.close(fig)

Norb=2
Ncell=1

klat=np.array([[ 6.283185,   3.627599,  -0.000000 ],
               [ 6.283185,  -3.627599,  -0.000000 ],
               [-0.000000,  -0.000000,   6.283185 ]]).T/(Ncell**0.5)
lat=np.array([[0.0,0.0,0.0],[-1.0,-1.0,0.0],[-1.0, 0.0,0.0],[ 0.0,-1.0,0.0],])


dat=np.loadtxt('fermidos_allproj.dat')
proj=dat[:,3:]
shape=proj.shape
shape=np.array([shape[0],Ncell,Norb])
proj.shape=shape
sum_orbs=np.sum(proj,axis=1)


kpts=unfold_other_Bz_k(lat,klat,dat[:,:3])
dos=unfold_other_Bz_dos(lat,sum_orbs)

dosarr=[dos[:,i] for i in range(0,dos.shape[1])]

plot_dos_arr(kpts,klat,dosarr,titlearr=None,plot_shape=[2,1],name='fermisurf_proj.png',title='projected fermi surfaces',vext=None,cmap='inferno',lc='white',clabel='')
