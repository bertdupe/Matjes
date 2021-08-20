import matplotlib.pyplot as plt
import numpy as np
import os


def unfold_other_Bz_k(lat,klat,k_in):
  #unfold the k-points and dos to other brillouinzones defined by lat
  k_out=[]
  for l in lat:
      k_out.append(k_in+klat@l)
  k_out=np.array(k_out).flatten()
  shape=k_out.shape
  k_out.shape=[len(k_out)//3,3]
  return k_out

def unfold_other_Bz_dos(lat,dos_in):
  #unfold the k-points and dos to other brillouinzones defined by lat
  shape=np.array(dos_in.shape)
  dos_out=[]
  for l in lat:
      dos_out.append(dos_in)
  dos_out=np.array(dos_out)
  shape[0]=shape[0]*len(lat)
  dos_out.shape=shape
  return dos_out

def plot_dos(k,klat,dos,name='dos.png',title='',vext=None,cmap='viridis',lc='white',clabel=''):
  fig,ax=plt.subplots()
  if vext:
      im=ax.tricontourf(k[:,0],k[:,1],dos,20,vmin=vext[0],vmax=vext[1],cmap=cmap)
  else:
      im=ax.tricontourf(k[:,0],k[:,1],dos,20,cmap=cmap)
  ax.set_aspect(1)
  ax.set_xlim([-2/3*klat[0,0],2/3*klat[0,0]])
  ax.set_ylim([-klat[1,0],klat[1,0]])
  ax.set_xticks([])
  ax.set_yticks([])
  ax.set_title(title)
  cbar=plt.colorbar(im,ax=ax)
  cbar.set_label(clabel)
  
  #add boundary 
  edges=np.array([[1,1,0],[2,-1,0],[1,-2,0],[-1,-1,0],[-2,1,0],[-1,2,0],[1,1,0]])/3
  edges= edges @ klat.T
  ax.plot(edges[:,0],edges[:,1],color=lc,lw=4)
  size_cm=[13.0,9.0] 
  fig.set_size_inches(size_cm[0]/2.54,size_cm[1]/2.54)
  fig.savefig(name,dpi=150)
  plt.close(fig)


Ncell=3*3
klat=np.array([[ 6.283185,   3.627599,  -0.000000 ],
               [ 6.283185,  -3.627599,  -0.000000 ],
               [-0.000000,  -0.000000,   6.283185 ]]).T/(Ncell**0.5)

lat=np.array([[0.0,0.0,0.0],[-1.0,-1.0,0.0],[-1.0, 0.0,0.0],[ 0.0,-1.0,0.0],])


dat=np.loadtxt('fermidos.dat')
kpts=unfold_other_Bz_k(lat,klat,dat[:,:3])
dos=unfold_other_Bz_dos(lat,dat[:,3])

plot_dos(kpts,klat,dos,name='fermisurf.png',title='fermi surface',vext=None,cmap='inferno',lc='white',clabel='')
