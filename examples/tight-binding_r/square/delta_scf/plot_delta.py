import numpy as np
import matplotlib.pyplot as plt

delta_dat=np.loadtxt("delta_onsite_final.dat",skiprows=2)
pos=np.loadtxt("positions.dat")
fig, ax =plt.subplots()
delta=np.linalg.norm(delta_dat,axis=1)*1e3
pltargs={'vmin':min(0,np.amin(delta)),'vmax':np.amax(delta),'cmap':'viridis_r','edgecolor':'black'}

c = ax.pcolor(np.reshape(delta[:],[11,11]),**pltargs)
cbar=fig.colorbar(c,ax=ax,label=r'|$\Delta_{i}$| [meV]')
#ax.set_title(r'$\Delta$')
ax.set_aspect(1)
ax.set_xticks([])
ax.set_yticks([])
fig.savefig('delta_spacedependence.png',dpi=150)

