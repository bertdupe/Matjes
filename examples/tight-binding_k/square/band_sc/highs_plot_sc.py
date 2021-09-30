import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

c="red"
cmap = colors.LinearSegmentedColormap.from_list("incr_alpha", [(0, (*colors.to_rgb(c),0)), (1, c)])
dat=np.loadtxt("highs_plot_sc.dat")
fig,ax =plt.subplots()
im=ax.scatter(dat[:,0],dat[:,1],c=dat[:,2],s=1.0,edgecolor=None,cmap=cmap,vmin=0.0,vmax=1.0,zorder=10)
ax.set_xlim(np.amin(dat[:,0]),np.amax(dat[:,0]))
ax.set_xticks([
  0.00000000E+00,
  0.31415927E+01,
  0.62831853E+01,
  0.10726068E+02,
])
ax.set_xticklabels([
'g', 
'X', 
'M', 
'g', 
])
ax.set_ylabel(r"E-E$_F$ [eV]")
ax.grid(axis="x",color="black",zorder=1)
cbar=plt.colorbar(im,ax=ax)
cbar.set_label("electron weight")
plt.show()
fig.savefig("highs_plot_sc.png")
