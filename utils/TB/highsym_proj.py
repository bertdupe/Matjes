import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors

c="red"
ylim=[-1.0,1.0]
dat=np.loadtxt('highs_plot_proj.dat')
subplot_size=[4,5]

dat=dat[np.logical_and(dat[:,1]>ylim[0],dat[:,1]<ylim[1]),:]

kdist=dat[:,0]
bnd=dat[:,1]
weights=dat[:,2:]
vmax=np.amax(weights)


cmap = colors.LinearSegmentedColormap.from_list("incr_alpha", [(0, (*colors.to_rgb(c),0)), (1, c)])
cmap='Reds'

fig,axes =plt.subplots(*subplot_size,sharex=True,sharey=True)
ax_rav=np.ravel(axes)
for i in range(0,weights.shape[1]):
    ax=ax_rav[i]
    im=ax.scatter(kdist,bnd,s=2.0,c=weights[:,i],vmin=0.0,vmax=vmax,cmap=cmap)
xlim=[np.amin(kdist),np.amax(kdist)]
#for ax in ax_rav:
ax.set_xlim(xlim)
ax.set_ylim(ylim)
ax.set_xticks([
  0.00000000E+00,
  0.36275986E+01,
  0.57219939E+01,
  0.99107841E+01,
  0.13052377E+02,
  0.16679975E+02,
  0.18774371E+02,
  0.22963161E+02,
])
ax.set_xticklabels([
"G",
"M",
"K",
"G",
"A",
"H",
"L",
"G",
])
for ax in axes[:,0]:
    ax.set_ylabel(r"E-E$_F$ [eV]")
for ax in ax_rav:
    ax.grid(axis="x",color="black",zorder=1)
for ax in ax_rav[weights.shape[1]:]:
    ax.set_axis_off()
cbar=fig.colorbar(im,ax=axes)
cbar.set_label('wannier state projection')
plt.show() 
#ax.set_xlim(np.amin(dat[:,0]),np.amax(dat[:,0]))
