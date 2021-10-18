import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

fig,ax =plt.subplots()

dat=np.loadtxt('EM.dat')

ax.errorbar(dat[:,0],dat[:,1]*1e3,yerr=dat[:,2]*1e3,c='C0',fmt='o',markersize=0.5, capsize=3,label='energy',zorder=100)
ax_r=ax.twinx()
ax_r.scatter(dat[:,0],dat[:,3],c='C2',s=20.0,marker='o',label='C',zorder=50)
ax_r.plot   (dat[:,0],dat[:,3],c='C2',zorder=40)
ax.set_xlim(np.amin(dat[:,0]),np.amax(dat[:,0]))
ax_r.set_ylim(bottom=0.0)

ax.spines['left'].set_color('C0')
ax.yaxis.label.set_color('C0')
ax.tick_params(axis='y', colors='C0')

ax_r.spines['left'].set_color('C2')
ax_r.yaxis.label.set_color('C2')
ax_r.tick_params(axis='y', colors='C2')

ax.set_xlabel("temperature [K]")
ax.set_ylabel("energy per site [meV]")
ax_r.set_ylabel("C")
ax.set_zorder(1)
ax.set_frame_on(False)
fig.set_size_inches(6, 4)
fig.savefig("MC.png",dpi=150)
