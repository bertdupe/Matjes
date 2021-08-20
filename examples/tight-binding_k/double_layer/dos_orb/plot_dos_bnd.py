import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

fig,ax =plt.subplots()
dat_sum=np.loadtxt("dos_k_nc.dat")
im=ax.scatter(dat_sum[:,0],dat_sum[:,1],edgecolor=None,c='black',s=16.0,label='total DOS')
dat=[]
dat.append(np.loadtxt("dos_k_nc_bnd_001.dat"))
dat.append(np.loadtxt("dos_k_nc_bnd_002.dat"))
legend=['Atom 1 orb. 1', 'Atom 2 orb. 1']
for i,d in enumerate(dat):
    im=ax.scatter(d[:,0],d[:,1],s=5.0,edgecolor=None,label=legend[i])
ax.set_xlim(np.amin(dat[0][:,0]),np.amax(dat[0][:,0]))
ax.set_xlabel(r"E-E$_F$ [eV]")
ax.set_ylabel(r"proj. DOS [arb. units]")
ax.set_yticks([])
ax.set_ylim(bottom=0.0)
ax.legend(loc='upper left')
fig.savefig("dos_orbital.png",dpi=150)
