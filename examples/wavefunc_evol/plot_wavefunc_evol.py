import numpy as np
import matplotlib.pyplot as plt

dat=np.loadtxt('wavefunc_evol.dat')
dat_norm=dat[:,5:]
time=dat[:,0]
fig,ax=plt.subplots()
im=ax.imshow(dat_norm,aspect='auto',extent=[0.0,1.0,np.amin(time),np.amax(time)],origin='lower',interpolation='nearest')
ax.set_ylabel('time [fs]')
ax.set_xlabel('position')
ax.set_xticks([])
cbar=plt.colorbar(im,ax=ax)
cbar.set_label(r'$|\phi_i|^2$')
fig.set_size_inches(4, 6)
fig.savefig("wavefunc_evol.png",dpi=150)

