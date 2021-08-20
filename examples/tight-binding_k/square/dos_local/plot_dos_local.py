import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import os

fig,ax =plt.subplots()
#dat_sum=np.loadtxt("dos_k_nc.dat")
#im=ax.scatter(dat_sum[:,0],dat_sum[:,1],edgecolor=None,c='black',s=16.0,label='total DOS')
files=[]
for file in os.listdir("./"):
   if 'dos_k_nc_bnd_' in file:
    files.append(file)
files.sort()

legend=[
" site 3 3",
" site 4 3",
" site 5 3"]

dat=[]
for f in files:
    dat.append(np.loadtxt(f))
#legend=['Atom 1 orb. 1', 'Atom 2 orb. 1']
dos_tot=np.loadtxt('dos_k_nc.dat')
im=ax.scatter(dos_tot[:,0],dos_tot[:,1]/25,s=16.0,color='black',edgecolor=None,label='total DOS/ number states')
for i,d in enumerate(dat):
    im=ax.scatter(d[:,0],d[:,1],s=5.0,edgecolor=None,label=legend[i])
ax.set_xlim(np.amin(dat[0][:,0]),np.amax(dat[0][:,0]))
ax.set_xlabel(r"E-E$_F$ [eV]")
ax.set_ylabel(r"proj. DOS [arb. units]")
ax.set_yticks([])
ax.set_ylim(bottom=0.0)
ax.legend(loc='upper left')
fig.set_size_inches(10, 6)
fig.savefig("dos_local.png",dpi=150)
