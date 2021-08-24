import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import os

fig,ax =plt.subplots()
#dat_sum=np.loadtxt("dos_k_nc.dat")
#im=ax.scatter(dat_sum[:,0],dat_sum[:,1],edgecolor=None,c='black',s=16.0,label='total DOS')
files=[]
for file in os.listdir("./"):
   if 'dos_k_sc_bnd_' in file:
    files.append(file)
files.sort()

legend=[
"site 6 6 spin up",
"site 6 6 spin down",
"site 7 6 spin up",
"site 7 6 spin down",
"site 7 7 spin up",
"site 7 7 spin down"]

dat=[]
for f in files:
    dat.append(np.loadtxt(f))
#legend=['Atom 1 orb. 1', 'Atom 2 orb. 1']
dos_tot=np.loadtxt('dos_k_sc.dat')
im=ax.scatter(dos_tot[:,0]*1e3,dos_tot[:,1]/11/11,s=5.0,color='black',edgecolor=None,label='total DOS/ number states')
ax.plot(dos_tot[:,0]*1e3,dos_tot[:,1]/11/11,color='black')
for i,d in enumerate(dat):
    im=ax.scatter(d[:,0]*1e3,d[:,1],s=5.0,edgecolor=None,label=legend[i])
    ax.plot(d[:,0]*1e3,d[:,1])
ax.set_xlim(np.amin(dat[0][:,0]),np.amax(dat[0][:,0]))
ax.set_xlim(-500,500)
ax.set_xlabel(r"E-E$_F$ [meV]")
ax.set_ylabel(r"proj. DOS [arb. units]")
ax.set_yticks([])
ax.set_ylim(bottom=0.0)
ax.legend(loc='upper left')
fig.set_size_inches(10, 6)
fig.savefig("dos_local.png",dpi=150)
