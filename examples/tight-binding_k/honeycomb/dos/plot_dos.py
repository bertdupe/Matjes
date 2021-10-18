import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import os

fig,ax =plt.subplots()

files=[]
for file in os.listdir("./"):
   if file[0:4]=='dos_' and file[-4:]=='.dat':
    files.append(file)
files.sort()

label=[]
dats=[]
for f in files:
    dats.append(np.loadtxt(f))
    #label.append(f[4:-4])
    label.append(f)

for i,dat in enumerate(dats):
    im=ax.scatter(dat[:,0],dat[:,1],edgecolor=None,s=10.0,label=label[i])
    ax.plot(dat[:,0],dat[:,1])
ax.set_xlim(np.amin(dat[:,0]),np.amax(dat[:,0]))
ax.set_xlabel(r"E-E$_F$ [eV]")
ax.set_ylabel(r"DOS [arb. units]")
ax.set_yticks([])
ax.set_ylim(bottom=0.0)
fig.set_size_inches(6, 4)
ax.legend()
fig.savefig("dos.png",dpi=150)
