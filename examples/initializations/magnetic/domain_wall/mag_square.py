import numpy as np
import matplotlib.pyplot as plt
import os
def plot_spin(position,spin,shape,name):
    fig,axes=plt.subplots(3,1,sharex='all',sharey='all')
    pltdict={'vmin':-1.0,'vmax':1.0,'cmap':'coolwarm','origin':'lower','interpolation':'nearest'}
    title=[r'$m_x$',r'$m_y$',r'$m_z$']
    
    text_args=dict(x=0.02, y=0.90, size=8, rotation=0.,color='black',
                   ha="left", va="top", multialignment='center',
                   bbox=dict(boxstyle="round",pad=0.15,lw=1.0,
                    ec=(0.0, 0.0, 0.0),
                    fc=(1.0, 1.0, 1.0),
                    )
               )
    spin.shape=shape
    for i in range(0,3):
        ax=axes[i]
        im=ax.imshow(spin[:,:,i],**pltdict)
    for i,ax in enumerate(axes):
        ax.text(s=title[i], transform=ax.transAxes,**text_args)
        ax.set_aspect(1)
        ax.set_xticks([])
        ax.set_yticks([])
    fig.subplots_adjust(wspace=0.05,hspace=0.05,top=0.95,bottom=0.05,left=0.01,right=0.99)
    plt.colorbar(im,ax=axes)
    fig.set_size_inches(10,2.5)
    fig.savefig(name,dpi=150)

with open('input', 'r') as f:
    for line in f:
        if 'Nsize' == line[0:5]:
            shape=list(map(int,(line[5:].split())))

shape=[shape[1],shape[0],3]

position=np.loadtxt('positions.dat')
spindats=[]
for file in os.listdir("./"):
    if 'magnetic_' in file and file.endswith(".dat"):
        spindats.append(file)

for f in spindats:
    spin=np.loadtxt(f)
    name=f[:-4]+'.png'
    plot_spin(position,spin,shape,name)

