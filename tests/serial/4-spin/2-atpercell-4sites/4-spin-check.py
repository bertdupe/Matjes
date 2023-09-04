import math
import sys
import numpy as np
import os
import linecache
import subprocess
import time
import array

print("4-spin check script")

Nmax = 20
counter = 0
step = math.pi / (Nmax-1)/2
#Matjes=os.system("/Users/bertranddupe/Documents/GitHub/Matjes/build_cmake/Matjes")
#Matjesprocess = subprocess.Popen(['/Users/bertranddupe/Documents/GitHub/Matjes/build_cmake/Matjes'])
#rm=os.system("rm -f magnetic_init.dat EM.dat")

Energy=np.empty((Nmax,2),dtype=float,order='C')

original_stdout = sys.stdout # Save a reference to the original standard output

while counter < Nmax:
    os.system("rm -f magnetic_init.dat EM.dat")
    Energy[counter,0]=counter*step
    
    c = math.cos(counter*step)
    s = math.sin(counter*step)
    Mi = np.array( [s*math.sqrt(3.0)/2.0  , s/2.0                  ,c ] )
    Mj = np.array( [s*0.5                 ,-s*math.sqrt(3.0)/2.0   ,c ] )
    Mk = np.array( [-s*0.5                , s*math.sqrt(3.0)/2.0   ,c ] )
    Ml = np.array( [-s*math.sqrt(3.0)/2.0 ,-s/2.0                  ,c ] )

    maginit = open('magnetic_init.dat', 'w')
    for i in Mi:
        print(i,end=" ",file=maginit)
    for i in Mj:
        print(i,end=" ",file=maginit)
    print(file=maginit)
    for i in Mk:
        print(i,end=" ",file=maginit)
    for i in Ml:
        print(i,end=" ",file=maginit)
    print(file=maginit)
    maginit.close()

    os.system("/Users/bertranddupe/Documents/GitHub/Matjes/build_cmake/Matjes")

    f = open('EM.dat', 'r')
    header = f.readline()
    data = f.readline()
    line = data.strip()
    value = line.split()
    f.close()

    Energy[counter,1] = float(value[2])

    print()
    print(counter,Energy[counter])
    print()
    time.sleep(1)

    counter = counter+1

resultats=open('resultats.dat','w')
i=0
while i < Nmax:
    for j in Energy[i]:
        print(j,end=" ",file=resultats)
    print(file=resultats)
    i=i+1
