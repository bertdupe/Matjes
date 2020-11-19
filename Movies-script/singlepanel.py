#input files
inp = "inp"
dyna = "dyna.in"
lattice = "lattice.in"
EM = "EM.dat"




#movie
      		                       
movie_start=input("At what timestep do you want to start? Default is 1: ")		#movie starts at this timestep         
lattice_shift_x = input("Give the x-coodinates with respect to the center of the lattice. For a shift to the left choose positive sign, for the right choose negative sign, for no shift enter 0: ")
lattice_shift_y = input("Give the y-coordinates with respect to the center of the lattice(at y=0). For a downward shift choose positive sign, for upwards choose negative sign, for no shift enter 0: ")



#settings
bottom_panel_plot = raw_input("Choose what to plot for the bottom panel. Enter abbreviations for Vorticity in x/y/z direction(VX,VY,VZ), Topological charge(TC), Spin lattice(SL) or B effective(BE): ")
top_panel_plot_1 = raw_input("What is the first thing you want to plot on the top panel? Enter  Total topological charge, Torque in eV, Magnetization, Magnetic field in T, Energy (meV), Energy or Temperature: ")
top_panel_plot_2 = raw_input("What is the second thing you want to plot? Choose one other thing from above. For no second axis enter none:" )
size = float(raw_input("Zoom size? 1 equals to 100% and 0 to 0%: "))                                                                                           
cb_max_start=1						                                                                                                      #taking average for bottom panel from this timestep until the end	#default: 1
						

#definition of the libraries
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import pylab



#######################################################
# definition of usefull function
#######################################################
def get_shape(filename):                                                                                                 #reads the required data for the spin lattice from lattice.in
    f = open(filename, 'r')
    for line in f:
        if "Nsize" in line:
            a = line.split()
            row = int(a[1])
            column = int(a[2])
    shape = (row, column)
    lines = row*column

    return shape, lines




def load_topo(filename, plot_data, shape):                                                                              # reads the desired data (see input bottom_panel_plot) from topomap.dat
    try:
        data = np.loadtxt(filename, unpack=True)
    except:
        print 'Bad topomap data given'

    xValue = data[0]
    yValue = data[1]

    density = data[5]
    if plot_data == "TC":
        density = data[5]
    elif plot_data == "VX":
        density = data[2]
    elif plot_data == "VY":
        density = data[3]
    elif plot_data == "VZ":
        density = data[4]
   

    return xValue.reshape(shape), yValue.reshape(shape), density.reshape(shape)





def top_panel(filename, plot_data):                                                                                   # reads the desired data for the top panel (see input top_panel_plot_1/2) from EM.dat
    try:
        data = np.loadtxt(filename, unpack=True)
    except:
        print 'Something is wrong with EM.dat'
   
   
    top_value=data[10] 
    if plot_data == "Total topological charge":
        top_value = data[10]        ###Total topological charge
    elif plot_data == "Torque in eV":
        top_value = data[11]    ###searching
    elif plot_data == "Magnetic field in T":
        top_value = data[18] ###searching
    elif plot_data == "Energy (meV)":
	    top_value = data[1]*1000.0
    elif plot_data == "Energy":
	    top_value = data[1]
    elif plot_data == "Magnetization":
	    top_value = data[2]
    elif plot_data == "Temperature":
	    top_value = data[13]
      


    return top_value.tolist()
    



def get_input(filename):                                                                                              #reads required data from dyna.in
    f = open(filename, 'r')
    for line in f:
        if "timestep" in line:
            a = line.split()
            timestep = float(a[1])
        elif "duration" in line:
            a = line.split()
            duration = int(a[1])

    return timestep, duration



##########three functions for spinse_dat##################
def load_spins(filename, shape):                                                                                       #reads Spinse.dat files (position and angle orientation) and creates data arrays
    try:
        data = open(filename, 'r')
    except:
        print 'Bad spin data given'

    theta, phi = [],[]
    X,Y = [],[]
    for line in data:
        a = line.split()
        theta.append(float(a[1][:-1]))
        phi.append(float(a[2][:-1]))
        X.append(float(a[3][:-1]))
        Y.append(float(a[4][:-1]))

    X, Y = np.array(X), np.array(Y)

    return np.array(theta), np.array(phi), X.reshape(shape), Y.reshape(shape)




def degree_to_rad(angle1, angle2):                                                                                     # converts degree to radiant                                                                                 
    rad1 = angle1 * np.pi/180
    rad2 = angle2 * np.pi/180

    return rad1, rad2




def spin_components(theta, phi, shape):                                                                               # creates arrays of spin components 
    S_x = np.cos(phi) * np.sin(theta)
    S_y = np.sin(phi) * np.sin(theta)
    S_z = np.cos(theta)

    return S_x.reshape(shape), S_y.reshape(shape), S_z.reshape(shape)

##########################################    



##########three functions for Ffield_dat(has rarely been used...but the functions are the same as above for the spinse files above)###################################################################################################

def load_Ffield(filename, shape):
    try:
       data = open(filename, 'r')
    except:
       print 'Bad Ffield data given'

    theta1, phi1 = [], []
    X1, Y1 = [], []
    for line in data:
	a = line.split()
	theta1.append(float(a[2][:-1]))
	phi1.append(float(a[3][:-1]))
	X1.append(float(a[4][:-1]))
	Y1.append(float(a[5][:-1]))

    X1, Y1 = np.array(X1), np.array(Y1)

    return np.array(theta1), np.array(phi1), X1.reshape(shape), Y1.reshape(shape)



def degree_to_rad2(angle1, angle2):
    rad1 = angle1 * np.pi/180
    rad2 = angle2 * np.pi/180

    return rad1, rad2



def spin_components2(theta1, phi1, shape):
    S_x1 = np.cos(phi1) * np.sin(theta1)
    S_y1 = np.sin(phi1) * np.sin(theta1)
    S_z1 = np.cos(theta1)

    return S_x1.reshape(shape), S_y1.reshape(shape), S_z1.reshape(shape)

##########################################################################################################################################    




def resize(x, y, rho, X, Y, S_x, S_y, S_z, size, shape):                                                                        
    qm = int(size*shape[0])
    nshape = (qm, qm)

    xar = np.zeros(nshape)
    yar = np.zeros(nshape)
    rhoar = np.zeros(nshape)
    Xar = np.zeros(nshape)
    Yar = np.zeros(nshape)
    S_xar = np.zeros(nshape)
    S_yar = np.zeros(nshape)
    S_zar = np.zeros(nshape)

    for k in range(qm):
        for l in range(qm):
            id1 = (len(x)-qm)/2 + k
            id2 = (len(x)-qm)/2 + l
            xar[k][l] = x[id1][id2]
            yar[k][l] = y[id1][id2]
            rhoar[k][l] = rho[id1][id2]
            Xar[k][l] = X[id1][id2]
            Yar[k][l] = Y[id1][id2]
            S_xar[k][l] = S_x[id1][id2]
            S_yar[k][l] = S_y[id1][id2]
            S_zar[k][l] = S_z[id1][id2]

    return xar, yar, rhoar, Xar, Yar, S_xar, S_yar, S_zar




#program starts
shape, lines = get_shape(lattice)
timestep, duration = get_input(dyna)


#counter
number_of_files = int(duration)
steps = input("Choose steps (plot frequency) between the data. Choose 1, 10, 20 ...: ")                                                #this function skips lots of files depending on what number has been typed and plots only the desired
img=1                                                                                                                                  #data


#search for maximum value ( bottom panel)
if bottom_panel_plot =="TC" or bottom_panel_plot =="VX" or bottom_panel_plot=="VY" or bottom_panel_plot=="VZ":                         #creates a list of the maximum values from the topomap data for the contourf level function;
  cb_max = []
  j = cb_max_start
 
  while j < duration:
     x, y, density = load_topo("topomap"+str(j)+".dat", bottom_panel_plot, shape)
     cb_max.append(np.max(density))
     print str(j)+" of "+str(duration)+" files checked."
     j = j + steps   #optional
	    
     scale = np.mean(cb_max)                                                                                                           #takes the average from the maximum values for the contourf levels



#top value

top_value1 = top_panel(EM, top_panel_plot_1)                                                                                           #forms a timelist for the top panel data
top_value2 = top_panel(EM, top_panel_plot_2)
timelist = []
for time in range(duration):
    timelist.append((time+1)*timestep)




#loop for all data
n = movie_start
while n < number_of_files:

 #create figure
       
	fig =plt.figure(figsize=plt.figaspect(1.20))                                                                                  #sets the aspect ratio of the whole figure
	
	
# syntax pos_x pos_y size_x size_y
        if bottom_panel_plot=="SL" or bottom_panel_plot=="BE":
           ax1 = plt.axes([0.10, 0.1, 0.90, 0.60])                                                                                    # generates axes and size for the top and bottom panel
	   plt.yticks(fontsize=10)
	   plt.xticks(fontsize=10)
           ax3 = plt.axes([0.10, 0.8, 0.80, 0.15])
	   plt.yticks(fontsize=10)
           plt.xticks(fontsize=10)
	  
	  

        elif bottom_panel_plot=="TC" or bottom_panel_plot=="VX" or bottom_panel_plot=="VY" or bottom_panel_plot=="VZ":   
	   ax1 = plt.axes([0.15, 0.1, 0.88, 0.6])
	   plt.yticks(fontsize=10)
	   plt.xticks(fontsize=10)
	   ax3 = plt.axes([0.10, 0.8, 0.85, 0.15])
	   plt.yticks(fontsize=10)
	   plt.xticks(fontsize=10)
       

        #bottom values
        if bottom_panel_plot=="TC" or bottom_panel_plot=="VX" or bottom_panel_plot=="VY" or bottom_panel_plot=="VZ":
	  x_array, y_array, density_array = load_topo("topomap"+str(n)+".dat", bottom_panel_plot, shape)                             #loads the desired data for the bottom panel
	if bottom_panel_plot=="SL":
          theta_deg, phi_deg, X, Y = load_spins('Spinse_'+str(n)+'.dat', shape)
          theta, phi = degree_to_rad(theta_deg, phi_deg)
	  S_x, S_y, S_z = spin_components(theta, phi, shape)
        elif bottom_panel_plot =="BE":
	  theta1_deg, phi1_deg, X1, Y1 = load_Ffield('Ffield_'+str(n)+'.dat', shape)
	  theta1, phi1 = degree_to_rad2(theta1_deg, phi1_deg)
	  S_x1, S_y1, S_z1 = spin_components2(theta1, phi1, shape)
      
	       
      
         


        #top panel
	ax3.set_xlim(0, duration*timestep)                                                                                           #settings for the axes of the top panel
	ax3.set_xlabel("time in fs", fontsize=10)
	ax3.set_ylabel(top_panel_plot_1, fontsize=10)
        if top_panel_plot_1 == "Total topological charge":
	   ax3.set_ylim(min(top_value1), max(top_value1))
        elif top_panel_plot_1 =="Torque in eV":
	   ax3.set_yticks([min(top_value1), max(top_value1)])
        elif top_panel_plot_1 =="Magnetic field in T":
	   ax3.set_yticks([min(top_value1), max(top_value1)])	
        elif top_panel_plot_1 =="Energy":
	   ax3.set_ylim(min(top_value1), max(top_value1))
        elif top_panel_plot_1 =="Magnetization":
	   ax3.set_ylim(min(top_value1), max(top_value1))
        elif top_panel_plot_1 =="Energy (meV)": 
	   ax3.set_ylim(min(top_value1), max(top_value1))
	elif top_panel_plot_1 =="Temperature":
	   ax3.set_ylim(min(top_value1), max(top_value1))
	   
        plt.plot(timelist, top_value1, 'r', linewidth=1.5)
	plt.plot(timelist[n-1], top_value1[n-1], 'ko', markersize=4)

     
        ax4= ax3.twinx()                                                                                                             #generates the second y axis on the right side 
	ax4.set_ylabel(top_panel_plot_2, fontsize=10)
	plt.yticks(fontsize=10)
        if top_panel_plot_2 == "Total topological charge":
           ax4.set_ylim(min(top_value2), max(top_value2))
	   plt.plot(timelist, top_value2, 'b', linewidth=1.5)
	   plt.plot(timelist[n-1], top_value2[n-1], 'ko', markersize=4)
        elif top_panel_plot_2 =="Torque in eV":
	   ax4.set_yticks([min(top_value2), max(top_value2)])
	   plt.plot(timelist, top_value2, 'b', linewidth=1.5)
	   plt.plot(timelist[n-1], top_value2[n-1], 'ko', markersize=4)
        elif top_panel_plot_2 =="Magnetic field in T":
	   ax4.set_yticks([min(top_value2), max(top_value2)])
	   plt.plot(timelist, top_value2, 'b', linewidth=1.5)
	   plt.plot(timelist[n-1], top_value2[n-1], 'ko', markersize=4)
        elif top_panel_plot_2 =="Energy":
	   ax4.set_yticks([min(top_value2), max(top_value2)])
	   plt.plot(timelist, top_value2, 'b', linewidth=1.5)
	   plt.plot(timelist[n-1], top_value2[n-1], 'ko', markersize=4)
        elif  top_panel_plot_2=="Magnetization":
	   ax4.set_ylim(min(top_value2), max(top_value2))
	   plt.plot(timelist, top_value2, 'b', linewidth=1.5)
	   plt.plot(timelist[n-1], top_value2[n-1], 'ko', markersize=4)
        elif top_panel_plot_2 =="Energy (meV)":
            ax4.set_ylim(min(top_value2), max(top_value2))
	    plt.plot(timelist, top_value2, 'b', linewidth=1.5)
	    plt.plot(timelist[n-1], top_value2[n-1], 'ko', markersize=4)
        elif top_panel_plot_2=="Temperature":
	    ax4.set_ylim(min(top_value2), max(top_value2))
	    plt.plot(timelist, top_value2, 'b', linewidth=1.5)
	    plt.plot(timelist[n-1], top_value2[n-1], 'ko', markersize=4)
        elif top_panel_plot_2=="none":  
	    ax4.set_frame_on(False)
	    ax4.axes.get_yaxis().set_visible(False)

     
        
		

        		




		


	#bottom panel (for topomap data)
	if bottom_panel_plot=="TC" or bottom_panel_plot=="VX" or bottom_panel_plot=="VY" or bottom_panel_plot=="VZ":                                 
	  ax1.set_xlim([(1-size)*(np.max(x_array)/2) + lattice_shift_x, (1+size)*(np.max(x_array)/2) + lattice_shift_x])                                 #shifts and zooms into the lattice according to the input
	  ax1.set_ylim([np.min(y_array)*(size) + lattice_shift_y, np.max(y_array)*(size) + lattice_shift_y])
          ax1.set_xlabel('x-position in nm')
	  ax1.set_ylabel('y-position in nm')
	  ax1.set_title(bottom_panel_plot, fontsize=12)

	  levels = MaxNLocator(nbins=100).bin_boundaries(-scale, scale)                                                                                  #contour function
	  cmap = plt.get_cmap('seismic')
	  contour = ax1.contourf(x_array, y_array, density_array, levels=levels, cmap=cmap)
									
	  cb = fig.colorbar(contour, ax=ax1)
	  cb.set_ticks([0])

          #bottom panel (for Spinse and Ffield data)
        elif bottom_panel_plot== "SL":
          ax1.set_xlim([(1-size)*(np.max(X)/2) + lattice_shift_x ,(1+size)*(np.max(X)/2) + lattice_shift_x ])
	  ax1.set_ylim([np.min(Y)*(size) + lattice_shift_y, np.max(Y)*(size)+ lattice_shift_y])
	  ax1.set_xlabel("X-position in nm")
	  ax1.set_title("Magnetization (Mz)", fontsize=12)


#	andere colorbar, die bei mara funktioniert, aber nicht bei stephan
#          jetmap = plt.get_cmap('jet')                                                                                                                  #generates a colorbar for the spin direction
#	  Mesh = ax1.pcolormesh(S_z, cmap=jetmap, vmin=-180, vmax=0)
#	  cm = fig.colorbar(Mesh, ax = ax1)
#	  cm.set_ticks([0, -90, -180])
#	  plt.cla()

          

          for sp1 in range(len(X)):                                                                                                                    #colorcode for the spins
	      for sp2 in range(len(X)):

	        r = 0.5*np.sqrt(S_x[sp1][sp2]**2+S_y[sp1][sp2]**2)
		if S_z[sp1][sp2] >= 0:
		   if r == 0.0:
			r = 0.01
		   if S_x[sp1][sp2] == 0.0:
		      S_x[sp1][sp2] = 0.01
		   if S_y[sp1][sp2] == 0.0:
		      S_y[sp1][sp2] = 0.01
		   ax1.arrow(X[sp1][sp2], Y[sp1][sp2], r*S_x[sp1][sp2], r*S_y[sp1][sp2], length_includes_head=True, head_width=r*2/3, head_length=r ,color=(S_z[sp1][sp2], np.sqrt(S_x[sp1][sp2]**2+S_y[sp1][sp2]**2), 0))
		else:
		   ax1.arrow(X[sp1][sp2], Y[sp1][sp2], r*S_x[sp1][sp2], r*S_y[sp1][sp2], length_includes_head=True, head_width=r*2/3, head_length=r ,color=(0, np.sqrt(S_x[sp1][sp2]**2+S_y[sp1][sp2]**2), -S_z[sp1][sp2]))

#########################################same function for the Ffield data (has rarely been used...)##############################################################

	elif bottom_panel_plot=="BE":
	  ax1.set_xlim([(1-size)*(np.max(X1)/2) + lattice_shift_x, (1+size)*(np.max(X1)/2) + lattice_shift_x])
          ax1.set_ylim([np.min(Y1)*(size) + lattice_shift_y, np.max(Y1)*(size) + lattice_shift_y])
	  ax1.set_xlabel("x-position in nm")
	  ax1.set_title("Magnetization (M-z)", fontsize=12)

	  jetmap = plt.get_cmap('jet')
	  Mesh = ax1.pcolormesh(S_z1, cmap=jetmap, vmin=-180, vmax=0)
	  cm = fig.colorbar(Mesh, ax = ax1)
	  cm.set_ticks([0, -90, -180])
	  plt.cla()

	  for sp1 in range(len(X1)):
	      for sp2 in range(len(X1)):

	        r = 0.5*np.sqrt(S_x1[sp1][sp2]**2+S_y1[sp1][sp2]**2)
		if S_z1[sp1][sp2] >= 0:
	           ax1.arrow(X1[sp1][sp2], Y1[sp1][sp2], r*S_x1[sp1][sp2], r*S_y1[sp1][sp2], length_includes_head=True, head_width=r*2/3, head_length=r ,color=(S_z1[sp1][sp2], np.sqrt(S_x1[sp1][sp2]**2+S_y1[sp1][sp2]**2), 0))
	        else:   
		   ax1.arrow(X1[sp1][sp2], Y1[sp1][sp2], r*S_x1[sp1][sp2], r*S_y1[sp1][sp2], length_includes_head=True, head_width=r*2/3, head_length=r ,color=(0, np.sqrt(S_x1[sp1][sp2]**2+S_y1[sp1][sp2]**2), -S_z1[sp1][sp2]))

###############################################################################################################################################################

		  
		  
		  
	pc = plt.pcolor(S_z, cmap='jet',vmin=-180, vmax=0)
	cm = fig.colorbar(pc, ax=ax1)
	cm.set_ticks([0,-90,-180])
	plt.cla()	  
		  
		  
		  
		  
		  
		  
		  
	if bottom_panel_plot=="SL" or bottom_panel_plot=="BE":                                                                                #settings for image saving
	   fig.savefig('spin_'+str(img).zfill(6)+'.png', dpi= 100, orientation='portrait')         
           print ('spin_'+str(img)+'.dat was converted.')
        elif bottom_panel_plot=="TC" or bottom_panel_plot=="VX" or bottom_panel_plot=="VY" or bottom_panel_plot=="VZ":
	   fig.savefig('topomap_'+str(img).zfill(6)+'.png', dpi=100, orientation='portrait', close=True)
	   print ('topomap_'+str(img)+'.dat was converted.')

        plt.close()
        #counter
        img = img+1
        n = n + steps
