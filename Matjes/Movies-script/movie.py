#input files
inp = "inp"
dyna = "dyna.in"
lattice = "lattice.in"
EM = "EM.dat"




#movie
fps = 12      		#Caution: This must be changed in both scripts: movie.py and topomovie.sh
movie_start=1		#movie starts at this timestep          #default: 1
length = 120            #seconds        length of complete movie




#settings
left_panel_plot = "Topological charge" 			#other possibilities: "Vorticity in x/y/z direction" -> compare with load_topo function in line 70 etc
top_panel_plot = "Total topological charge"		#other possibilities: "Magnetic field" or "Torque" -> compare with top_panel function in line 90 etc
size = 0.25                                      	# 1 equals to 100% and 0 to 0%
cb_max_start=1						#taking average for bottom panel from this timestep until the end	#default: 1
						



import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import pylab




def get_shape(filename):
    f = open(filename, 'r')
    for line in f:
        if "Nsize" in line:
            a = line.split()
            row = int(a[1])
            column = int(a[2])
    shape = (row, column)
    lines = row*column

    return shape, lines




def load_topo(filename, plot_data, shape):
    try:
        data = np.loadtxt(filename, unpack=True)
    except:
        print 'Bad topomap data given'

    xValue = data[0]
    yValue = data[1]

    density = data[5]
    if plot_data == "Topological charge":
        density = data[5]
    elif plot_data == "Vorticity in x direction":
        density = data[2]
    elif plot_data == "Vorticity in y direction":
        density = data[3]
    elif plot_data == "Vorticity in z direction":
        density = data[4]
    else:
        print "No plot given for left panel... Plotting Topological charge..."

    return xValue.reshape(shape), yValue.reshape(shape), density.reshape(shape)




def top_panel(filename, plot_data):
    try:
        data = np.loadtxt(filename, unpack=True)
    except:
        print 'Something is wrong with EM.dat'

    top_value = data[10]
    if plot_data == "Total topological charge":
        top_value = data[10]        ###Total topological charge
    elif plot_data == "Torque in eV":
        top_value = data[11]    ###searching
    elif plot_data == "Magnetic field in T":
        top_value = data[18]    ###searching
    else:
        print "No correct top_value found."

    return top_value.tolist()




def get_input(filename):
    f = open(filename, 'r')
    for line in f:
        if "timestep" in line:
            a = line.split()
            timestep = float(a[1])
        elif "duration" in line:
            a = line.split()
            duration = int(a[1])

    return timestep, duration




def load_spins(filename, shape):
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




def degree_to_rad(angle1, angle2):
    rad1 = angle1 * np.pi/180
    rad2 = angle2 * np.pi/180

    return rad1, rad2




def spin_components(theta, phi, shape):
    S_x = np.cos(phi) * np.sin(theta)
    S_y = np.sin(phi) * np.sin(theta)
    S_z = np.cos(theta)

    return S_x.reshape(shape), S_y.reshape(shape), S_z.reshape(shape)




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
steps = number_of_files/(fps*length)
img=1



#search for maximum value (left bottom panel)
cb_max = []
j = cb_max_start

while j < duration:
    x, y, density = load_topo("topomap"+str(j)+".dat", left_panel_plot, shape)
    cb_max.append(np.max(density))
    print str(j)+" of "+str(duration)+" files checked."
    j = j + steps   #optional

scale = np.mean(cb_max)




#top value
top_value = top_panel(EM, top_panel_plot)
timelist = []
for time in range(duration):
    timelist.append((time+1)*timestep)




#loop for all data
n = movie_start
while n < number_of_files:


        #create figure
        fig = pylab.figure()
        ax1 = plt.axes([0.1, 0.1, 0.425, 0.60])
        ax2 = plt.axes([0.55, 0.1, 0.425, 0.60])
        ax3 = plt.axes([0.1, 0.8, 0.85, 0.15])




        #bottom values
        x_array, y_array, density_array = load_topo("topomap"+str(n)+".dat", left_panel_plot, shape)
        theta_deg, phi_deg, X, Y = load_spins('Spinse_'+str(n)+'.dat', shape)
        theta, phi = degree_to_rad(theta_deg, phi_deg)
        S_x, S_y, S_z = spin_components(theta, phi, shape)

        if size != 1:

            x_array, y_array, density_array, X, Y, S_x, S_y, S_z = resize(x_array, y_array, density_array, X, Y, S_x, S_y, S_z, size, shape)




        #get colorbar
        pc = plt.pcolor(S_z, cmap='jet', vmin=-180, vmax=0)
        cm = fig.colorbar(pc, ax=ax2)
        cm.set_ticks([0,-90,-180])
        plt.cla()




        #top panel
        ax3.set_xlim(0, duration*timestep)
        ax3.set_xlabel("time in fs")
        ax3.set_title(top_panel_plot, fontsize=12)
        ax3.set_ylim(min(top_value)-1, max(top_value)+1)
	if min(top_value) == 0 and max(top_value) == 0:
	    ax3.set_yticks([0])
	else:
	    ax3.set_yticks([min(top_value), 0, max(top_value)])

        plt.plot(timelist, top_value, 'r', linewidth=1.5)
        plt.plot(timelist[n-1], top_value[n-1], 'ko', markersize=4)




	#left bottom panel
	ax1.set_xlim([np.min(x_array), np.max(x_array)])
	ax1.set_ylim([np.min(y_array), np.max(y_array)])
	ax1.set_xlabel('x-position in nm')
	ax1.set_ylabel('y-position in nm')
	ax1.set_title(left_panel_plot, fontsize=12)

	levels = MaxNLocator(nbins=100).tick_values(-scale, scale)
	cmap = plt.get_cmap('seismic')
	contour = ax1.contourf(x_array, y_array, density_array, levels=levels, cmap=cmap)
	
	cb = fig.colorbar(contour, ax=ax1)
	cb.set_ticks([0])




	#right bottom panel
	ax2.set_xlim([np.min(X), np.max(X)])
	ax2.set_ylim([np.min(Y), np.max(Y)])
	ax2.set_xlabel("x-position in nm")
	ax2.set_title("Magnetization (M-z)", fontsize=12)

	for sp1 in range(len(X)):
	    for sp2 in range(len(X)):

	        r = 0.5*np.sqrt(S_x[sp1][sp2]**2+S_y[sp1][sp2]**2)
	        if S_z[sp1][sp2] >= 0:
        	    ax2.arrow(X[sp1][sp2], Y[sp1][sp2], r*S_x[sp1][sp2], r*S_y[sp1][sp2], length_includes_head=True, head_width=r*2/3, head_length=r ,color=(S_z[sp1][sp2], np.sqrt(S_x[sp1][sp2]**2+S_y[sp1][sp2]**2), 0))
	        else:
	            ax2.arrow(X[sp1][sp2], Y[sp1][sp2], r*S_x[sp1][sp2], r*S_y[sp1][sp2], length_includes_head=True, head_width=r*2/3, head_length=r ,color=(0, np.sqrt(S_x[sp1][sp2]**2+S_y[sp1][sp2]**2), -S_z[sp1][sp2]))




        #saving
        fig.savefig('spin_'+str(img)+'.png', dpi= 100)          #dpi means resolution of images
        plt.close('all')
        print ('topomap'+str(n)+'.dat was converted.')

        #counter
        img = img+1
        n = n + steps
