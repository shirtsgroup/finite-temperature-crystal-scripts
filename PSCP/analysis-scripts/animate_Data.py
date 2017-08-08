#!/usr/bin/python
#
# Animate a set of data
#
# Copyright Eric Dybeck, University of Virginia, 2015
#
import numpy # numerical array library
import pymbar # multistate Bennett acceptance ratio
from pymbar import timeseries # timeseries analysis
from optparse import OptionParser # for parsing command-line options
from scipy.stats import norm
import os
import pdb
import matplotlib # for making plots, version 'matplotlib-1.1.0-1'; errors may pop up when using earlier versions
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.font_manager import FontProperties as FP
import matplotlib.animation as ani

######################################################################################################
#################################### ANIMATION SETUP #################################################
######################################################################################################

#Create the figure and subplot objects
#         subplots(rows, columns)
f,a = plt.subplots(1,1)
'''
Anything you want to animate, you should define as an empty object first, as is shown below.
Note how the x and y data is input as an empty list, but everything else about the line is defined. We do this since the 'animation' will be simply updating the data.
Notes:
-You MUST assign the objects to a variable with a comma after the variable name. If you don't include the comma, it assigns a tuple with 1 object to the variable, and manipulating that requires something like "var[0]" every time. Adding the comma at assignment is easier.
-Animating scatter plots: Do NOT use '.scatter(...)', this is a pain to animate. instead use '.plot([], [], linestyle='', markersize=x, color=y, marker=z, markeredgewidth=t)' where x,y,z,t are your choice, just make linestyle an empty string to not plot the line.
-Animating 2d surfaces: Dont set empty data, set some baseline data first. Then use .set_array(...), I'll show an example of this later.
'''
#Define Empty Templates
#EACH LINE BEING ANIMATED
#att, = a.plot([], [], linestyle='solid', color='#FFA500', linewidth=3, label='Attractive Basis') #R
#rep, = a.plot([], [], '-r', linewidth=3, label='Repulsive Basis') #A
#ele, = a.plot([], [], '-b', linewidth=3, label='Electrostatic Basis') #E
#pot, = a.plot([], [], '--g', linewidth=7, label='Effective Potential') #Full

Potentials=["AMOEBA","Point-Charge"]

#hist1,bins1,patches1= a.hist([],20,alpha=0.3)
#hist2,bins2,patches2= a.hist([],20,alpha=0.3)
scat1, = a.plot([1,2,3,4,5],[1,2,3,4,5],linestyle='none',marker='.',alpha=0.3,color='b',label=Potentials[0])
scat2, = a.plot([10,11,12,13,14],[10,11,12,13,14],linestyle='none',marker='.',alpha=0.3,color='g',label=Potentials[1])
#plt.show()
#Load in the raw data
u_kln_save=numpy.load('Ukln.npy')
N_k_save=numpy.load('Nk.npy')
Independent=4
#Determine the spacing for the animation
N_states=len(N_k_save)
lam_range = numpy.arange(1,2*N_states+1)


#Blank text box
#This is how you animate text, you first define a string with formatting things you can feed into, but then dont set the formatters.
#text_template = r'$\lambda_{R\,} = %.2f$' + '\n' + r'$\lambda_{A\,} = %.2f$' + '\n' + r'$\lambda_{LJ} = %.2f$' + '\n' + r'$\lambda_{E\,} = %.2f$'
atxt = a.text(1.6,3.75, '', fontsize=20) #Position text box

#Constants of graph, things NOT being animated
#a.axvline(linewidth=2,color='k', ymin=ymin, ymax=ymax)
#a.axhline(linewidth=2,color='k', xmin=xmin, xmax=xmax)
#a.set_xlabel(r'Radial Distance in $\sigma$', fontsize=20)
#a.set_ylabel(r'Potential in $\epsilon$', fontsize=20)
#legend
#plt.show()

######################################################################################################
################################## ANIMATION EXECUTION ###############################################
######################################################################################################

#Define function to reset the curves.
#This function blanks the curves and helps create a smoother animation. Its also a good way for you to keep track of what your animation function should return and keep track of how each functions data changes
def cleanup():
    #Returns blank copies of everyhing being animated
    #Note the format used for each data type:
    #.plot() uses .set_data([], [])
    #2d-surfaces will use .set_array([])
    #text uses .set_text('')

    #This is the list of all objects being animated
    #hist1.set_data([], [])
    #hist2.set_data([], [])
    #a.cla()
    scat1.set_data([], [])
    scat2.set_data([], [])
    atxt.set_text('')
    #Return all of the objects being animated, order does not mater
    #return hist1, bins1, patches1, hist2, bins2, patches2, atxt
    return scat1, scat2, atxt

#animation function, returns all functions being animated
#This is the workhorse function
def movegraph(k): #Function accepts a single input, although the FuncAnimation can pass in more options

    #Update the objects being animated.
    #Use .set_data(x, y) to create new data

    filename='frame' + str(k+6) + '.jpg'
    if k<1:
	k=1
    elif k>=N_states:
	k=N_states-1

    ##########USE THIS FOR CREATING HISTOGRAMS OF THE ENERGY DIFFERENCE#############################
    """
    a.cla()
 
    Energy_Differences1 = (u_kln_save[0,k,:N_k_save[0]-1] - u_kln_save[0,0,:N_k_save[0]-1])/float(Independent)
    Energy_Differences2 = (u_kln_save[k,k,:N_k_save[k]-1] - u_kln_save[k,0,:N_k_save[k]-1])/float(Independent)
    
    # Fit a normal distribution to the data:
    mu_1, std_1 = norm.fit(Energy_Differences1)
    mu_2, std_2 = norm.fit(Energy_Differences2)
    #Determine the domain of the graph (if it has not already been defined)
    if k==1:
	if (mu_2-mu_1 > 0):
	    xlength=5*(mu_2-mu_1)
	else:
	    xlength=5*(mu_1-mu_2)
    else:
	xlength = float(a.get_xlim()[1]) - float(a.get_xlim()[0])

    xmin=0.5*(mu_2+mu_1-xlength)
    xmax=0.5*(mu_2+mu_1+xlength)
    a.set_xlim([xmin,xmax])
    a.set_ylim([0,3.0])
    x=numpy.linspace(xmin, xmax, 1000)
    dU1=norm.pdf(x, mu_1, std_1)
    dU2=norm.pdf(x, mu_2, std_2)	
    
    scat1.set_data(x,dU1)
    scat2.set_data(x,dU2)
    scat1.set_linestyle('-')
    scat2.set_linestyle('-')
    scat1.set_marker('None')
    scat2.set_marker('None')
    scat1.set_linewidth(3)
    scat2.set_linewidth(3)
    a.fill_between(x,0,dU1,alpha=0.3,edgecolor='k',facecolor=scat1.get_color(),label=Potentials[0])
    a.fill_between(x,0,dU2,alpha=0.3,edgecolor='k',facecolor=scat2.get_color(),label=Potentials[1])
    a.set_xlabel('Energy Difference Between ' + Potentials[1] + ' and ' + Potentials[0],fontsize=20)
    a.set_ylabel('Normalized Probability',fontsize=20)
    #weights = numpy.ones_like(Energy_Differences1)/len(Energy_Differences1)
    #hist1,bins1,patches1 = a.hist(Energy_Differences1,20,alpha=0.3,label=Potentials[0], weights=weights)
    #weights = numpy.ones_like(Energy_Differences2)/len(Energy_Differences2)
    #hist2,bins2,patches2 = a.hist(Energy_Differences2,20,alpha=0.3,label=Potentials[1], weights=weights)
    return a,scat1, scat2, atxt 
    """
    ################################################################################################

    ##########USE THIS FOR CREATING SCATTER PLOTS OF THE CONFIGURATION SPACE########################
    
    #a.cla()
    
    #Update the domain and range by fitting both clouds to a line with a slope of 1 and setting the domain to the interquartile range of the largest cloud
    #A1 = numpy.ones(N_k_save[0])
    #w1 = numpy.linalg.lstsq(A1.T,(u_kln_save[0,k,:N_k_save[0]]-u_kln_save[0,0,:N_k_save[0]]))[0]1
    #A2 = numpy.ones(N_k_save[0])
    #w2 = numpy.linalg.lstsq(A1.T,(u_kln_save[k,k,:N_k_save[k]]-u_kln_save[k,0,:N_k_save[k]]))[0]
    b1 = -1.0*numpy.average(u_kln_save[0,k,:N_k_save[0]]-u_kln_save[0,0,:N_k_save[0]])
    b2 = -1.0*numpy.average(u_kln_save[k,k,:N_k_save[k]]-u_kln_save[k,0,:N_k_save[k]])
    b=0.5*(b1+b2)
    q75_1, q25_1 = numpy.percentile(u_kln_save[0,k,:N_k_save[0]], [75, 25])
    q75_2, q25_2 = numpy.percentile(u_kln_save[k,k,:N_k_save[k]], [75, 25])
    IQR1 = q75_1-q25_1
    IQR2 = q75_2-q25_2
    #xmin1 = numpy.min(u_kln_save[0,k,:N_k_save[0]])
    #xmin2 = numpy.min(u_kln_save[k,k,:N_k_save[k]])
    #xmax1 = numpy.max(u_kln_save[0,k,:N_k_save[0]])
    #xmax2 = numpy.max(u_kln_save[k,k,:N_k_save[k]])
    if(IQR1>IQR2):
	xmin=q25_1-3*IQR1
	xmax=q75_1+3*IQR1
    else:
	xmin=q25_2-3*IQR2
        xmax=q75_2+3*IQR2
    ymin=xmin+b
    ymax=xmax+b
    a.set_xlim([xmin,xmax])
    a.set_ylim([ymin,ymax])
    a.set_xlabel(Potentials[1] + ' Total Energy',fontsize=20)
    a.set_ylabel(Potentials[0] + ' Total Energy',fontsize=20)
    scat1.set_data(u_kln_save[0,k,:N_k_save[0]],u_kln_save[0,0,:N_k_save[0]])
    scat2.set_data(u_kln_save[k,k,:N_k_save[k]],u_kln_save[k,0,:N_k_save[k]])

    #Text update, note that we are using the non-formatted text string from before, we now just set the formatting
    #atxt.set_text(text_template % (lamR, lamA, lamLJ, lamE))
    #If we have a 2d array we are updating, we would call set_array(...) on the object and fill in "..." with the array object.
    #Return all the objects we animate
    #plt.show()
    #plt.savefig(filename)
    return a, scat1, scat2, atxt
    
    ###################################################################################################    

#Here is the animation command in all its glory. Quite simple really
#I use the FuncAnimation since it animates a function, the rest of the API is a bit hard to understand (I haven't figured it out yet)
'''
ani.FuncAnimation(fig, fun, iter, interval=xx, blit=bool, init_func=cleanfun)
fig: The figure object where all the objects you want to iterate exist.
fun: The function we will be using to update the objects in fig we want to update
iter: An iterable object which will be passed into fun to update the figures, there will be a number of frames equal to the number of objects interated over
interval=xx: xx is the number of ms between each frame. Increase to slow down your animation
blit=bool: Helps with the backend animation. Not quite sure what it does, but 2d animations can sometimes have issues if its True
init_func=cleanfun: cleanfun is the function you use to blank the figure, in my case, its called cleanup
'''
#movegraph(1)
aniU = ani.FuncAnimation(f, movegraph, lam_range, interval=250, blit=True, init_func=cleanup)

#Save the animation, note that nothing is rendered until .save or .show is called, so up unto this point, the function should be very fast
#I commmented these out so you can choose how you want the output
aniU.save('Scatters.mp4', dpi=100)
plt.show()

