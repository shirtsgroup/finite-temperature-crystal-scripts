#!/usr/bin/python
#
# =============================================================================================
#=======================HAmiltonian Reweighting VISualizer Toolkit (HARVIST)===================
# =============================================================================================
#
#
#==============================================================================================
#                                       DESCRIPTION
# =============================================================================================
# This toolkit allows one to visualize the overlap between two hamiltonians used for hamiltonian
# reweighting. The visialization can be done with either Histograms or Scatterplots, and the
# visualization can be over a single snapshot or animated over a set of hamiltonians. The visualizer uses a matrix
# Of energy values from samples collected and re-evaluated in each of 2 states.


#
# =============================================================================================
# COPYRIGHT NOTICE
#
# Written by Eric Dybeck <ecd4bd@virginia.edu>.
#
# Copyright (c) 2015 The University of Virginia. All Rights Reserved.
#
# This program is free software; you can redistribute it and/or modify it under the terms of
# the GNU General Public License as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301, USA.
# =============================================================================================

from __future__ import print_function
import numpy
import numpy.random
import scipy
import scipy.optimize
import scipy.stats
from scipy.stats import norm
import matplotlib # for making plots, version 'matplotlib-1.1.0-1'; errors may pop up when using earlier versions
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.font_manager import FontProperties as FP
import matplotlib.animation as ani
import os
import pdb

"""
#==============================================================================================
#                                     MAIN (For Testing...)
#==============================================================================================
#Create the figure and subplot objects
#         subplots(rows, columns)
f,a = plt.subplots(1,1)
potentials=['Amoeba','DesignedA']
scat1, = a.plot([1,2,3,4,5],[1,2,3,4,5],linestyle='none',marker='.',alpha=0.3,color='b',label=potentials[0])
scat2, = a.plot([10,11,12,13,14],[10,11,12,13,14],linestyle='none',marker='.',alpha=0.3,color='g',label=potentials[1])
#time = numpy.arange(5)+1
#aniU = ani.FuncAnimation(fig=f, func=UpdateTime, frames=time, interval=250, blit=True)#, init_func=cleanup)#fargs=(a, scat1, scat2, u_klnT, Terms_l, k1t, k2t, Terms1, Terms2, N_k, xlabel, ylabel, potentials), blit=True)#, init_func=cleanup)
#Load in the raw data
u_kln_save=numpy.load('Ukln.npy')
N_k_save=numpy.load('Nk.npy')
Independent=4
#Determine the spacing for the animation
N_states=len(N_k_save)
lam_range = numpy.arange(1,2*N_states+1)


"""
#==============================================================================================
#                                       MAKE FRAME
#==============================================================================================

#This function takes in a large u_klnT matrix of all the energy terms and returns a [2x2xN] matrix
# with a single energy term for 2 of the states of interest
#
# @Param u_klnT - The large matrix of all energy terms from all states in all other states
# @Param N_k	- The vector of number of samples collected from each state k
# @Param Terms_l- A vector of term strings for each state 'l' e.g. ['Bonds', 'Angles', 'LJ-14'...]
# @Param k1 	- The first state of interest
# @param k2 	- The second state of interest
# @Param Terms_1- A vector of term strings to include from state 1
# @Param Terms_2- A vector of term strings to include from state 2
#
# @Return u_kln - The [2,2,N] Matrix of combined energy terms in each of the 2 states of interest

def MakeFrame(u_klnT, N_k, Terms_l, k1, k2, Terms_1, Terms_2):
    Nk1 = N_k[k1]
    Nk2 = N_k[k2]
    u_kln = numpy.zeros([2,2,numpy.max(Nk1,Nk2)], float)

    #Create the binary vectors of terms. Each index 'i' will be 1 or 0 depending on if the 'i' energy
    #term in the term vector is included in this frame
    binary_term_vect1 = numpy.zeros(len(Terms_l[k1,:]),int)
    binary_term_vect2 = numpy.zeros(len(Terms_l[k2,:]),int)
    for j,term in enumerate(Terms_l[k1,:]):
	binary_term_vect1[j] = term in Terms_1
    for j,term in enumerate(Terms_l[k2,:]):
	binary_term_vect2[j] = term in Terms_2 
    
    #Fill in the energy terms
    u_kln[0,0,:Nk1] = numpy.dot(u_klnT[k1,k1,:Nk1,:], binary_term_vect1)
    u_kln[1,0,:Nk2] = numpy.dot(u_klnT[k2,k1,:Nk2,:], binary_term_vect1)
    u_kln[0,1,:Nk1] = numpy.dot(u_klnT[k1,k2,:Nk1,:], binary_term_vect2)
    u_kln[1,1,:Nk2] = numpy.dot(u_klnT[k2,k2,:Nk2,:], binary_term_vect2)
    
    return u_kln


#==============================================================================================
#                                     CREATE TERMS
#==============================================================================================
# This function takes in an xvg file (or flat file) of energy terms from a trajectory and returns a list of 
# strings corresponding with each term in the file
#
# @Param energyfile	- the energy file being parsed 
#
# @Return Terms 	- a list of terms contained within this energy file

def CreateTerms(energyfile):

    infile = open(energyfile,'r')
    lines = infile.readlines()
    infile.close()

    Terms=[]

    #If this is a gromacs file, loop through the lines and store everything with the preheader 'legend'
    tokens = lines[0].split()
    if tokens[0] == '#' and tokens[1] == 'This' and tokens[2] == 'file' and tokens[3] == 'was' and tokens[4] == 'created':
        for line in lines:
	    tokens=line.split()
	    
	    if tokens[0] == '@' and tokens[2] == 'legend':
		#Pull the term from the energy file
		Terms.append(' '.join(tokens[3:]))
	    if tokens[0] != '#' and tokens[0] != '@' and tokens[0] != '@TYPE':
		break
	#Remove quotation marks from each term
	Terms = [term.translate(None, '"') for term in Terms]
    
    #If this is not a gromacs file, loop through the first line of the energy file
    else:
        for i,token in enumerate(tokens):
	    if i==0:
	        continue
	    Terms.append(token)

    return Terms

#==============================================================================================
#                                     GRAB TERMS
#==============================================================================================

# This function takes in an xvg file of energy terms from a trajectory and returns a matrix of
# all the energy terms for each configuration (optional: only the supplied terms will be grabbed)
#
# @Param energyfile    	- the energy file being parsed
# @Param Terms		- (optional) an array of terms to be grabbed
# @Param ignoreframes   - (optional) a number of frames to ignore from the beginning of the energy file
#
# @Return u_nT		- a matrix of energy terms from the supplied energy file

def GrabTerms(energyfile,Terms='All',ignoreframes=0):

    All_Terms = CreateTerms(energyfile)

    #If no terms have been specified, grab all terms from the energy file
    if Terms == 'All':
	Terms=All_Terms

    #Ensure that all specified terms are in the   

    Term_names = [j for j in Terms if j in All_Terms] 
    Term_indicies = [All_Terms.index(j) for j in Terms if j in All_Terms] 

    infile = open(energyfile,'r')
    lines = infile.readlines()
    infile.close()
    print("loading terms from " + energyfile)

    u_nT = numpy.zeros([50000,len(Term_names)])
    ignore_symbols = ['#', '@', '@TYPE', 'STEP', '#Time', 'TIME(ps)', '=====================']; #Lines to ignore when reading in energies
    N=0
    frames=0
    for line in lines:
        #Skip header lines
        Tokens=line.split()
        if Tokens[0] in ignore_symbols:
	    continue

	frames+=1

        if frames <= ignoreframes:
            continue

        u_nT[N,:] = [float(Tokens[j+1]) for j in Term_indicies] #+1 ignores the time in the energy file
        N+=1
	#print N
    
    return [u_nT[:N,:],Term_names]

#==============================================================================================
#                                 DETERMINE LENGTH
#==============================================================================================

# This function takes in an xvg file of energy terms from a trajectory and returns the number
# of frames in the trajectory
#
# @Param energyfile	- the energy file being parsed
# @Param ignoreframes	- (optional) a number of frames to ignore from the beginning of the energy file
#
# @Return N          	- the number of frames in the trajectory

def DetermineLength(energyfile, ignoreframes=0):

    infile = open(energyfile,'r')
    lines = infile.readlines()
    infile.close()
    
    ignore_symbols = ['#', '@', '@TYPE', 'STEP', 'TIME(ps)', '=====================']; #Lines to ignore when reading in energies
    N=0
    frames=0
    for line in lines:
        #Skip header lines
        Tokens=line.split()
        if Tokens[0] in ignore_symbols:
            continue

	frames+=1

	if frames <= ignoreframes:
	    continue

        N+=1

    return N

#==============================================================================================
#                                   MAKE FRAME
#==============================================================================================

# This function takes in a complete U_klnT matrix and returns a [2,2,N] U_kln matrix representing
# the data that will be used to make a histogram or scatterplot.
#
# @Param u_klnT     	- the complete matrix of all energy data. 'k' is the sampled state, 'l' is 
#			  the evaluated state, 'n' is the sample number,  and T is the energy term
# @Param Terms_l	- the list of terms in the u_klnT matrix, in order
# @Param k1		- the first state in the frame
# @Param k2		- the second state in the frame
# @Param Terms1		- the list of terms to combine from the first state
# @Param Terms2		- the list of terms to combine from the second state
# @Param Nk1		- the total number of frames collected from state 1
# @Param Nk2		- the total number of frames collected from state 2
#
# @Return u_kln         - the matrix of energy data for the two states

def MakeFrame(u_klnT, Terms_l, k1, k2, Terms1, Terms2, Nk1, Nk2):

    u_kln=numpy.zeros([2,2,numpy.max([Nk1,Nk2])],float)

    #Use predefined terms if desired

    #Determine which terms in Terms_l[k] are included in Termsk i.e. [0 1 1 0 1]
    Terms_included1=numpy.zeros([len(u_klnT[0,0,0,:]),1],int)
    Terms_included2=numpy.zeros([len(u_klnT[0,0,0,:]),1],int)
    for i,term in enumerate(Terms_l[k1]):
	if term in Terms1:
	    Terms_included1[i]=1
    for i,term in enumerate(Terms_l[k2]):
        if term in Terms2:
            Terms_included2[i]=1

    #Fill in the new u_kln matrix by summing the terms from above
    u_kln[0,0,:Nk1] = numpy.dot(u_klnT[k1, k1, :Nk1, :], Terms_included1)[:,0]
    u_kln[1,0,:Nk2] = numpy.dot(u_klnT[k2, k1, :Nk2, :], Terms_included1)[:,0]
    u_kln[0,1,:Nk1] = numpy.dot(u_klnT[k1, k2, :Nk1, :], Terms_included2)[:,0]
    u_kln[1,1,:Nk2] = numpy.dot(u_klnT[k2, k2, :Nk2, :], Terms_included2)[:,0]

    return u_kln

#==============================================================================================
#                                   PLOT HISTOGRAM
#==============================================================================================

# This function takes in a u_kln matrix of 2 states and plots the energy difference as histograms
#
# @Param a              - a matplotlib axis handle
# @Param u_kln          - the matrix of energy data for the two states
# @Param N_k            - the number of configurations collected from state k
# @Param xlabel         - the xlabel of the plot
# @Param ylabel         - the ylabel of the plot
# @Param potentials     - the name of the two potentials
# @Param scat1          - a matplotlib scatterplot handle
# @Param scat2          - a matplotlib scatterplot handle
# @Param rescale        - determine whether the histogram axes need to be rescaled
#
# @Return a             - the updated axis handle with the new data

def PlotHistogram(a, u_kln, N_k, xlabel, ylabel, potentials, scat1='None', scat2='None', rescale='yes'):


    a.cla()

    #Create the scatterplot objects (if necessary)
    if scat1=='None' and scat2=='None':
        scat1, = a.plot([],[],linestyle='none',marker='.',alpha=0.3,color='b',label=potentials[0])
        scat2, = a.plot([],[],linestyle='none',marker='.',alpha=0.3,color='g',label=potentials[1])

    #Set the x and y axis labels
    a.set_xlabel(xlabel,fontsize=20)
    a.set_ylabel(ylabel,fontsize=20)

    #Create an array of energy differences sampled in each state
    Energy_Differences1 = (u_kln[0,1,:N_k[0]-1] - u_kln[0,0,:N_k[0]-1])
    Energy_Differences2 = (u_kln[1,1,:N_k[1]-1] - u_kln[1,0,:N_k[1]-1]) 	

    # Fit a normal distribution to the data:
    mu_1, std_1 = norm.fit(Energy_Differences1)
    mu_2, std_2 = norm.fit(Energy_Differences2)

    #Determine the domain of the graph (if it has not already been defined)
    if rescale=='yes':
        if (mu_2-mu_1 > 0):
            xlength=10*(mu_2-mu_1)
        else:
            xlength=10*(mu_1-mu_2)
    else:
        xlength = float(a.get_xlim()[1]) - float(a.get_xlim()[0])

    xmin=0.5*(mu_2+mu_1-xlength)
    xmax=0.5*(mu_2+mu_1+xlength)
    a.set_xlim([xmin,xmax])
    #a.set_ylim([0,0.25])
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
    a.fill_between(x,0,dU1,alpha=0.3,edgecolor='k',facecolor=scat1.get_color(),label=potentials[0])
    a.fill_between(x,0,dU2,alpha=0.3,edgecolor='k',facecolor=scat2.get_color(),label=potentials[1])

    ##Create histograms of the data
    #a.hist(u_kln[0,1,:N_k[0]]-u_kln[0,0,:N_k[0]], alpha=0.3, color='b', label=potentials[0])
    #a.hist(u_kln[1,1,:N_k[1]]-u_kln[1,0,:N_k[1]], alpha=0.3, color='g', label=potentials[1])

    ##Set the x and y axis labels
    #a.set_xlabel(xlabel,fontsize=20)
    #a.set_ylabel(ylabel,fontsize=20)


    return a, scat1, scat2



#==============================================================================================
#                                   PLOT FRAME
#==============================================================================================

# This function takes in a u_kln matrix of 2 states and plots the data as either a set of 
# histograms or as a scatterplot
#
# @Param a		- a matplotlib axis handle
# @Param u_kln		- the matrix of energy data for the two states
# @Param N_k		- the number of configurations collected from state k
# @Param xlabel		- the xlabel of the plot
# @Param ylabel		- the ylabel of the plot
# @Param potentials	- the name of the two potentials
# @Param scat1          - a matplotlib scatterplot handle
# @Param scat2          - a matplotlib scatterplot handle
# @Param rescale	- determine whether the scatterplot axes need to be rescaled
#
# @Return a,scat1,scat2	- the updated axis handle with the new data as well as the individual scatterplots

def PlotFrame(a, u_kln, N_k, xlabel, ylabel, potentials, scat1='None', scat2='None', rescale='yes'):


    #Create the scatterplot objects (if necessary)
    if scat1=='None' and scat2=='None':
        scat1, = a.plot([],[],linestyle='none',marker='.',alpha=0.3,color='b',label=potentials[0])
        scat2, = a.plot([],[],linestyle='none',marker='.',alpha=0.3,color='g',label=potentials[1])

    #Clear the old data
    scat1.set_data([], [])
    scat2.set_data([], [])

    #Fill in the data from the u_kln matrix
    scat1.set_data(u_kln[0,1,:N_k[0]],u_kln[0,0,:N_k[0]])
    scat2.set_data(u_kln[1,1,:N_k[1]],u_kln[1,0,:N_k[1]])

    #Set the x and y axis labels
    a.set_xlabel(xlabel,fontsize=20)
    a.set_ylabel(ylabel,fontsize=20)

    #Set the x and y axis limits by fitting both clouds to a line with a slope of 1 and setting the domain to the interquartile range of the largest cloud
    b1 = -1.0*numpy.average(u_kln[0,1,:N_k[0]]-u_kln[0,0,:N_k[0]])
    b2 = -1.0*numpy.average(u_kln[1,1,:N_k[1]]-u_kln[1,0,:N_k[1]])
    b=0.5*(b1+b2)
    q75_1_x, q25_1_x = numpy.percentile(u_kln[0,1,:N_k[0]], [75, 25])
    q75_2_x, q25_2_x = numpy.percentile(u_kln[1,1,:N_k[1]], [75, 25])
    q75_1_y, q25_1_y = numpy.percentile(u_kln[0,0,:N_k[0]], [75, 25])
    q75_2_y, q25_2_y = numpy.percentile(u_kln[1,0,:N_k[1]], [75, 25])
    IQR1_x = q75_1_x-q25_1_x
    IQR2_x = q75_2_x-q25_2_x
    IQR1_y = q75_1_y-q25_1_y
    IQR2_y = q75_2_y-q25_2_y
    IQR=numpy.amax([IQR1_x, IQR2_x, IQR1_y, IQR2_y])
    xmin=q25_1_x-9*IQR
    xmax=q75_1_x+9*IQR
    #if(IQR1>IQR2):
    #    xmin=q25_1-7*IQR1
    #    xmax=q75_1+7*IQR1
    #else:
    #    xmin=q25_2-7*IQR2
    #    xmax=q75_2+7*IQR2
    ymin=xmin+b
    ymax=xmax+b


    #if rescale=='yes':
    #    Domain1 = numpy.max(u_kln[0,1,:N_k[0]]) - numpy.min(u_kln[0,1,:N_k[0]])
    #    Domain2 = numpy.max(u_kln[1,1,:N_k[1]]) - numpy.min(u_kln[1,1,:N_k[1]])
    #    Range1 = numpy.max(u_kln[0,0,:N_k[0]]) - numpy.min(u_kln[0,0,:N_k[0]])
    #    Range2 = numpy.max(u_kln[1,0,:N_k[0]]) - numpy.min(u_kln[1,0,:N_k[0]])
    #    Length=numpy.max([Domain1, Domain2, Range1, Range2])
    #else:
    #    xmin,xmax=a.get_xlim()
    #    Length=xmax-xmin - (numpy.max(u_kln[0,1,:N_k[0]]) - numpy.min(u_kln[0,1,:N_k[0]]))
    #xmin=numpy.min(u_kln[0,1,:N_k[0]])-0.5*Length
    #xmax=numpy.max(u_kln[0,1,:N_k[0]])+0.5*Length
    #ymin=xmin+b
    #ymax=xmax+b

    ##Set the x and y axis limits by fitting the first cloud to a line with a slope of 1 and setting the domain to the max range of the cloud
    #b=-1.0*numpy.average(u_kln[0,1,:N_k[0]]-u_kln[0,0,:N_k[0]])
    #Domain1 = numpy.max(u_kln[0,1,:N_k[0]]) - numpy.min(u_kln[0,1,:N_k[0]])
    #Range1 = numpy.max(u_kln[0,0,:N_k[0]]) - numpy.min(u_kln[0,0,:N_k[0]])
    #Length=numpy.max([Domain1, Range1])

    #xmin=numpy.min(u_kln[0,1,:N_k[0]])-1.5*Length
    #xmax=numpy.max(u_kln[0,1,:N_k[0]])+1.5*Length
    #ymin=xmin+b
    #ymax=xmax+b

    a.set_xlim([xmin,xmax])
    a.set_ylim([ymin,ymax])

    #plt.show()

    #Return the updated axis handle
    return a, scat1, scat2

    #Text update, note that we are using the non-formatted text string from before, we now just set the formatting
    #atxt.set_text(text_template % (lamR, lamA, lamLJ, lamE))
    #If we have a 2d array we are updating, we would call set_array(...) on the object and fill in "..." with the array object.
    #Return all the objects we animate
    #plt.show()
    #plt.savefig(filename)
    #return a, scat1, scat2, atxt

#==============================================================================================
#                                  ANIMATE FRAMES
#==============================================================================================

#This function takes in a large u_klnT matrix and animates a series of specified frames
#
# @Param f		- a matplotlib figure handle
# @Param a              - a matplotlib axis handle
# @Param u_klnT         - the complete matrix of all energy data. 'k' is the sampled state, 'l' is
#                         the evaluated state, 'n' is the sample number,  and T is the energy term
# @Param Terms_l        - the list of terms in the u_klnT matrix, in order
# @Param k1t            - the vector indicating which state will be used as k1 at time t
# @Param k2t            - the vector indicating which state will be used as k2 at time t
# @Param Terms1         - the list of terms to combine from the first state
# @Param Terms2         - the list of terms to combine from the second state
# @Param N_k            - the number of configurations collected from state k
# @Param xlabel         - the xlabel of the plot
# @Param ylabel         - the ylabel of the plot
# @Param potentials     - the name of the two potentials
# @Param type		- the type of figure to animate. Either 'scat' or 'hist'.
#
# @Return a             - the updated axis handle with the new data

def AnimateFrames(f, a, u_klnT, Terms_l, k1t, k2t, Terms1, Terms2, N_k, xlabel='Xaxis', ylabel='Yaxis', potentials=['None', 'None'],type='scat'):
    time=numpy.arange(len(k1t))

    #f2,a2 = plt.subplots(1,1)

    #Create the scatterplot objects
    scat1, = a.plot([],[],linestyle='none',marker='.',alpha=0.3,color='b',label=potentials[0])
    scat2, = a.plot([],[],linestyle='none',marker='.',alpha=0.3,color='g',label=potentials[1])
    #atxt = a.text(1.6,3.75, '', fontsize=20) #Position text box

    aniU = ani.FuncAnimation(fig=f, func=UpdateTimeScat, frames=time, interval=250, fargs=(a, scat1, scat2, u_klnT, Terms_l, k1t, k2t, Terms1, Terms2, N_k, xlabel, ylabel, potentials,type), blit=True)#, init_func=cleanup)
    #aniU.save('Test.mp4', dpi=100)
    return aniU

def UpdateTimeScat(t, a, scat1, scat2, u_klnT, Terms_l, k1t, k2t, Terms1, Terms2, N_k,xlabel,ylabel,potentials,type):
    u_klnt = MakeFrame(u_klnT, Terms_l, k1t[t], k2t[t], Terms1, Terms2, N_k[k1t[t]], N_k[k1t[t]])
    N_kt = [N_k[k1t[t]], N_k[k2t[t]]]
    resize='yes'
    if t>0:
        resize='no'
    if type=='scat':
        a=PlotFrame(a, u_klnt, N_kt, xlabel, ylabel, potentials,scat1=scat1,scat2=scat2,rescale=resize)
    elif type=='hist':
	a=PlotHistogram(a, u_klnt, N_kt, xlabel, ylabel, potentials,scat1=scat1,scat2=scat2,rescale=resize)
    #ChargeSeparation=[0.0550, 0.0700, 0.0850, 0.1000, 0.1150, 0.1300, 0.1450, 0.1600, 0.1750, 0.1900];
    #atxt.set_text('Charge Separation =' + str(ChargeSeparation[t]))

    
    return a#,atxt





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
    #Return all of the objects being animated, order does not mater
    #return hist1, bins1, patches1, hist2, bins2, patches2, atxt
    return scat1, scat2


