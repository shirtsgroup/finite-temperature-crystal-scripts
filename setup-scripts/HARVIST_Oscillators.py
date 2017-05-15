#!/bin/python
#This script generates movies based on a set of 1D oscillators

import numpy
import numpy.random
import scipy
import scipy.optimize
import scipy.stats
import Harvist
import matplotlib # for making plots, version 'matplotlib-1.1.0-1'; errors may pop up when using earlier versions
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.font_manager import FontProperties as FP
import matplotlib.animation as ani
from pymbar import testsystems, EXP, EXPGauss, BAR, MBAR
import os
import pdb


###################################################################################
#			OFFSET OSCILLATORS
###################################################################################



#K_k = numpy.ones(12,int) # spring constants for each state
#O_k = numpy.array([0, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5]) # offsets for spring constants

#Offset Oscillators
K=40
K_k = 0.1*numpy.ones(2*K,float)
K_k[K:]=0.1
O_k = numpy.zeros(2*K,float)
O_k[K:] = 2.5

N_k=1000*numpy.ones(2*K, int)
N_max = numpy.max(N_k)

beta=1
randomsample = testsystems.harmonic_oscillators.HarmonicOscillatorsTestCase(O_k=O_k, K_k=K_k, beta=beta)
[x_kn,u_kln,N_k] = randomsample.sample(N_k,mode='u_kln')

u_klnT = numpy.zeros([len(numpy.nonzero(N_k)[0]), len(N_k), N_max, 1]) # u_klnT is the complete matrix of all energy data. 'k' is the sampled state, 'l' is

                                                                        # the evaluated state, 'n' is the sample number,  and T is the energy term
u_klnT_temp = numpy.zeros([len(numpy.nonzero(N_k)[0]), len(N_k), N_max, 1])

Terms_l=[]
for i in range(len(N_k)):
    Terms_l.append(['Potential'])

u_klnT_temp[:,:,:,0]=u_kln
K=len(K_k)/2
for j in range(K):
    u_klnT[0,0,:,:] += u_klnT_temp[j,j,:,:]
    u_klnT[0,1,:,:] += u_klnT_temp[j,j+K,:,:]
    u_klnT[1,1,:,:] += u_klnT_temp[j+K,j+K,:,:]
    u_klnT[1,0,:,:] += u_klnT_temp[j+K,j,:,:]

#u_klnT[0,1,:,0] += 0.15*(x_kn[1,:] - x_kn[2,:])**2 #Tak on a repulsion penalty
#u_klnT[0,1,:,0] += -0.15*(x_kn[1,:] - x_kn[0,:])**2 #Tak on an attraction bonus


#Plot the Oscillators
f,a = plt.subplots(1,1)
xlabel='Test Oscillators Energy'
ylabel='Target Oscillators Energy'
u_kln = Harvist.MakeFrame(u_klnT, Terms_l, 0, 1, ['Potential'], ['Potential'], N_k[0], N_k[1])
a=Harvist.PlotFrame(a, u_kln, numpy.array([N_k[0], N_k[len(N_k)-1]]), 'Sampled Oscillator', 'Target Oscillator', ['Sampled', 'Target'])[0]
a.legend(loc='upper right')
f2,a2 = plt.subplots(1,1)
a2=Harvist.PlotHistogram(a2, u_kln, numpy.array([N_k[0], N_k[len(N_k)-1]]), 'Sampled Oscillator', 'Target Oscillator', ['Sampled', 'Target'])[0]
a2.legend(loc='upper right')

###################################################################################
#                       STIFF OSCILLATORS
###################################################################################

#Stiff Oscillators
K=80
K_k = 0.1*numpy.ones(2*K,float)
K_k[K:] = 0.20
#K_k[K*1.5:]=0.05
O_k = numpy.zeros(2*K,float)
O_k[K:] = 0.5

N_k=1000*numpy.ones(2*K, int)
N_max = numpy.max(N_k)

beta=1
randomsample = testsystems.harmonic_oscillators.HarmonicOscillatorsTestCase(O_k=O_k, K_k=K_k, beta=beta)
[x_kn,u_kln,N_k] = randomsample.sample(N_k,mode='u_kln')

u_klnT = numpy.zeros([len(numpy.nonzero(N_k)[0]), len(N_k), N_max, 1]) # u_klnT is the complete matrix of all energy data. 'k' is the sampled state, 'l' is

                                                                        # the evaluated state, 'n' is the sample number,  and T is the energy term
u_klnT_temp = numpy.zeros([len(numpy.nonzero(N_k)[0]), len(N_k), N_max, 1])

Terms_l=[]
for i in range(len(N_k)):
    Terms_l.append(['Potential'])

u_klnT_temp[:,:,:,0]=u_kln
K=len(K_k)/2
for j in range(K):
    u_klnT[0,0,:,:] += u_klnT_temp[j,j,:,:]
    u_klnT[0,1,:,:] += u_klnT_temp[j,j+K,:,:]
    u_klnT[1,1,:,:] += u_klnT_temp[j+K,j+K,:,:]
    u_klnT[1,0,:,:] += u_klnT_temp[j+K,j,:,:]


#Plot the Oscillators
f,a = plt.subplots(1,1)
xlabel='Test Oscillators Energy'
ylabel='Target Oscillators Energy'
u_kln = Harvist.MakeFrame(u_klnT, Terms_l, 0, 1, ['Potential'], ['Potential'], N_k[0], N_k[1])
a=Harvist.PlotFrame(a, u_kln, numpy.array([N_k[0], N_k[len(N_k)-1]]), 'Sampled Oscillator', 'Target Oscillator', ['Sampled', 'Target'])[0]
a.legend(loc='upper right')
f2,a2 = plt.subplots(1,1)
a2=Harvist.PlotHistogram(a2, u_kln, numpy.array([N_k[0], N_k[len(N_k)-1]]), 'Sampled Oscillator', 'Target Oscillator', ['Sampled', 'Target'])[0]
a2.legend(loc='upper right')


###################################################################################
#                       BOTH OFFSET AND STIFF OSCILLATORS
###################################################################################

#Stiff Oscillators
K=16
K_k = 0.1*numpy.ones(2*K,float)
K_k[K:] = 0.15
O_k = numpy.zeros(2*K,float)
O_k[K:] = 3.5

N_k=1000*numpy.ones(2*K, int)
N_max = numpy.max(N_k)

beta=1
randomsample = testsystems.harmonic_oscillators.HarmonicOscillatorsTestCase(O_k=O_k, K_k=K_k, beta=beta)
[x_kn,u_kln,N_k] = randomsample.sample(N_k,mode='u_kln')

u_klnT = numpy.zeros([len(numpy.nonzero(N_k)[0]), len(N_k), N_max, 1]) # u_klnT is the complete matrix of all energy data. 'k' is the sampled state, 'l' is

                                                                        # the evaluated state, 'n' is the sample number,  and T is the energy term
u_klnT_temp = numpy.zeros([len(numpy.nonzero(N_k)[0]), len(N_k), N_max, 1])

Terms_l=[]
for i in range(len(N_k)):
    Terms_l.append(['Potential'])

u_klnT_temp[:,:,:,0]=u_kln
K=len(K_k)/2
for j in range(K):
    u_klnT[0,0,:,:] += u_klnT_temp[j,j,:,:]
    u_klnT[0,1,:,:] += u_klnT_temp[j,j+K,:,:]
    u_klnT[1,1,:,:] += u_klnT_temp[j+K,j+K,:,:]
    u_klnT[1,0,:,:] += u_klnT_temp[j+K,j,:,:]


#Plot the Oscillators
f,a = plt.subplots(1,1)
xlabel='Test Oscillators Energy'
ylabel='Target Oscillators Energy'
u_kln = Harvist.MakeFrame(u_klnT, Terms_l, 0, 1, ['Potential'], ['Potential'], N_k[0], N_k[1])
a=Harvist.PlotFrame(a, u_kln, numpy.array([N_k[0], N_k[len(N_k)-1]]), 'Sampled Oscillator', 'Target Oscillator', ['Sampled', 'Target'])[0]
a.legend(loc='upper right')
f2,a2 = plt.subplots(1,1)
a2=Harvist.PlotHistogram(a2, u_kln, numpy.array([N_k[0], N_k[len(N_k)-1]]), 'Sampled Oscillator', 'Target Oscillator', ['Sampled', 'Target'])[0]
a2.legend(loc='upper right')






plt.show()
