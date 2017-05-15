#!/usr/bin/python
#
# PYTHON SCRIPT TO MINIMIZE A CRYSTALS BOX VECTORS USING GROMACS

import numpy # numerical array library
from scipy.optimize import minimize
import resize_gro
import calculate_gro_volume
from optparse import OptionParser # for parsing command-line options
import os
import pdb

def xtalmin(fname):
  
    infile = open(fname, 'r')
    lines = filter(None, (line.rstrip() for line in infile))
    infile.close() 
    x0=numpy.array([float(j) for j in lines[len(lines)-1].split()])
    if len(x0)>3:
    	x0[3]=x0[5]
	x0[4]=x0[7]
	x0[5]=x0[8]
    	x0=x0[0:6]
    #x0 = float(calculate_gro_volume.Volume(fname))
    res = minimize(potenergy, x0, method='Nelder-Mead', options={'xtol': 1e-2, 'disp': True})

def potenergy(Vector):

    #Resize the gro file to the new vector
    if len(Vector)==3:
        final_vector = Vector
        vectstr = ' '.join([str(j) for j in final_vector])
    elif len(Vector)>3:
        final_vector = numpy.zeros(9,float)
        final_vector[0:3]=Vector[0:3]
        final_vector[5]=Vector[3]
        final_vector[7]=Vector[4]
        final_vector[8]=Vector[5]
        vectstr = ' '.join([str(j) for j in final_vector])
    
    #resize_gro.changeBoxVector(fname="MIN.gro",volume=Vector)
    resize_gro.changeBoxVector(fname="MIN.gro",boxvect=vectstr)

    #Rerun the minimization subroutine
    os.system('./submit_minimization_local.sh')

    energyfile=open('energy.xvg', 'r')
    lines = filter(None, (line.rstrip() for line in energyfile))
    energyfile.close()

    print float(lines[len(lines)-1].split()[1])
    return float(lines[len(lines)-1].split()[1])
