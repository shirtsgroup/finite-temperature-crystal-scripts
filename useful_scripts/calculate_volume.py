#!/usr/bin/python
#
# Calculate the volume of a .gro file
#
# Copyright Michael R. Shirts, University of Virginia, 2014
#

import numpy # numerical array library
import pymbar # multistate Bennett acceptance ratio
from pymbar import timeseries # timeseries analysis
from optparse import OptionParser # for parsing command-line options
import pdb

parser = OptionParser()
parser.add_option('-g', '--gro', dest = 'grofile', help = 'Gromacs File') #.gro file to be resized
(options, args) = parser.parse_args()
fname = options.grofile

#=============================================================================================
# CALCULATE THE GRO VOLUME
#=============================================================================================
# Read in input files
fname = options.grofile
infile = open(fname, 'r')
lines = filter(None, (line.rstrip() for line in infile))
infile.close()
print "loading " + fname

#Read in the crystal basis matrix (the last line of the .gro file)
crystal_basis = numpy.zeros([3,3],float)        #Matrix to take a crystal vector into xyz coordinates
xyz_to_crystal = numpy.zeros([3,3],float)       #Matrix to take an xyz vector into crystal coordinates
tokens = lines[len(lines)-1].split()
oldvect=[]
for i,token in enumerate(tokens):
    if i == 0:
	crystal_basis[0,0]=float(token)
	oldvect.append(token)
    elif i==1:
	crystal_basis[1,1]=float(token)
	oldvect.append(token)
    elif i==2:
	crystal_basis[2,2]=float(token)
	oldvect.append(token)
    elif i==3:
	crystal_basis[0,1]=float(token)
	oldvect.append(token)
    elif i==4:
	crystal_basis[0,2]=float(token)
	oldvect.append(token)
    elif i==5:
	crystal_basis[1,0]=float(token)
	oldvect.append(token)
    elif i==6:
	crystal_basis[1,2]=float(token)
	oldvect.append(token)
    elif i==7:
	crystal_basis[2,0]=float(token)
	oldvect.append(token)
    elif i==8:
	crystal_basis[2,1]=float(token)
	oldvect.append(token)
xyz_to_crystal = numpy.linalg.inv(crystal_basis)

#Calculate the initial volume of the gro file
Volume = float(numpy.linalg.det(crystal_basis))

print "Volume: " + str(Volume)

