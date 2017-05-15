#/usr/bin/python
#PYTHON SCRIPT TO ADD THE CORRECT BOX VECTORS TO A GRO FILE GIVEN A TINKER XYZ FILE
#param x - The name of the .gro file
#param n - The scaling factor to apply to each xyz coordinate

import numpy # numerical array library
from optparse import OptionParser # for parsing command-line options
import os
import pdb

#=============================================================================================
# READ IN USER INPUTS
#=============================================================================================
parser = OptionParser()
parser.add_option('-g', dest = 'grofile', help = 'Input gro file', default='frame.gro')
parser.add_option('-x', dest = 'xyzfile', help = 'Input xyz file', default='frame.xyz')
parser.add_option('-s', dest = 'num' , help = 'scaling factor', default='0.1')
(options, args) = parser.parse_args()
scaling_factor = float(options.num)

fname = options.xyzfile

infile = open(fname, 'r')
lines = infile.readlines()
infile.close()

#Read in the second line of the xyz file
tokens = lines[1].split()

A=float(tokens[0])
B=float(tokens[1])
C=float(tokens[2])
Alpha=float(tokens[3])
Beta=float(tokens[4])
Gamma=float(tokens[5])

#Calculate the box vectors
V1x=float(A)
V1y=0.0
V1z=0.0
V2x=B*numpy.cos(Gamma*numpy.pi/180)
V2y=B*numpy.sin(Gamma*numpy.pi/180)
V2z=0.0
V3x=C*numpy.cos(Beta*numpy.pi/180)
V3y=C*numpy.cos(Alpha*numpy.pi/180)*numpy.sin(Gamma*numpy.pi/180)
V3z=float(C**2 - V3x**2 - V3y**2)**0.5

#Apply the scaling factor to go from A to nm
V1x=V1x*scaling_factor
V2x=V2x*scaling_factor
V2y=V2y*scaling_factor
V3x=V3x*scaling_factor
V3y=V3y*scaling_factor
V3z=V3z*scaling_factor

#Write the last line to the file
if Alpha==90.0 and Beta==90.0 and Gamma==90.0:
    newstr=str("%.8f" % V1x) + " " + str("%.8f" % V2y) + " " + str("%.8f" % V3z)
else: 
    newstr=str("%.8f" % V1x) + " " + str("%.8f" % V2y) + " " + str("%.8f" % V3z) + " " + str("%.8f" % V1y) + " " + str("%.8f" % V1z) + " " + str("%.8f" % V2x) + " " + str("%.8f" % V2z) + " " + str("%.8f" % V3x) + " " + str("%.8f" %V3y)
commandstring = "sed -i 's/0.00000000   0.00000000   0.00000000/ " + newstr + "/g' " + options.grofile
os.system(commandstring)

