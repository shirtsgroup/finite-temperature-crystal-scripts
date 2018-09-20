# /usr/bin/python
# PYTHON SCRIPT TO ADD THE CORRECT BOX VECTORS TO A GRO FILE GIVEN A TINKER XYZ FILE
# param x - The name of the .gro file
# param n - The scaling factor to apply to each xyz coordinate

import numpy as np
from optparse import OptionParser
import os

#=============================================================================================
# READ IN USER INPUTS
#=============================================================================================
parser = OptionParser()
parser.add_option('-g', dest='grofile', help='Input Gromacs file (.gro)', default='frame.gro')
parser.add_option('-x', dest='xyzfile', help='Input Tinker file (.xyz)', default='frame.xyz')
#parser.add_option('-s', dest='scaling_factor', help='scaling factor to change the units', default='0.1')
(options, args) = parser.parse_args()
#scaling_factor = float(options.scaling_factor)

# If we add in another program than this may need to be changed somehow
scaling_factor = 0.1

xyzfile = options.xyzfile
grofile = options.grofile

# Reading in the Tinker coordinate file
infile = open(xyzfile, 'r')
lines = infile.readlines()
infile.close()

# Read in the second line of the xyz file
tokens = lines[1].split()

# Tinker lattice vectors
A = float(tokens[0])
B = float(tokens[1])
C = float(tokens[2])

# Tinker lattice angles
Alpha = float(tokens[3])
Beta = float(tokens[4])
Gamma = float(tokens[5])

# Calculate the box vectors in tensor form and applying the scaling factor to convert from Ang. to nm
V1x = float(A) * scaling_factor
V1y = 0.0
V1z = 0.0
V2x = B * np.cos(np.radians(Gamma)) * scaling_factor
V2y = B * np.sin(np.radians(Gamma)) * scaling_factor
V2z = 0.0
V3x = C * np.cos(np.radians(Beta)) * scaling_factor
V3y = C * np.cos(np.radians(Alpha)) * np.sin(np.radians(Gamma)) * scaling_factor
V3z = np.sqrt(C ** 2 - V3x ** 2 - V3y ** 2) * scaling_factor

# Write the last line to the .gro file
if Alpha == 90.0 and Beta == 90.0 and Gamma == 90.0:
    newstr = str("%.8f" % V1x) + " " + str("%.8f" % V2y) + " " + str("%.8f" % V3z)
else: 
    newstr = str("%.8f" % V1x) + " " + str("%.8f" % V2y) + " " + str("%.8f" % V3z) + " " + str("%.8f" % V1y) + " " + \
             str("%.8f" % V1z) + " " + str("%.8f" % V2x) + " " + str("%.8f" % V2z) + " " + str("%.8f" % V3x) + " " + \
             str("%.8f" % V3y)
commandstring = "sed -i 's/0.00000000   0.00000000   0.00000000/ " + newstr + "/g' " + grofile
os.system(commandstring)

