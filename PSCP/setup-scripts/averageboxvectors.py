# /usr/bin/python
# PYTHON SCRIPT TO TAKE IN AN XVG FILE OF BOX VECTORS AT EACH STEP AND RETURN AN AVERAGE
# NOTE: The Box ZX number is taken as the absolute value due to flipping in Benzene 3
import numpy as np
from optparse import OptionParser
import subprocess

#=============================================================================================
# READ IN USER INPUTS
#=============================================================================================
parser = OptionParser()
parser.add_option('-f', dest='xvgfile', help='Input xvg file')
parser.add_option('-g', dest='defaultgro', help='Default box vector used to match the sign', default='none')
(options, args) = parser.parse_args()

ignore_symbols = ['#', '@', '@TYPE', 'STEP', '=====================', 'TIME(ps)']
defaultvect = []

# Returning the sign of the values for the the box vectors
if options.defaultgro != 'none':
    # Reading in the default Gromacs coordinate file
    grofile = str(options.defaultgro)
    infile = open(grofile, 'r')
    lines = infile.readlines()
    infile.close()

    # Grabbing the last line of the Gromacs file, which contains the lattice parameters
    defaultvect = lines[len(lines) - 1].split()
    defaultvect = [float(j) for j in defaultvect]

    for i, X in enumerate(defaultvect):
        if X != 0:
            defaultvect[i] = np.round(np.abs(X)/X, 0)
    # Remove the extraneous terms
    if len(defaultvect) == 9:
        defaultvect = [defaultvect[j] for j in [0, 1, 2, 5, 7, 8]]
    print defaultvect


# Opening the xvgfile containing the box vectors
xvgfile = str(options.xvgfile)
infile = open(xvgfile, 'r')
lines = infile.readlines()
infile.close()
print "loading " + str(xvgfile)

# Making a place to store the box vectors
vectors = np.zeros([6, 100000])

# Setting tolerance for the value of the box vector
tol = 0.01

# Moving through each line of the xvg file
counter = 0
for line in lines:
    tokens = line.split()
    # Checking to see if the line can be ignored
    if tokens[0] in ignore_symbols:
        continue

    # Converting the string to a float
    boxvect = [float(j) for j in tokens[1:]]

    # If the value is close to zero, set it to exactly zero
    for i, X in enumerate(boxvect):
        if np.abs(X) < tol:
            boxvect[i] = 0.0000

    # If a default vector was given, ensure that the sign is the same
    if not defaultvect == []:
        boxvect = [np.abs(j) for j in boxvect]
        for i in range(len(defaultvect)):
            boxvect[i] = boxvect[i] * defaultvect[i]
    vectors[:len(boxvect), counter] = boxvect
    counter += 1

# Averaging the box vectors
averagevect = np.zeros(len(boxvect))
for i in range(len(boxvect)):
    averagevect[i] = np.average(vectors[i, :counter])

# Ouputing the box vectors to a string and storing in a
averagevect = [str(j) for j in averagevect]
averagevectstr = '   '.join(averagevect)
subprocess.call(['echo "' + averagevectstr + '" > avgboxvect.xvg'], shell=True)
