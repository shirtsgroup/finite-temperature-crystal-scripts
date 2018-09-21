#!/usr/bin/python
#
# Remove spurious COM motion by pinning the COM in the center of the box
#
# Copyright Michael R. Shirts, University of Virginia, 2014
#

import numpy as np
from optparse import OptionParser
import sys


# =============================================================================================
# INPUT VARIABLES
# =============================================================================================
# Gro file to be adjusted
# mol_name is the name of the molecule in the .gro files
mol_name = 'BNZ'
# nam is the number of atoms per molecule
nam = int(12)
# nmu is the number of molecules per unit cell
nmu = 4

# =============================================================================================
# OUTPUT VARIABLES
# =============================================================================================
# This script produces a gro file with all molecules translated so that the COM is at the box center
parser = OptionParser()
parser.add_option('-f', '--gro', dest='grofile', help='gro file to adjust', default='frame1.gro')
(options, args) = parser.parse_args()
groname = options.grofile
molecule = groname[0:6]
if molecule == "benzen":
    nam = int(12)
elif molecule == "formam":
    nam = int(6)
elif molecule == "acetac":
    nam = int(8)
elif molecule == "imazol":
    nam = int(9)
elif molecule == "glycin":
    nam = int(10)
elif molecule == "hxacan":
    nam = int(20)


# =============================================================================================
# ENSURE THAT ALL INPUTS ARE CORRECT
# =============================================================================================
if groname == "":
    print "Please enter a gro file name"
    sys.exit()

# =============================================================================================
# READ IN THE GRO FILE AND DETERMINE THE CENTERS OF MASSES
# =============================================================================================

# Read in input files
infile = open(groname, 'r')
lines = infile.readlines()
infile.close()

nam = 0
for i, line in enumerate(lines):
    tokens = line.split()
    if i < 2:
        continue
    if tokens[0][0] != "1":
        break
    else:
        nam += 1

# Read in the total number of atoms in the system (should be the only item on the second line of the .gro file)
if len(lines[1].split()) == 1:
    # na is the number of total atoms
    na = int(lines[1].split()[0])
    # nm is the number of total molecules
    nm = na / nam
else:
    sys.exit('Unexpected .gro file format')

# Read in atom coordinate data (starts on 3rd line of the .gro file)
# acoordA and acoordB are the atom coordinates for polymorphs A and B
line = 2
mol = 0
xyz_position = np.zeros((nm, 3, nam))  # Original xyz positions of all the atoms

while mol < nm:
    acounter = 0
    while acounter < nam:
        xyz_position[mol, 0, acounter] = float(lines[line].split()[3])
        xyz_position[mol, 1, acounter] = float(lines[line].split()[4])
        xyz_position[mol, 2, acounter] = float(lines[line].split()[5])
        line += 1
        acounter += 1
    mol += 1

# =============================================================================================
# CALCULATE THE CENTROID OF EACH MOLECULE IN THE INPUT FILE
# =============================================================================================

w = np.ones((nam, 1)) * (1.0 / nam)
centroid_grofile = np.zeros((nm, 3, 1))
for mol in range(nm):
    centroid_grofile[mol, :, :] = np.dot(xyz_position[mol, :, :], w)

# Now calculate the centroid of the entire system
w = np.ones((1, nm)) * (1.0 / nm)
OldCOM = np.dot(w, centroid_grofile[:, :, 0])

tokens = lines[len(lines) - 1].split()  # Last line of the gro file
crystal_basis = np.zeros([3, 3], float)
newvect = []
# Read in the crystal basis matrix (the last line of the .gro file)
for i, token in enumerate(tokens):
    if i == 0:
        crystal_basis[0, 0] = float(token)
    elif i == 1:
        crystal_basis[1, 1] = float(token)
    elif i == 2:
        crystal_basis[2, 2] = float(token)
    elif i == 3:
        crystal_basis[0, 1] = float(token)
    elif i == 4:
        crystal_basis[0, 2] = float(token)
    elif i == 5:
        crystal_basis[1, 0] = float(token)
    elif i == 6:
        crystal_basis[1, 2] = float(token)
    elif i == 7:
        crystal_basis[2, 0] = float(token)
    elif i == 8:
        crystal_basis[2, 1] = float(token)
    newvect.append(token)

# Determine the center of the box
TrueCOM = np.dot(0.5 * np.ones(3, float), crystal_basis)

# Loop through all molecules in the gro file and determine the translation
#   that puts the center of masses closest together
new_grofile = np.zeros((nm, 3, 1))
new_xyz_position = np.zeros((nm, 3, nam))  # New xyz positions of all the atoms
for n in range(len(centroid_grofile)):
    for a in range(nam):
        new_xyz_position[n, :, a] = xyz_position[n, :, a] + (TrueCOM-OldCOM)  # Apply the Translation Operation

# =============================================================================================
# Output the new gro file
# =============================================================================================
# Write output .gro file
# tacount is the total atom count
outfile = open(groname, 'w')
outfile.write(lines[0])
outfile.write(lines[1])
tacount = 0
for mol in range(nm):
    # mcount is the current molecule tracker
    for atom in range(int(nam)):
        mcount = mol + 1
        tacount += 1
        x = round(new_xyz_position[mol, 0, atom], 8)
        y = round(new_xyz_position[mol, 1, atom], 8)
        z = round(new_xyz_position[mol, 2, atom], 8)
    if x < 0.0:
        xstr = '   ' + "%.8f" % x
    elif x == 0.0:
        xstr = '    0.00000'
    else:
        xstr = '    ' + "%.8f" % x
    if y < 0.0:
        ystr = '   ' + "%.8f" % y
    elif y == 0.0:
        ystr = '    0.00000'
    else:
        ystr = '    ' + "%.8f" % y
    if z < 0.0:
        zstr = '   ' + "%.8f" % z
    elif z == 0.0:
        zstr = '    0.00000'
    else:
        zstr = '    ' + "%.8f" % z

    #xstr = '%(x) 10g' %vars()
    #ystr = '%(y) 10g' %vars()
    #zstr = '%(z) 10g' %vars()
    line = str(mcount).rjust(5) + mol_name + lines[tacount + 1].split()[1].rjust(7) + str(tacount).rjust(5) + xstr \
           + ystr + zstr + '\n'
    outfile.write(line)
outfile.write("   " + '   '.join(newvect))
outfile.write('\n')  # make sure the file ends on a newline
outfile.close()
