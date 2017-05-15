#!/usr/bin/python
#
# Adjust the molecules in a gro file (by jumping them across the box periodic boundaries) such 
# that the molecules have centers of masses closest to a reference gro file
#
# Copyright Michael R. Shirts, University of Virginia, 2014
#

import numpy # numerical array library
from optparse import OptionParser # for parsing command-line options
import pdb

#=============================================================================================
# INPUT VARIABLES
#=============================================================================================
#Gro file to be resized
# polymorph B'~ polymorph B
# mol_name is the name of the molecule in the .gro files
mol_name = 'BNZ'
# nam is the number of atoms per molecule
nam = int(12)
# nmu is the number of molecules per unit cell
nmu = 4

#=============================================================================================
# OUTPUT VARIABLES
#=============================================================================================
# This script produces a gro file with all molecules translated to be in their original starting periodic copy


parser = OptionParser()

parser.add_option('-f', '--gro', dest = 'grofile', help = 'gro file to adjust', default = 'benzene_PROD.gro')
parser.add_option('-r', '--ref', dest = 'reffile', help = 'gro file used as a reference', default = 'benzene_restraint.gro')
(options, args) = parser.parse_args()
groname = options.grofile
refname = options.reffile

#=============================================================================================
# ENSURE THAT ALL INPUTS ARE CORRECT
#=============================================================================================
if groname == "":
    print "Please enter a gro file name"
    sys.exit()
if refname == "":
    print "Please enter a reference file name"
    sys.exit()

#=============================================================================================
# READ IN THE REFERENCE GRO FILE AND DETERMINE CENTERS OF MASSES
#=============================================================================================

# Read in input files
infile = open(refname, 'r')
lines = infile.readlines()
infile.close()

# Read in the total number of atoms in the system (should be the only item on the second line of the .gro file)
num_atoms=12
if len(lines[1].split()) == 1:
    # na is the number of total atoms
    na = int(lines[1].split()[0])
    # nm is the number of total molecules
    nm = na/(num_atoms)
else:
    sys.exit('Unexpected .gro file format')

# Read in atom coordinate data (starts on 3rd line of the .gro file)
# acoordA and acoordB are the atom coordinates for polymorphs A and B
line = 2
mol = 0
xyz_position = numpy.zeros((nm,3,nam))      #Original xyz positions of all the atoms

while mol < nm:
    acounter = 0
    while acounter < nam:
        xyz_position[mol,0,acounter] = float(lines[line].split()[3])
        xyz_position[mol,1,acounter] = float(lines[line].split()[4])
        xyz_position[mol,2,acounter] = float(lines[line].split()[5])
        line += 1
        acounter +=1
    mol += 1

#=============================================================================================
# CALCULATE THE CENTROID OF EACH MOLECULE IN THE REFERENCE FILE
#=============================================================================================

w = numpy.ones((nam,1))*(1.0/nam)
centroid_ref = numpy.zeros((nm,3,1))
centroid_crystal=numpy.zeros((nm,3,1))
for mol in range(nm):
    centroid_ref[mol,:,:] = numpy.dot(xyz_position[mol,:,:], w)

#=============================================================================================
# READ IN THE NEW GRO FILE AND DETERMINE THE CENTERS OF MASSES
#=============================================================================================

# Read in input files
infile = open(groname, 'r')
lines = infile.readlines()
infile.close()

# Read in the total number of atoms in the system (should be the only item on the second line of the .gro file)
if len(lines[1].split()) == 1:
    # na is the number of total atoms
    na = int(lines[1].split()[0])
    # nm is the number of total molecules
    nm = na/(num_atoms)
else:
    sys.exit('Unexpected .gro file format')

# Read in atom coordinate data (starts on 3rd line of the .gro file)
# acoordA and acoordB are the atom coordinates for polymorphs A and B
line = 2
mol = 0
xyz_position = numpy.zeros((nm,3,nam))      #Original xyz positions of all the atoms

while mol < nm:
    acounter = 0
    while acounter < nam:
        xyz_position[mol,0,acounter] = float(lines[line].split()[3])
        xyz_position[mol,1,acounter] = float(lines[line].split()[4])
        xyz_position[mol,2,acounter] = float(lines[line].split()[5])
        line += 1
        acounter +=1
    mol += 1

#=============================================================================================
# CALCULATE THE CENTROID OF EACH MOLECULE IN THE INPUT FILE
#=============================================================================================

w = numpy.ones((nam,1))*(1.0/nam)
centroid_grofile = numpy.zeros((nm,3,1))
for mol in range(nm):
    centroid_grofile[mol,:,:] = numpy.dot(xyz_position[mol,:,:], w)


tokens = lines[len(lines)-1].split()				#Last line of the gro file
crystal_basis=numpy.zeros([3,3],float)
newvect=[]                          
#Read in the old crystal basis matrix (the last line of the .gro file)
for i,token in enumerate(tokens):
    if i == 0:
        crystal_basis[0,0]=float(token)
    elif i==1:
        crystal_basis[1,1]=float(token)
    elif i==2:
        crystal_basis[2,2]=float(token)
    elif i==3:
        crystal_basis[0,1]=float(token)
    elif i==4:
        crystal_basis[0,2]=float(token)
    elif i==5:
        crystal_basis[1,0]=float(token)
    elif i==6:
        crystal_basis[1,2]=float(token)
    elif i==7:
        crystal_basis[2,0]=float(token)
    elif i==8:
        crystal_basis[2,1]=float(token)
    newvect.append(token)


#Loop through all molecules in the gro file and determine the translation
# that puts the center of masses closest together
new_grofile = numpy.zeros((nm,3,1))
new_xyz_position = numpy.zeros((nm,3,nam))      #New xyz positions of all the atoms
for n in range (len(centroid_grofile)):
    distance_Min = 1000000
    for i in range(-3,3):
        for j in range(-3,3):
	    for k in range(-3,3):	
	        distanceVect = centroid_grofile[n,:,0] - centroid_ref[n,:,0] + i*crystal_basis[0,:] + j*crystal_basis[1,:] + k*crystal_basis[2,:]
		distance=distanceVect.dot(distanceVect)
		if distance < distance_Min:
		    distance_Min = distance
		    translationVect = i*crystal_basis[0,:] + j*crystal_basis[1,:] + k*crystal_basis[2,:]
    for a in range(12):
        new_xyz_position[n,:,a] = xyz_position[n,:,a] + translationVect   #Apply the Translation Operation

#=============================================================================================
# Output the new gro file
#=============================================================================================
# Write output .gro file
# tacount is the total atom count
outfile = open('fixed_' + groname, 'w')
outfile.write(lines[0])
#outfile.write('\n')
outfile.write(lines[1])
#outfile.write('\n')
tacount = 0
for mol in range(nm):
    for atom in range(int(nam)):
    # mcount is the current molecule tracker
        mcount = mol + 1
        tacount += 1
        x = round(new_xyz_position[mol,0,atom],5)
        y = round(new_xyz_position[mol,1,atom],5)
        z = round(new_xyz_position[mol,2,atom],5)
	if x < 0.0:
	    xstr = '   ' + "%.5f" % x
	elif x == 0.0:
	    xstr = '    0.00000'
	else:
	    xstr = '    ' + "%.5f" % x
	if y < 0.0:
            ystr = '   ' + "%.5f" % y
	elif y == 0.0:
	    ystr = '    0.00000'
        else:
            ystr = '    ' + "%.5f" % y
	if z < 0.0:
            zstr = '   ' + "%.5f" % z
	elif z == 0.0:
	    zstr = '    0.00000'
        else:
            zstr = '    ' + "%.5f" % z

        #xstr = '%(x) 10g' %vars()
        #ystr = '%(y) 10g' %vars()
        #zstr = '%(z) 10g' %vars()
        line = str(mcount).rjust(5) + mol_name + lines[tacount+1].split()[1].rjust(7) + str(tacount).rjust(5) + xstr + ystr + zstr + '\n'
        outfile.write(line)
outfile.write("   " + '   '.join(newvect))
outfile.write('\n')  # make sure the file ends on a newline
outfile.close()
