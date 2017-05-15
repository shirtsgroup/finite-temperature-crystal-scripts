#!/usr/bin/python
#
# PYTHON SCRIPT TO TRANSFORM ALL BENZENE MOLECULES IN A GIVEN GRO FILE INTO PLANAR
# MOLECULES WITH THE EQUILIBRIUM C-C AND C-H BOND LENGTHS
#
# Copyright Michael R. Shirts and Eric Dybeck, University of Virginia, 2015
#

import numpy # numerical array library
from optparse import OptionParser # for parsing command-line options
import pdb

#=============================================================================================
# INPUT VARIABLES
#=============================================================================================
#Gro file to be resized
#C-C Bond Length
#C-H Bond Length

#=============================================================================================
# OUTPUT VARIABLES
#=============================================================================================
# This script produces a gro file with the transformed molecules


parser = OptionParser()
parser.add_option('-g', '--gro', dest = 'grofile', help = 'Gromacs File') #.gro file to be resized
parser.add_option('-C', '--carbon', dest = 'CC', help = 'C-C bond length') # Equilibrium Carbon-Carbon bond length
parser.add_option('-H', '--hydrogen', dest = 'CH', help = 'C-H bond length') # Equilibrium Carbon-Hydrogen bond length
(options, args) = parser.parse_args()
CC_length = float(options.CC)
CH_length = float(options.CH)
mol_name = 'BNZ'
num_atoms = 12
num_mol = 4
fname = options.grofile

#=============================================================================================
# Ensure that all inputs are correct
#=============================================================================================
if num_mol < 0:
    print "Invalid number of molecules per unit cell: " + str(num_mol)
    sys.exit()
if fname == "":
    print "Please enter a gro file name"
    sys.exit()

#=============================================================================================
# Generate atom location arrays with all polymorph atom coordinates
#=============================================================================================
# Read in input files
fname = options.grofile
infile = open(fname, 'r')
lines = filter(None, (line.rstrip() for line in infile))
infile.close()
print "loading " + fname

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
nam = num_atoms
xyz_position = numpy.zeros((nm,3,nam))	#Original xyz positions of all the atoms
new_xyz_position = numpy.zeros((nm,3,nam))  #New xyz positions of all the atoms

while mol < nm:
    acounter = 0
    while acounter < nam:
        xyz_position[mol,0,acounter] = float(lines[line].split()[3])
        xyz_position[mol,1,acounter] = float(lines[line].split()[4])
        xyz_position[mol,2,acounter] = float(lines[line].split()[5])
        line += 1
        acounter +=1
    mol += 1

#Store the final line
finalline = lines[len(lines)-1]

#=============================================================================================
# DEFINE VARIABLES AND PERFORM THE NECESSARY TRANSFORMATIONS
#=============================================================================================
COM = numpy.zeros(3,float) #XYZ position of the Center of Mass
C1 = numpy.zeros(3,float) #XYZ position of the first carbon
C3 = numpy.zeros(3,float) #XYZ position of the third carbon
C5 = numpy.zeros(3,float) #XYZ position of the fifth carbon
C_1_3 = numpy.zeros(3,float) #Vector connecting the first and third carbon
C_1_5 = numpy.zeros(3,float) #vector connecting the first and fifth carbon
C_z = numpy.zeros(3,float) #Normal vector to the plane containing the carbon
C_x = numpy.zeros(3,float) #vector connecting the COM to the first carbon
C_y = numpy.zeros(3,float) #vector normal to C_z and C_y 
Internal_coords = numpy.zeros([nam,3]) #Matrix containing the position of the atoms in internal coordinates
Basis = numpy.zeros([3,3]) #Basis matrix for the XYZ coordinates of the internal coordinate vectors
acounter = 0



for mol in range(nm):
    #=============================================================================================
    # CALCULATE THE CENTER OF MASS OF EACH MOLECULE
    #=============================================================================================
    w = numpy.ones((nam,1))*(1.0/nam)
    COM_flipped = numpy.dot(xyz_position[mol,:,:], w)
    COM[0] = COM_flipped[0]
    COM[1] = COM_flipped[1]
    COM[2] = COM_flipped[2]

    #=============================================================================================
    # DETERMINE THE XYZ POSITION OF THE 3 PRINCIPLE CARBONS
    #=============================================================================================
    C1 = xyz_position[mol,:,acounter]
    C3 = xyz_position[mol,:,acounter+4]
    C5 = xyz_position[mol,:,acounter+8]

    #=============================================================================================
    # DETERMINE THE 'INTERNAL COORDINATE VECTORS'
    #=============================================================================================
    C_1_3 = C3 - C1
    C_1_5 = C5 - C1
    C_z = numpy.cross(C_1_3, C_1_5)
    C_z = C_z / (C_z[0]**2 + C_z[1]**2+C_z[2]**2)**0.5
    COM_135 = COM + C_z*(numpy.dot(C_z,C1-COM))
    C_x = C1 - COM_135
    C_x = C_x / (C_x[0]**2 + C_x[1]**2+C_x[2]**2)**0.5
    C_y = numpy.cross(C_z,C_x)
    Basis[0,:] = C_x
    Basis[1,:] = C_y
    Basis[2,:] = C_z

    #=============================================================================================
    # FILL IN THE INTERNAL COORDINATES OF THE BENZENE ATOMS
    #=============================================================================================
    Internal_coords[0,:] = [CC_length,0,0]
    Internal_coords[1,:] = [CC_length+CH_length,0,0]
    Internal_coords[2,:] = [0.5*CC_length, 0.5*(3**0.5)*CC_length,0]
    Internal_coords[3,:] = [0.5*(CC_length+CH_length), 0.5*(3**0.5)*(CC_length+CH_length),0] 
    Internal_coords[4,:] = [-0.5*CC_length, 0.5*(3**0.5)*CC_length, 0] 
    Internal_coords[5,:] = [-0.5*(CC_length+CH_length), 0.5*(3**0.5)*(CC_length+CH_length), 0]
    Internal_coords[6,:] = [-1*CC_length,0,0]
    Internal_coords[7,:] = [-1*(CC_length+CH_length),0,0]
    Internal_coords[8,:] = [-0.5*CC_length, -0.5*(3**0.5)*CC_length,0]
    Internal_coords[9,:] = [-0.5*(CC_length+CH_length), -0.5*(3**0.5)*(CC_length+CH_length),0]
    Internal_coords[10,:] = [0.5*CC_length, -0.5*(3**0.5)*CC_length,0]
    Internal_coords[11,:] = [0.5*(CC_length+CH_length), -0.5*(3**0.5)*(CC_length+CH_length),0]

    #=============================================================================================
    # CONVERT THE INTERNAL COORDINATES INTO XYZ COORDINATES
    #=============================================================================================
    #print (numpy.dot(Internal_coords, Basis) + COM).transpose()
    #pdb.set_trace()
    new_xyz_position[mol,:,acounter:(acounter+12)] = (numpy.dot(Internal_coords, Basis) + COM).transpose()
acounter += nam

#=============================================================================================
# Output the new gro file
#=============================================================================================
# Write output .gro file
# tacount is the total atom count
outfile = open('relax_' + fname, 'w')
outfile.write(lines[0])
outfile.write('\n')
outfile.write(lines[1])
outfile.write('\n')
tacount = 0
for mol in range(nm):
    for atom in range(int(nam)):
    # mcount is the current molecule tracker
        mcount = mol + 1
        tacount += 1
        x = round(new_xyz_position[mol,0,atom],9)
        y = round(new_xyz_position[mol,1,atom],9)
        z = round(new_xyz_position[mol,2,atom],9)
	if x < 0.0:
	    xstr = '   ' + "%.9f" % x
	elif x == 0.0:
	    xstr = '    0.00000'
	else:
	    xstr = '    ' + "%.9f" % x
	if y < 0.0:
            ystr = '   ' + "%.9f" % y
	elif y == 0.0:
	    ystr = '    0.00000'
        else:
            ystr = '    ' + "%.9f" % y
	if z < 0.0:
            zstr = '   ' + "%.9f" % z
	elif z == 0.0:
	    zstr = '    0.00000'
        else:
            zstr = '    ' + "%.9f" % z

        #xstr = '%(x) 10g' %vars()
        #ystr = '%(y) 10g' %vars()
        #zstr = '%(z) 10g' %vars()
        line = str(mcount).rjust(5) + mol_name + lines[tacount+1].split()[1].rjust(7) + str(tacount).rjust(5) + xstr + ystr + zstr + '\n'
        outfile.write(line)
outfile.write(finalline)
outfile.write('\n')  # make sure the file ends on a newline
outfile.close()
