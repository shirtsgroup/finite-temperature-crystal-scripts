#!/usr/bin/python
#
# CALCULATE THE AVERAGE POSITION OF MOLECULES IN A SERIES OF GRO FILES
# NOTE: THIS SHOULD ONLY BE DONE ON SOLID PHASE CRYSTALS
#
#
# Copyright Michael R. Shirts, University of Virginia, 2014
#

import numpy # numerical array library
import pymbar # multistate Bennett acceptance ratio
from pymbar import timeseries # timeseries analysis
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

parser = OptionParser()
parser.add_option('-G', '--gro', dest = 'grofile', help = 'Gromacs file name', default='frame') #.gro file to be resized
parser.add_option('-M', '--molecule', dest = 'molname', help = 'Name of the molecule in the gro file', default='BNZ')
parser.add_option('-a', '--atoms', dest = 'numatoms', help = 'number of atoms in each molecule', default='12')
parser.add_option('-u', '--unit', dest = 'nummol', help = 'number of molecules in unit cell',default='4')
parser.add_option('-n', '--frames', dest = 'frames', help = 'number of frames to be averaged',default='1000')
(options, args) = parser.parse_args()
mol_name = options.molname
num_atoms = int(options.numatoms)
num_mol = int(options.nummol)
fname = options.grofile
frames=int(options.frames)

#=============================================================================================
# Ensure that all inputs are correct
#=============================================================================================
if mol_name == "":
    print "Please enter an abbreviated molecule name"
    sys.exit()
if num_atoms < 0:
    print "Invalid number of atoms per molecule: " + str(num_atoms)
    sys.exit()
if num_mol < 0:
    print "Invalid number of molecules per unit cell: " + str(num_mol)
    sys.exit()
if fname == "":
    print "Please enter a gro file name"
    sys.exit()

#=============================================================================================
# Generate atom location arrays with all polymorph atom coordinates
#=============================================================================================
#Open the first frame to determine the number of molecules in the system
fname = str(options.grofile)+".gro"
infile = open(fname, 'r')
lines = filter(None, (line.rstrip() for line in infile))
infile.close()
print "loading " + fname
if len(lines[1].split()) == 1:
    # na is the number of total atoms
    na = int(lines[1].split()[0])
    # nm is the number of total molecules
    nm = na/(num_atoms)
else:
    sys.exit('Unexpected .gro file format')

line = 2
mol = 0
boxvect = numpy.zeros(9,float)
old_xyz_position = numpy.zeros((frames,nm,3,nam),float)      #Reduced xyz positions of all the atoms
avg_xyz_position = numpy.zeros((nm,3,nam))      #Final averaged xyz positions of all the atoms
old_crystal_position = numpy.zeros((frames,nm,3,nam),float)      #Original crystal positions of all the atoms
avg_crystal_position = numpy.zeros((nm,3,nam))      #Averaged crystal positions of all the atoms
crystal_basis = numpy.zeros([3,3],float)        #Matrix to take a crystal vector into xyz coordinatestoms


#Loop through each frame of the directory
for i in range(frames+1):
    if i==0:
	continue	
    fname = str(options.grofile)+str(i)+".gro"
    infile = open(fname, 'r')
    lines = filter(None, (line.rstrip() for line in infile))
    infile.close()
    print "loading " + fname
    #Read in the crystal basis matrix (the last line of the .gro file)
    tokens = lines[len(lines)-1].split()
    for j,token in enumerate(tokens):
        if j == 0:
            crystal_basis[0,0]=float(token)
        elif j==1:
            crystal_basis[1,1]=float(token)
        elif j==2:
            crystal_basis[2,2]=float(token)
        elif j==3:
            crystal_basis[0,1]=float(token)
        elif j==4:
            crystal_basis[0,2]=float(token)
        elif j==5:
            crystal_basis[1,0]=float(token)
        elif j==6:
            crystal_basis[1,2]=float(token)
        elif j==7:
            crystal_basis[2,0]=float(token)
        elif j==8:
            crystal_basis[2,1]=float(token)

    crystal_basis = numpy.transpose(crystal_basis)
    xyz_to_crystal = numpy.linalg.inv(crystal_basis)

    #Read in the xyz positions of all the atoms in this gro file
    mol = 0
    line = 2
    while mol < nm:
        acounter = 0
        while acounter < nam:
            old_xyz_position[i-1,mol,0,acounter] = float(lines[line].split()[3])
            old_xyz_position[i-1,mol,1,acounter] = float(lines[line].split()[4])
            old_xyz_position[i-1,mol,2,acounter] = float(lines[line].split()[5])
            line += 1
            acounter +=1
	#Now convert the xyz positions into crystals positions
	#pdb.set_trace()
        old_crystal_position[i-1,mol,:,:]=xyz_to_crystal.dot(old_xyz_position[i-1,mol,:,:])	
	avg_crystal_position[mol,:,:]+=(1.0/frames)*old_crystal_position[i-1,mol,:,:]
        mol += 1
    
#=============================================================================================
# READ IN THE FINAL CRYSTAL BASIS MATRIX
#=============================================================================================
fname = str(options.grofile)+".gro"
infile = open(fname, 'r')
lines = filter(None, (line.rstrip() for line in infile))
infile.close()
print "loading " + fname
#Read in the crystal basis matrix (the last line of the .gro file)
final_vect = lines[len(lines)-1]
tokens = lines[len(lines)-1].split()
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

#=============================================================================================
# CONVERT THE AVERAGED CRYSTAL POSITIONS INTO XYZ POSITIONS AND OUTPUT THE GRO FILE
#=============================================================================================

mol=0
while mol < nm:
    avg_xyz_position[mol,:,:]=crystal_basis.dot(avg_crystal_position[mol,:,:])
    mol += 1
# Write output .gro file
# tacount is the total atom count
fname = str(options.grofile)+"_avg.gro"
outfile = open(fname, 'w')
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
        x = round(avg_xyz_position[mol,0,atom],5)
        y = round(avg_xyz_position[mol,1,atom],5)
        z = round(avg_xyz_position[mol,2,atom],5)
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
outfile.write(final_vect)
outfile.write('\n')  # make sure the file ends on a newline
outfile.close()
