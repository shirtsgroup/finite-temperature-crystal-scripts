# /usr/bin/python
# PYTHON SCRIPT TO CUT THE DEPENDENT BENZENES FROM A GRO FILE
# param f - The name of the .gro file
# param I - The number of independent benzenes

import numpy
from optparse import OptionParser
import os
import subprocess
import pdb

# =============================================================================================
# READ IN USER INPUTS
# =============================================================================================
parser = OptionParser()
parser.add_option('-f', dest='grofile', help='Input gro file', default='frame.gro')
parser.add_option('-n', dest='index', help='index file', default='index.ndx')
(options, args) = parser.parse_args()

# Create the space buffer
spaces = []
Independent_atoms = []
Farthest_copies = []

#rjustvect = numpy.ones(10,int);
#rjustvect = ['8', '7', '5', '15', '15', '15', '6', '6', '6']

fname_gro = options.grofile
fname_xyz = fname_gro[:len(fname_gro)-4] + ".xyz"
fname_key = fname_gro[:len(fname_gro)-4] + "_AMOEBA09.key"
fname_index = options.index

infile = open(fname_index, 'r')
lines = infile.readlines()
infile.close()

#outfile = open(fname, 'w')

#Loop through the index file until the force-averaged lists begin

linecount=1
for i,line in enumerate(lines):
    
    tokens=line.split()

    if i==0 or len(tokens)==0:
        continue #Skip the initial system header and linebreaks

    #Grab the first number and last number from each row and save only those atoms
    if tokens[0]=='[':
	tokens=lines[i+1].split()
	Mirror_atoms=[int(j)+1 for j in tokens]
	break
	Independent_atoms.append(int(tokens[0]))
	Farthest_copies.append(int(tokens[len(tokens)-1]))

print Independent_atoms
print Farthest_copies

#Store the original box vectors
infile = open(fname_gro, 'r')
lines = infile.readlines()
infile.close()

basis = numpy.zeros([3,3],float)
boxV = lines[len(lines)-1].split()
for i,token in enumerate(boxV):
    if i == 0:
        basis[0,0]=float(token)
    elif i==1:
        basis[1,1]=float(token)
    elif i==2:
        basis[2,2]=float(token)
    elif i==3:
        basis[0,1]=float(token)
    elif i==4:
        basis[0,2]=float(token)
    elif i==5:
        basis[1,0]=float(token)
    elif i==6:
        basis[1,2]=float(token)
    elif i==7:
        basis[2,0]=float(token)
    elif i==8:
        basis[2,1]=float(token)

#Store the atomic positions of all the independent atoms and all the farthest copies
Mirror_atoms_xyz=numpy.zeros([len(Mirror_atoms),3], float)
Independent_atoms_xyz=numpy.zeros([len(Independent_atoms),3], float)
Farthest_copies_xyz=numpy.zeros([len(Farthest_copies),3], float)
for i,index in enumerate(Mirror_atoms):
        Mirror_atoms_xyz[i,:] = [float(j) for j in (lines[index+1].split())[3:6]]
for i,index in enumerate(Independent_atoms):
	Independent_atoms_xyz[i,:] = [float(j) for j in (lines[index+1].split())[3:6]]
for i,index in enumerate(Farthest_copies):
        Farthest_copies_xyz[i,:] = [float(j) for j in (lines[index+1].split())[3:6]]

#Calculate the matrix to take the xyz coordinates into crystal basis coordinates.
xyz_to_crystal = numpy.linalg.inv(basis)

#Find the smallest nonzero distance between the first atom and each image in the X,Y,Z direction
Xdist=100
Ydist=100
Zdist=100


#Loop through all mirror atoms and find the smallest distance in each dimension
tol=0.1
for i,index in enumerate(Mirror_atoms):
	v=(Mirror_atoms_xyz[i,:] - Mirror_atoms_xyz[0,:])
	v_crystal=numpy.dot(v,xyz_to_crystal);
	X=numpy.absolute(v_crystal[0])
	Y=numpy.absolute(v_crystal[1])
	Z=numpy.absolute(v_crystal[2])
	#X=numpy.absolute(numpy.dot(v,basis[0,:])/(numpy.dot(basis[0,:],basis[0,:])))
	#Y=numpy.absolute(numpy.dot(v,basis[1,:])/(numpy.dot(basis[1,:],basis[1,:])))
	#Z=numpy.absolute(numpy.dot(v,basis[2,:])/(numpy.dot(basis[2,:],basis[2,:])))


	#print v
	v_crystal=[numpy.round(j,2) for j in v_crystal]
	print v_crystal
	#print "X: " + str(X)
	#pdb.set_trace()
	#If the component is not zero and is smaller than the current smallest value, replace
	if X < Xdist and X > tol:
	    Xdist = X
	if Y < Ydist and Y > tol:
            Ydist = Y
	if Z < Zdist and Z > tol:
            Zdist = Z

#Now flip the distance and round to the nearest integer
if Xdist == 100:
    Xdist = 1.0
else:
    Xdist=numpy.round(1.0/Xdist,0)
if Ydist == 100:
    Ydist = 1.0
else:
    Ydist=numpy.round(1.0/Ydist,0)
if Zdist == 100:
   Zdist = 1.0
else:
   Zdist=numpy.round(1.0/Zdist,0)

print "X: " + str(Xdist) + " Y: " + str(Ydist) + " Z: " + str(Zdist)

#Use linux to cut out the unnecessary lines
cutline = Mirror_atoms[1]-Mirror_atoms[0]+3
sedcommand = "sed -i '" + str(cutline) + ", $d' " + fname_gro
os.system(sedcommand)

#Change the number of atoms at the top of the gro file
oldnumatoms=lines[1].split()[0]
sedcommand = "sed -i '2s/" + oldnumatoms + "/" + str(cutline-3) + "/g' "  + fname_gro
os.system(sedcommand)

#Paste in the revised box vectors based on the supercell dimensions
if len(boxV)==3:
    newVect=numpy.zeros([3],float)
    newVect[0] = str(numpy.round(float(boxV[0])/Xdist,8))
    newVect[1] = str(numpy.round(float(boxV[1])/Ydist,8))
    newVect[2] = str(numpy.round(float(boxV[2])/Zdist,8))
else:
    newVect=numpy.zeros([9],float)
    newVect[0] = str(numpy.round(float(boxV[0])/Xdist,8))
    newVect[1] = str(numpy.round(float(boxV[1])/Ydist,8))
    newVect[2] = str(numpy.round(float(boxV[2])/Zdist,8))
    newVect[3] = str(numpy.round(float(boxV[3])/Xdist,8))
    newVect[4] = str(numpy.round(float(boxV[4])/Xdist,8))
    newVect[5] = str(numpy.round(float(boxV[5])/Ydist,8))
    newVect[6] = str(numpy.round(float(boxV[6])/Ydist,8))
    newVect[7] = str(numpy.round(float(boxV[7])/Zdist,8))
    newVect[8] = str(numpy.round(float(boxV[8])/Zdist,8))
newVect=[str(j) for j in newVect]
newVectstr="   " + '   '.join(newVect)
os.system("echo '" + newVectstr + "' >> " + fname_gro)

#If an xyz file exists for this configuration, cut these atoms out as well
if os.path.isfile(fname_xyz):
    sedcommand = "sed -i '" + str(cutline-1) + ", $d' " + fname_xyz
    os.system(sedcommand)
    oldnumatoms=lines[1].split()[0]
    sedcommand = "sed -i '0,/" + oldnumatoms + "/{s/" + oldnumatoms + "/"  + str(cutline-3) + "/}' "  + fname_xyz
    os.system(sedcommand)

#If a key file exists for this configuration, change the box vectors in this file as well
if os.path.isfile(fname_key):
    infile = open(fname_key, 'r')
    lines = infile.readlines()
    infile.close()    

    for line in lines:
	tokens = line.split()
	if len(tokens) == 0:
	    continue
	if tokens[0] == "A-AXIS":
	    OldX = tokens[1]
	elif tokens[0] == "B-AXIS":
	    OldY = tokens[1]
	elif tokens[0] == "C-AXIS":
	    OldZ = tokens[1]
    sedcommand = "sed -i '0,/" + OldX + "/{s/" + OldX + "/" + str(numpy.round(float(OldX)/Xdist,5)) + "/}' " + fname_key
    os.system(sedcommand)
    sedcommand = "sed -i '0,/" + OldY + "/{s/" + OldY + "/" + str(numpy.round(float(OldY)/Ydist,5)) + "/}' " + fname_key
    os.system(sedcommand)
    sedcommand = "sed -i '0,/" + OldZ + "/{s/" + OldZ + "/" + str(numpy.round(float(OldZ)/Zdist,5)) + "/}' " + fname_key
    os.system(sedcommand)

