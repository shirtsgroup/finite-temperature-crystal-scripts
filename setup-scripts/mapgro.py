#!/usr/bin/python
"""
Python Script to map unit cell of polymorph A to unit cell of polymorph B of a certain crystal utilizing the Kabsch algorithm
References:
1. Kabsch W. A solution for the best rotation to relate two sets of vectors. Acta Cryst A. 1976; 32:922-923
2. Kabsch W. A discussion of the solution for the best rotation to relate two sets of vectors.  Acta Cryst A. 1978; 34: 827-828
3. Wikipedia contributors. 'Kabsch algorithm.' Wikipedia, The Free Encyclopedia.
"""

#==============================================================================
# TO DO
# HIGH PRIORITY
# 
# LOW PRIORITY
# 1. Reduce the number of inputs needed so that all information is gathered from the two .gro files
# (Don't have to enter # of molecules, atoms)
# 2. Instead of writing huge rotation array, write small rotation arrays and dictionary to log file
# 3. Make figures prettier.
#==============================================================================

#==============================================================================
# IMPORTS
#==============================================================================
import sys
import os
import subprocess
import numpy as np
import pdb

#==============================================================================
# INPUT VARIABLES
#==============================================================================
# Two .gro files (polymorph A and polymorph B) each containing the same number of molecules
# Information about the crystal, entered by the user
# pmA is the .gro file for polymorph A.  This is the polymorph that polymorph B' will be mapped from.
#'inputs/benzene_p1_256_toprocessnvt.gro'
pmA = sys.argv[1]
# pmB is the .gro file for polymorph B.  This is the desired end result of the mapped polymorph.
# polymorph B'~ polymorph B
#'inputs/benzene_p3_256_toprocessnvt.gro'
pmB = sys.argv[2]
# pmC is the true .gro file for polymorph A.  This is the desired end result of the mapped polymorph.
# polymorph B'~ polymorph B
#'inputs/benzene_p3_256_toprocessnvt.gro'
pmC = sys.argv[3]
# mol_name is the name of the molecule in the .gro files
mol_name = 'BNZ'
# nam is the number of atoms per molecule
nam = int(12)
# nmu is the number of molecules per unit cell
nmu = 4

#==============================================================================
# OUTPUT VARIABLES
#==============================================================================
# This script produces three output files
#runid = 'p1_p3_prod'
runid = sys.argv[4]

output = runid + '-mapped.gro'

#==============================================================================
# AUXILIARY FUNCTION
#==============================================================================

def computeboxvol(boxv):
    if len(boxv) == 3:

        # volume is simple for rectangular box
        volume = boxv[0]*boxv[1]*boxv[2]

    elif len(boxv) == 9:
        # more complicated for a more general box

        # from gromacs v1(x) v2(y) v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y)
        boxm = np.zeros([3,3],float)
        boxm[0,0] = boxv[0]
        boxm[0,1] = boxv[3]
        boxm[0,2] = boxv[4]
        boxm[1,0] = boxv[5]
        boxm[1,1] = boxv[1]
        boxm[1,2] = boxv[6]
        boxm[2,0] = boxv[7]
        boxm[2,1] = boxv[8]
        boxm[2,2] = boxv[2]
        volume = np.dot(boxm[0,:],np.cross(boxm[1,:],boxm[2,:]))
    else:
        print 'Vector not defined correctly!'
        volume = 0
    return volume

#==============================================================================
# Generate atom location arrays with all polymorph atom coordinates
#==============================================================================
# Read in input files
fileA = open(pmA, 'r')
linesA = filter(None, (line.rstrip() for line in fileA))
fileA.close()
fileB = open(pmB, 'r')
linesB = filter(None, (line.rstrip() for line in fileB))
fileB.close()
fileC = open(pmC, 'r')
linesC = filter(None, (line.rstrip() for line in fileC))
fileC.close()

# Read in the total number of atoms in the system (should be the only item on the second line of the .gro file)
if len(linesA[1].split()) == 1:
    # na is the number of total atoms
    na = int(linesA[1].split()[0])
    # nm is the number of total molecules
    nm = na/nam
else:
    sys.exit('Unexpected .gro file format')
if int(linesB[1].split()[0]) != na:
    sys.exit('Number of molecules changes between polymorphs')
else:
    pass

# Read in atom coordinate data (starts on 3rd line of the .gro file)
# acoordA and acoordB are the atom coordinates for polymorphs A and B
line = 2
mol = 0
acoordA = np.zeros((nm,3,nam))
acoordB = np.zeros((nm,3,nam))
acoordC = np.zeros((nm,3,nam))
while mol < nm:
    acounter = 0
    while acounter < nam:
        acoordA[mol,0,acounter] = float(linesA[line].split()[3])
        acoordA[mol,1,acounter] = float(linesA[line].split()[4])
        acoordA[mol,2,acounter] = float(linesA[line].split()[5])
        
        acoordB[mol,0,acounter] = float(linesB[line].split()[3])
        acoordB[mol,1,acounter] = float(linesB[line].split()[4])
        acoordB[mol,2,acounter] = float(linesB[line].split()[5])

	acoordC[mol,0,acounter] = float(linesC[line].split()[3])
        acoordC[mol,1,acounter] = float(linesC[line].split()[4])
        acoordC[mol,2,acounter] = float(linesC[line].split()[5])

        line += 1
        acounter +=1
    mol += 1

#==============================================================================
# Calculate the centroid of each molecule
#==============================================================================

w = np.ones((nam,1))*(1.0/nam)
centroidA = np.zeros((nm,3,1))
centroidB = np.zeros((nm,3,1))
centroidC = np.zeros((nm,3,1))
transC = np.zeros((nm,3,1)) #Translate the benzene configuration to most closely match the original mapping configuration
boxA = np.array(linesA[na+2].split()).astype(np.float)
box_vectors = np.zeros([3,3],float)
if len(boxA)==3:
    box_vectors[0,0] = boxA[0]
    box_vectors[0,1] = 0
    box_vectors[0,2] = 0
    box_vectors[1,0] = 0
    box_vectors[1,1] = boxA[1]
    box_vectors[1,2] = 0
    box_vectors[2,0] = 0
    box_vectors[2,1] = 0
    box_vectors[2,2] = boxA[2]
else:
    box_vectors[0,0] = boxA[0]
    box_vectors[0,1] = boxA[3]
    box_vectors[0,2] = boxA[4]
    box_vectors[1,0] = boxA[5]
    box_vectors[1,1] = boxA[1]
    box_vectors[1,2] = boxA[6]
    box_vectors[2,0] = boxA[7]
    box_vectors[2,1] = boxA[8]
    box_vectors[2,2] = boxA[2]
for mol in range(nm):
    centroidA[mol,:,:] = np.dot(acoordA[mol,:,:], w)
    centroidB[mol,:,:] = np.dot(acoordB[mol,:,:], w)
    centroidC[mol,:,:] = np.dot(acoordC[mol,:,:], w)
    rmsd=100
    
    #translate the molecule to find the closest periodic copy to the original mapping configuration
    for i in range(-1,2):
	for j in range(-1,2):
	    for k in range(-1,2):
		v=[[i],[j],[k]]
		centroid=np.dot(acoordC[mol,:,:], w) +np.dot(box_vectors,v)
		temp_rmsd = (centroid[0]-centroidA[mol,0,:])**2 + (centroid[1]-centroidA[mol,1,:])**2 + (centroid[2]-centroidA[mol,2,:])**2
		if temp_rmsd < rmsd:
		    rmsd = temp_rmsd
		    transC[mol,:,:]=np.dot(box_vectors,v)
    

#==============================================================================
# Calculate unit cell mapping
#==============================================================================
# Translate each atom to be centered at the origin
# cacoordA and cacoordB are the centered atom coordinate arrays
cacoordA = acoordA - np.ones((nm,3,nam))*centroidA
cacoordB = acoordB - np.ones((nm,3,nam))*centroidB


# Compute the covariance matrix going from Polymorph A to Polymorph B for each individual molecule
diag = np.identity(int(nam))*(1.0/nam)
cov = np.zeros((nm,3,3))
for mol in range(nm):
    cov[mol,:,:]=np.dot(cacoordB[mol,:,:], np.transpose(cacoordA[mol,:,:]))


# Computation of the optimal rotation matrix
     # This can be done using singular value decomposition (SVD)
     # Getting the sign of the det(V)*(W) to decide
     # whether we need to correct our rotation matrix to ensure a
     # right-handed coordinate system.
     # And finally calculating the optimal rotation matrix U
     # see http://en.wikipedia.org/wiki/Kabsch_algorithm

# singular value decomposition of the covariance matrices
##V = np.zeros((nm,3,3))
#S = np.zeros((nm,3,3))
#W = np.zeros((nm,3,3))
rot = np.zeros((nm,3,3))
trans = np.zeros((nm,3,nam))
for mol in range(nm):
    V, S, W = np.linalg.svd(cov[mol,:,:])
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0
    if d:
         S[-1] = -S[-1]
         V[:, -1] = -V[:, -1]
    rot[mol,:,:]=np.dot(V,W)

writerotA = np.zeros((nm,3,3))
writerotB = np.zeros((nm,3,3))
#Translational mapper
trans = np.zeros((nm,3,nam))
# Map unit cell A to unit cell B'
acoordBprime = np.zeros((nm,3,nam))
for mol in range(nm):
        trans[mol,:,:] = centroidB[mol,:,:]-np.dot(rot[mol,:,:],centroidA[mol,:,:])
        acoordBprime[mol,:,:] = np.dot(rot[mol,:,:],(acoordC[mol,:,:]+transC[mol,:,:]))+trans[mol,:,:]

#==============================================================================
# Write Output Files
#==============================================================================
logfile= open('jacobian.txt', 'w')
# Jacobian for each molecule is the box volume transformation, times the number of molecules 
boxA = np.array(linesA[na+2].split()).astype(np.float)
boxB = np.array(linesB[na+2].split()).astype(np.float)
volA = computeboxvol(boxA)
volB = computeboxvol(boxB)
Jacobian = nm*np.log(volB/volA)
logfile.write(str(Jacobian)+'\n') 
# Write output .gro file
outfile= open(output, 'w')
# tacount is the total atom count
outfile.write(linesC[0])
outfile.write('\n')
outfile.write(linesC[1])
outfile.write('\n')
tacount = 0
for mol in range(nm):
    for atom in range(int(nam)):
    # mcount is the current molecule tracker
        mcount = mol + 1
        tacount += 1
        x = round(acoordBprime[mol,0,atom],9)
        y = round(acoordBprime[mol,1,atom],9)
        z = round(acoordBprime[mol,2,atom],9)

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

        line = str(mcount).rjust(5) + mol_name + linesA[tacount+1].split()[1].rjust(7) + str(tacount).rjust(5) + xstr + ystr + zstr + '\n'
        outfile.write(line)	


outfile.write(linesB[len(linesB)-1])
outfile.write('\n')  # make sure the file ends on a newline
outfile.close()
logfile.close()

