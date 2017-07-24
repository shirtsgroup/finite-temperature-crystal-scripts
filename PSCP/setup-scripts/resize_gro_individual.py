#!/usr/bin/python
#
# Resize a gro file by isotropically resizing its box vectors and translating the center of mass of each molecule
# The input can either be a new total volume or a new box vector
# NOTE: THE VECTOR MUST BE ENCASED IN QUOTATION MARKS WHEN ENTERED 
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
# This script produces a gro file with the resized dimensions and molecules


parser = OptionParser()
parser.add_option('-V', '--volume', dest = 'volume', help = 'Final Volume', default='-1.0') #Final volume in cubic angstroms
parser.add_option('-f', '--gro', dest = 'grofile', help = 'Gromacs File') #.gro file to be resized
parser.add_option('-M', '--molecule', dest = 'molname', help = 'Name of the molecule in the gro file', default='benzene')
parser.add_option('-n', '--atoms', dest = 'numatoms', help = 'number of atoms in each molecule', default='12')
parser.add_option('-u', '--unit', dest = 'nummol', help = 'number of molecules in unit cell',default='4')
parser.add_option('-v', '--vector', dest = 'vector', help = 'New box vector', default = '0 0 0 0 0 0 0 0 0')
(options, args) = parser.parse_args()
final_volume = float(options.volume)
molecule = options.molname
num_atoms = int(options.numatoms)
num_mol = int(options.nummol)
boxvect = options.vector
fname = options.grofile


#=============================================================================================
# DETERMINE THE MOLECULE FROM THE INPUT FILES
#=============================================================================================
if [ -f "keyfile.key" ]; then
   #Get the molecule from the key file
   molecule=$(less keyfile.key | grep "Molecule: " | awk '{print $3}')
else
   #Get the molecule from the topology file
   molecule=$(less topology.top | grep "solid" | awk '{print $1}')
fi
if [ "$molecule" == "benzene" ]; then
    nam=12
elif [ "$molecule" == "glycin" ]; then
    nam=10
elif [ "$molecule" == "acetac" ]; then
    nam=8
elif [ "$molecule" == "hxacan" ]; then
    nam=20
elif [ "$molecule" == "formam" ]; then
    nam=6
elif [ "$molecule" == "imazol" ]; then
    nam=9
elif [ "$molecule" == "cafine" ]; then
    nam=24
elif [ "$molecule" == "zzzvye" ]; then
    nam=15
elif [ "$molecule" == "dohfem" ]; then
    nam=14
elif [ "$molecule" == "bismev" ]; then
    nam=20
elif [ "$molecule" == "cbmzpn" ]; then
    nam=30
elif [ "$molecule" == "pyrzin" ]; then
    nam=14
elif [ "$molecule" == "kobfud" ]; then
    nam=15
else
    echo "Unrecognized molecule: $molecule"
    exit 1
fi


#=============================================================================================
# RUN THE ALGORITHM TO CREATE THE NEW BOX VECTOR
#=============================================================================================
def changeBoxVector(fname="pre_EQ.gro",molecule='benzene',volume=0.0,boxvect='0 0 0 0 0 0 0 0 0', outfile='None'):

    #Ensure that all inputs are correct
    checkInputs(fname,molecule,volume,boxvect)

    if outfile == 'None':
	outfile=fname

    #nam = grabNAM(molecule)

    #Generate arary of old atom locations
    old_xyz_position = generateAtomLocations(fname,nam)

    infile = open(fname, 'r')
    lines = filter(None, (line.rstrip() for line in infile))
    infile.close()

    #Grab the old crystal basis and generate the new crystal basis
    old_Basis = vectorToBasis(lines[len(lines)-1].split())
    if volume > 0.0:
        initial_volume = float(numpy.linalg.det(old_Basis))
        scaling_param = (final_volume/initial_volume)**(1.0/3.0)
        new_Basis = scaling_param*old_Basis
    else:
	new_Basis = vectorToBasis(boxvect.split())
    newboxvect = basisToVector(new_Basis)

    #Create the transformation operator between basis sets
    transformation_Basis=numpy.transpose(numpy.dot(new_Basis,numpy.linalg.inv(old_Basis)))

    #Generate the new atom locations
    new_xyz_position = generateNewAtomLocations(old_xyz_position,transformation_Basis,nam)

    #Output the gro file
    newboxstring = [str(j) for j in newboxvect]
    outputGroFile(outfile,new_xyz_position,newboxstring,fname) 
     


#=============================================================================================
# ENSURE THAT ALL INPUTS ARE CORRECT
#=============================================================================================
def checkInputs(fname,molecule,volume,boxvect):
    moleculeList = ['benzene', 'acetac', 'formam', 'imazol', 'glycin', 'hxacan']
    if  volume <= 0.0 and boxvect == '':
        print "Invalid final volume: " + str(volume)
        sys.exit()
    if molecule not in moleculeList:
	print "Invalid molecule: " + molecule
        print "Please enter a valid molecule name"
        sys.exit()
    if fname == "":
        print "Please enter a gro file name"
        sys.exit()


#=============================================================================================
# DETERMINE THE NUMBER OF ATOMS PER MOLECULE
#=============================================================================================
def grabNAM(molecule):

    if molecule=="benzene":
	return 12
    elif molecule=="acetac":
	return 9
    elif molecule=="formam":
	return 6
    elif molecule=="imazol":
	return 9
    else:
	return -1

#=============================================================================================
# GENERATE ARRAY OF OLD ATOM LOCATIONS
#=============================================================================================
def generateAtomLocations(fname,nam):
    
    # Read in input files
    infile = open(fname, 'r')
    lines = filter(None, (line.rstrip() for line in infile))
    infile.close()
    print "loading " + fname
    
    # Read in the total number of atoms in the system (should be the only item on the second line of the .gro file)
    if len(lines[1].split()) == 1:
        # na is the number of total atoms
        na = int(lines[1].split()[0])
        # nm is the number of total molecules
        nm = na/(nam)
    else:
        sys.exit('Unexpected .gro file format')
    
    # Read in atom coordinate data (starts on 3rd line of the .gro file)
    # acoordA and acoordB are the atom coordinates for polymorphs A and B
    line = 2
    mol = 0
    old_xyz_position = numpy.zeros((nm,3,nam))	#Original xyz positions of all the atoms
    
    while mol < nm:
        acounter = 0
        while acounter < nam:
            old_xyz_position[mol,0,acounter] = float(lines[line].split()[3])
            old_xyz_position[mol,1,acounter] = float(lines[line].split()[4])
            old_xyz_position[mol,2,acounter] = float(lines[line].split()[5])
            line += 1
            acounter +=1
        mol += 1
    return old_xyz_position

#=============================================================================================
# CONVERT A BOX VECTOR INTO A BASIS MATRIX
#=============================================================================================
def vectorToBasis(vector):
    
    basis = numpy.zeros([3,3],float)        #Matrix to transform internal crystal coordinates into xyz coordinates
   
    for i,num in enumerate(vector):
        if i == 0:
            basis[0,0]=float(num)
        elif i==1:
            basis[1,1]=float(num)
        elif i==2:
            basis[2,2]=float(num)
        elif i==3:
            basis[0,1]=float(num)
        elif i==4:
            basis[0,2]=float(num)
        elif i==5:
            basis[1,0]=float(num)
        elif i==6:
            basis[1,2]=float(num)
        elif i==7:
            basis[2,0]=float(num)
        elif i==8:
            basis[2,1]=float(num)
	
    return basis

#=============================================================================================
# CONVERT A BASIS MATRIX INTO A BOX VECTOR
#=============================================================================================
def basisToVector(basis):
    vector = numpy.array([basis[0,0], basis[1,1], basis[2,2], basis[0,1], basis[0,2], basis[1,0], basis[1,2], basis[2,0], basis[2,1]])
    if all(j == 0 for j in vector[3:9]):
	vector = numpy.array(vector[0:3])

    return vector

#=============================================================================================
# GENERATE THE NEW ATOM LOCATIONS
#=============================================================================================
def generateNewAtomLocations(old_xyz_position, transformation_basis,nam):
    #Calculate the centroid of each molecule
    w = numpy.ones((nam,1))*(1.0/nam)
    nm = len(old_xyz_position[:,0,0])
    centroid_xyz = numpy.zeros((nm,3,1))
    for mol in range(nm):
        centroid_xyz[mol,:,:] = numpy.dot(old_xyz_position[mol,:,:], w)

    mol = 0
    new_xyz_position=old_xyz_position.copy()
    while mol < nm:
        new_xyz_position[mol,:,:] = (old_xyz_position[mol,:,:]-centroid_xyz[mol,:,:])+numpy.dot(transformation_basis,centroid_xyz[mol,:,:])
        mol += 1

    return new_xyz_position


#=============================================================================================
# OUTPUT THE NEW GRO FILE
#=============================================================================================
def outputGroFile(outname,new_xyz_position,newboxvect,fname):
    
    #Read in the original file
    infile = open(fname, 'r')
    lines = filter(None, (line.rstrip() for line in infile))
    infile.close()    

    # Write output .gro file
    # tacount is the total atom count
    outfile = open(outname, 'w')
    outfile.write(lines[0])
    outfile.write('\n')
    outfile.write(lines[1])
    outfile.write('\n')
    tacount = 0
    for mol in range(len(new_xyz_position[:,0,0])):
        for atom in range(len(new_xyz_position[0,0,:])):
            # mcount is the current molecule tracker
            mcount = mol + 1
            tacount += 1
            x = round(new_xyz_position[mol,0,atom],8)
            y = round(new_xyz_position[mol,1,atom],8)
            z = round(new_xyz_position[mol,2,atom],8)
    	    if x < 0.0:
    	        xstr = '   ' + "%.8f" % x
    	    elif x == 0.0:
    	        xstr = '    0.00000000'
    	    else:
    	        xstr = '    ' + "%.8f" % x
    	    if y < 0.0:
                    ystr = '   ' + "%.8f" % y
    	    elif y == 0.0:
    	        ystr = '    0.00000000'
            else:
                    ystr = '    ' + "%.8f" % y
    	    if z < 0.0:
                zstr = '   ' + "%.8f" % z
    	    elif z == 0.0:
    	        zstr = '    0.00000000'
            else:
                zstr = '    ' + "%.8f" % z
    
            #xstr = '%(x) 10g' %vars()
            #ystr = '%(y) 10g' %vars()
            #zstr = '%(z) 10g' %vars()
            line = str(mcount).rjust(5) + mol_name + lines[tacount+1].split()[1].rjust(7) + str(tacount).rjust(5) + xstr + ystr + zstr + '\n'
            outfile.write(line)
    outfile.write("   " + '   '.join(newboxvect))
    outfile.write('\n')  # make sure the file ends on a newline
    outfile.close()


changeBoxVector(fname,molecule,volume=final_volume,boxvect=boxvect)
