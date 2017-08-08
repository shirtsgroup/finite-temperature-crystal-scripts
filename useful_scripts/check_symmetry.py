#!/bin/python

#
# Replicate the unit cell of a benzene polymorph into a supercell and expand its contents if desired
#
# Copyright Eric Dybeck and Michael R. Shirts, University of Virginia, 2014
#
import numpy # numerical array library
import pymbar # multistate Bennett acceptance ratio
from pymbar import timeseries # timeseries analysis
from optparse import OptionParser # for parsing command-line options
import pdb

parser = OptionParser()
parser.add_option('-g', dest = 'gro', help = 'Gro file to analyze', default = 'benzene_PROD.gro')
parser.add_option('-l', dest = 'length', help = 'Supercell length', default = 3)
parser.add_option('-w', dest = 'width', help = 'Supercell width', default = 2)
parser.add_option('-z', dest = 'height', help = 'Supercell height', default = 3)
parser.add_option('-m', dest = 'apermol', help = 'Number of atoms per molecule or molecule group', default = 12)
parser.add_option('-n', dest = 'index', help = 'Index file', default = 'benzene_index.ndx')
parser.add_option('-u', dest = 'unit', help = 'Number of benzene molecules in the unit cell', default = 4)
nndx = 15 # number of ndx entries per line
#Ensure that user inputs are sensible
(options, args) = parser.parse_args()
apermol = options.apermol
Length = int(options.length)
Width = int(options.width)
Height = int(options.height)
Molecule_number = int(options.unit)

#Determine the total number of independent atoms
independent_atoms = Molecule_number * apermol

#Read in the gro file
fname = options.gro
infile = open(fname, 'r')
lines = infile.readlines()
infile.close()
print "loading " + fname

#Read in the index file
fname_index = options.index
infile = open(fname_index, 'r')
lines_index = infile.readlines()
infile.close()
print "loading " + fname_index

#Determine the unit dimensions from the last line of the .gro file
v1 = numpy.zeros(3,float)
v2 = numpy.zeros(3,float)
v3 = numpy.zeros(3,float)
tokens = lines[len(lines)-1].split()
stopper = ""	#Need to use this later during the loop
if len(tokens)==3:
    stopper = tokens[0]
    v1[0] = float(tokens[0])/Length
    v2[1] = float(tokens[1])/Width
    v3[2] = float(tokens[2])/Height
    boxvect_length=3
else:
    stopper = tokens[0]
    v1[0] = float(tokens[0])/Length
    v2[1] = float(tokens[1])/Width
    v3[2] = float(tokens[2])/Height
    v1[1] = float(tokens[3])/Length
    v1[2] = float(tokens[4])/Length
    v2[0] = float(tokens[5])/Width
    v2[2] = float(tokens[6])/Width
    v3[0] = float(tokens[7])/Height
    v3[1] = float(tokens[8])/Height
    boxvect_length=9


maxdeviation = 0.0;

#Loop over all the independent atoms
atom_counter = 0
for a in range(independent_atoms):
    tokens = lines[a+2].split() #Read in the atomic coordinates in the gro file

    #Find this atom's group in the index file
    for l,line in enumerate(lines_index):
	if l <= 5:
	    continue
	tokens_index = line.split()
	if tokens_index[0] == str(a+1) and (lines_index[l-1][2] == "C" or lines_index[l-1][2] == "H"):
    	    #Read in the vector of identical atoms
	    index_vect = line.split()

    #Collapse the initial atom down to the 1x1x1 unit cell.
    for i in range(len(index_vect)):
	tokens = lines[int(index_vect[i])+1].split()
	if i == 0: #Determine the initial xyz offset from (0,0,0)
	    xoffset = float(tokens[3])
	    yoffset = float(tokens[4])
	    zoffset = float(tokens[5])
	    #print xoffset
	    #print v1
	    #print v2
	    #print v3
	    #pdb.set_trace()
	    while xoffset > (v1[0] + v2[0] + v3[0]):
		xoffset -= (v1[0] + v2[0] + v3[0])
	    while xoffset < 0:
		xoffset += (v1[0] + v2[0] + v3[0])
	    while yoffset > (v1[1] + v2[1] + v3[1]):
		yoffset -= (v1[1] + v2 [1] + v3[1])
	    while yoffset < 0:
		yoffset += (v1[1] + v2 [1] + v3[1])
	    while zoffset > (v1[2] + v2[2] + v3[2]):
		zoffset -= (v1[2] + v2[2] + v3[2])
	    while zoffset < 0:
		zoffset += (v1[2] + v2[2] + v3[2])
	else:
	    tempxoffset = float(tokens[3])
            tempyoffset = float(tokens[4])
            tempzoffset = float(tokens[5])
            while tempxoffset > (v1[0] + v2[0] + v3[0]):
                tempxoffset -= (v1[0] + v2[0] + v3[0])
            while tempxoffset < 0:
                tempxoffset += (v1[0] + v2[0] + v3[0])
            while tempyoffset > (v1[1] + v2[1] + v3[1]):
                tempyoffset -= (v1[1] + v2[1] + v3[1])
            while tempyoffset < 0:
                tempyoffset += (v1[1] + v2[1] + v3[1])
            while tempzoffset > (v1[2] + v2[2] + v3[2]):
                tempzoffset -= (v1[2] + v2[2] + v3[2])
            while tempzoffset < 0:
                tempzoffset += (v1[2] + v2[2] + v3[2])

	    tempxoffset2 = tempxoffset + (v1[0] + v2[0] + v3[0])
	    tempxoffset3 = tempxoffset - (v1[0] + v2[0] + v3[0])
	    tempyoffset2 = tempyoffset + (v1[1] + v2[1] + v3[1])
            tempyoffset3 = tempyoffset - (v1[1] + v2[1] + v3[1])
	    tempzoffset2 = tempzoffset + (v1[2] + v2[2] + v3[2])
            tempzoffset3 = tempzoffset - (v1[2] + v2[2] + v3[2])
	    



	    if abs(tempxoffset2-xoffset) < abs(tempxoffset - xoffset):
		tempxoffset = tempxoffset2
	    if abs(tempxoffset3-xoffset) < abs(tempxoffset - xoffset):
		tempxoffset = tempxoffset3
	    if abs(tempyoffset2-yoffset) < abs(tempyoffset - yoffset):
                tempyoffset = tempyoffset2
            if abs(tempyoffset3-yoffset) < abs(tempyoffset - yoffset):
                tempyoffset = tempyoffset3
	    if abs(tempzoffset2-zoffset) < abs(tempzoffset - zoffset):
                tempzoffset = tempzoffset2
            if abs(tempzoffset3-zoffset) < abs(tempzoffset - zoffset):
                tempzoffset = tempzoffset3

	    if tempxoffset != xoffset or tempyoffset != yoffset or tempzoffset != zoffset:
		if ((tempxoffset - xoffset)**2)**0.5 > maxdeviation:
		    maxdeviation = ((tempxoffset - xoffset)**2)**0.5
		if ((tempyoffset - yoffset)**2)**0.5 > maxdeviation:
                    maxdeviation = ((tempyoffset - yoffset)**2)**0.5
		if ((tempzoffset - zoffset)**2)**0.5 > maxdeviation:
                    maxdeviation = ((tempzoffset - zoffset)**2)**0.5
		#if maxdeviation > 0.1:
		#    print "Max Deviation: " + str(maxdeviation)
		#    sys.exit()
		if (((tempxoffset - xoffset)**2)**0.5) > 0.0001 or (((tempyoffset - yoffset)**2)**0.5) > 0.0001 or (((tempzoffset - zoffset)**2)**0.5) > 0.0001:
		    print "Unsymetric Atoms!"
                    print str(index_vect[0]) + ": " + str(xoffset) + " " + str(yoffset) + " " + str(zoffset)
                    print str(index_vect[i]) + ": " + str(tempxoffset) + " " + str(tempyoffset) + " " + str(tempzoffset)
                    print lines[int(index_vect[0])+1]
                    print lines[int(index_vect[i])+1]
                    print str(v1[0]) + " " + str(v2[1]) + " " + str(v3[2])
                    print index_vect
                    print "DEVIATION: " + str(((tempxoffset - xoffset)**2)**0.5)
                    print "DEVIATION: " + str(((tempyoffset - yoffset)**2)**0.5)
                    print "DEVIATION: " + str(((tempzoffset - zoffset)**2)**0.5)

print "MAX DEVIATION: " + str(maxdeviation)

"""
#Open the output gro file
numout = molecule_number*Length*Width*Height
fname = "benzene_" + poly + "_" + str(numout) + ".gro"
open(fname, 'w').close() 	#delete any previous contents
outfile = open(fname, 'w')

#Open the output index file
fname_index = "benzene_" + poly + "_" + str(numout) + ".ndx"
open(fname_index, 'w').close()        #delete any previous contents
outfile_index = open(fname_index, 'w')

#Initialize the atom header container and the dictionary of atom numbers for each header
index_headers = [];
atom_indicies = dict();
atom_types = dict();

atom_types['C']=0
atom_types['H']=0

atom_count = 0
total_atom_count = 0
linecount = 0
moleculecount = 1
for i in range(0,Length):
    for j in range(0,Width):
	for k in range(0,Height):
	    for l,line in enumerate(lines):
		tokens = line.split()
		linecount+=1
		if tokens[0] == "benzene" or tokens[0] == "48" or tokens[0] == "24" or tokens[0] == stopper:
        	    continue
    		if l == 0 or l == 1:
		    continue
		if (float(tokens[3])+v1[0]*i+v2[0]*j+v3[0]*k) < 0.0:
		    tokens[3] = "  " + str("%.3f" % (float(tokens[3])+v1[0]*i+v2[0]*j+v3[0]*k))
		else:
		    tokens[3] = "   " + str("%.3f" % (float(tokens[3])+v1[0]*i+v2[0]*j+v3[0]*k))
		if (float(tokens[4])+v1[1]*i+v2[1]*j+v3[1]*k) < 0.0:
                    tokens[4] = "  " + str("%.3f" % (float(tokens[4])+v1[1]*i+v2[1]*j+v3[1]*k))
                else:
                    tokens[4] = "   " + str("%.3f" % (float(tokens[4])+v1[1]*i+v2[1]*j+v3[1]*k))
		if (float(tokens[5])+v1[2]*i+v2[2]*j+v3[2]*k) < 0.0:
                    tokens[5] = "  " + str("%.3f" % (float(tokens[5])+v1[2]*i+v2[2]*j+v3[2]*k))
                else:
                    tokens[5] = "   " + str("%.3f" % (float(tokens[5])+v1[2]*i+v2[2]*j+v3[2]*k))

		atom_count+=1
		total_atom_count+=1 
		if atom_count > apermol:	#We have reached another benzene molecule
		    atom_count = 1;
		    moleculecount+=1

		#Now format the outputs
		
		#If this is the first time encountering this atom, create a header for the index file
                if i==0 and j==0 and k==0:
                    #header = '[ ' + tokens[1] + '    ' + str((int(atom_types[tokens[1]]) + 1)) + ' ]'
                    header = '[ ' + tokens[1] + str((int(atom_types[tokens[1]]) + 1)) + ' ]'
                    atom_types[tokens[1]] = atom_types[tokens[1]] + 1
                    index_headers.append(header)
                    atom_indicies[header] = str(total_atom_count)
                else:
                    atom_indicies[index_headers[l-2]] = atom_indicies[index_headers[l-2]] + '    ' + str(total_atom_count)

		if moleculecount < 10:
		    tokens[0] = "    " + str(moleculecount) + "BNZ"
		elif moleculecount < 100:
		    tokens[0] = "   " + str(moleculecount) + "BNZ"
		elif moleculecount < 1000:
                    tokens[0] = "  " + str(moleculecount) + "BNZ"
		else:
		    tokens[0] = " " + str(moleculecount) + "BNZ"
		
		tokens[1] = "      " + tokens[1]

		if total_atom_count < 10:
                    tokens[2] = "    " + str(total_atom_count)
                elif total_atom_count < 100:
                    tokens[2] = "   " + str(total_atom_count)
                elif total_atom_count < 1000:
                    tokens[2] = "  " + str(total_atom_count)
                else:
                    tokens[2] = " " + str(total_atom_count)
		
		outfile.write(''.join(tokens))
		outfile.write('\n')

#write the last line to the output file
tokens = lines[len(lines)-1].split()
if len(tokens)==3:
    tokens[0] = str(v1[0]*Length)
    tokens[1] = str(v2[1]*Width)
    tokens[2] = str(v3[2]*Height)
else:
    tokens[0] = str(v1[0]*Length)
    tokens[1] = str(v2[1]*Width)
    tokens[2] = str(v3[2]*Height)
    tokens[3] = str(v1[1]*Length)
    tokens[4] = str(v1[2]*Length)
    tokens[5] = str(v2[0]*Width)
    tokens[6] = str(v2[2]*Width)
    tokens[7] = str(v3[0]*Height)
    tokens[8] = str(v3[1]*Height)
outfile.write("   " + '   '.join(tokens))
outfile.write('\n')
outfile.close()
print 'wrote final gro file to: ' + fname

#Write to the index file
# write [ system ] group 
outfile_index.write('[ System ]\n')
for i in range(1,numout*apermol+1):
    outfile_index.write("%5d" % (i))
    if (i % nndx) == 0:
        outfile_index.write('\n')
outfile_index.write('\n')
for header in index_headers:
    outfile_index.write(header)
    outfile_index.write('\n')
    outfile_index.write(atom_indicies[header])
    outfile_index.write('\n')
outfile_index.close()
print 'wrote final ndx file to: ' + fname_index
print 'put this section in the mdp:'
print '-------------------'
print 'symmetry-averaging = yes'
print 'ave-ngroups = %d' % (len(index_headers))
for i,header in enumerate(index_headers):
    grpname = header[2:-2]
    print "ave-group%d = %s" % (i,grpname)
print '-------------------'
"""
