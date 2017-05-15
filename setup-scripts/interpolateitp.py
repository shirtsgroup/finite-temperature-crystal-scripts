#!/bin/python
#
# Interpolate intramolecular parameters in an itp file 
#
# Copyright Eric Dybeck and Michael R. Shirts, University of Virginia, 2014
#
import numpy # numerical array library
from optparse import OptionParser # for parsing command-line options
import os
import pdb

parser = OptionParser()
parser.add_option('-f', dest = 'itpfile', help = 'itpfile to interpolate', default = 'itpfile')
parser.add_option('-d', dest = 'interpolation', help = 'Percentage to interpolate between state A and state B', default = 0.0)

#Ensure that user inputs are sensible
(options, args) = parser.parse_args()
fname = options.itpfile
atom = options.atomname

#Read in the original gromacs file
infile = open(fname, 'r')
lines = infile.readlines()
infile.close()
print "loading " + fname

#Write the first two lines to the new gro file
fname = fname[:len(fname)-4] + "_drude.gro"
outfile = open(fname, 'w')
print "writing to " + fname
outfile.write(lines[0])
tokens = lines[1].split()
outfile.write(" " + str(int(1.5*int(tokens[0]))))
outfile.write('\n')

#Add a drude particle to all specified atoms
atom_counter=1
for i,line in enumerate(lines):
    if i<2 or i==len(lines)-1:
	continue

    tokens=line.split()

    #Add the original atom (with the updated atom number)
    if (float(tokens[3])) < 0.0:
	tokens[3] = "  " + str("%.8f" % (float(tokens[3])))
    else:
	tokens[3] = "   " + str("%.8f" % (float(tokens[3])))
    if (float(tokens[4])) < 0.0:
	tokens[4] = "  " + str("%.8f" % (float(tokens[4])))
    else:
	tokens[4] = "   " + str("%.8f" % (float(tokens[4])))
    if (float(tokens[5])) < 0.0:
	tokens[5] = "  " + str("%.8f" % (float(tokens[5])))
    else:
	tokens[5] = "   " + str("%.8f" % (float(tokens[5])))

    moleculecount = int(tokens[0][:len(tokens[0])-3])
    if moleculecount < 10:
	tokens[0] = "    " + str(moleculecount) + "BNZ"
    elif moleculecount < 100:
	tokens[0] = "   " + str(moleculecount) + "BNZ"
    elif moleculecount < 1000:
	tokens[0] = "  " + str(moleculecount) + "BNZ"
    else:
	tokens[0] = " " + str(moleculecount) + "BNZ"
    
    tokens[1] = "      " + tokens[1]

    if atom_counter < 10:
	tokens[2] = "    " + str(atom_counter)
    elif atom_counter < 100:
	tokens[2] = "   " + str(atom_counter)
    elif atom_counter < 1000:
	tokens[2] = "  " + str(atom_counter)
    else:
	tokens[2] = " " + str(atom_counter)

    outfile.write(''.join(tokens))
    outfile.write('\n')
	
    atom_counter+=1

    #Add a drude oscillator particle if this is one of the specified atoms
    tokens[1] = line.split()[1]
    if tokens[1] == atom:
	tokens[1] = "     D" + tokens[1]

	if atom_counter < 10:
	    tokens[2] = "    " + str(atom_counter)
	elif atom_counter < 100:
	    tokens[2] = "   " + str(atom_counter)
	elif atom_counter < 1000:
	    tokens[2] = "  " + str(atom_counter)
	else:
	    tokens[2] = " " + str(atom_counter)

	outfile.write(''.join(tokens))
	outfile.write('\n')

	atom_counter+=1

#Add the last line of the gro file
outfile.write(lines[len(lines)-1])

#Finally, change the second line of the gro file to update the true number of particles in the system
#origatoms = int(lines[1].split()[0])
#sedcommand="sed -i " + "'0,/" + str(origatoms) + "/{s/" + str(origatoms) + "/" + str(atom_counter-1)+ "/}' " + fname
#os.system(sedcommand)



	


