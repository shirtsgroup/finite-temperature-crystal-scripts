#/usr/bin/python
#PYTHON SCRIPT TO ADD FULL PRECISION TO A GRO FILE GIVEN A TINKER XYZ FILE
#param x - The name of the .gro file
#param n - The scaling factor to apply to each xyz coordinate

import numpy # numerical array library
from optparse import OptionParser # for parsing command-line options
import pdb

#=============================================================================================
# READ IN USER INPUTS
#=============================================================================================
parser = OptionParser()
parser.add_option('-g', dest = 'grofile', help = 'Input gro file', default='frame.gro')
parser.add_option('-n', dest = 'num' , help = 'scaling factor', default='0.1')
(options, args) = parser.parse_args()
scaling_factor = float(options.num)
#Create the space buffer
spaces = [];
spaces.append('     ')
spaces.append('  ')
spaces.append('     ')
spaces.append('    ')
spaces.append('    ')
spaces.append('     ')
spaces.append('    ')
spaces.append('    ')
spaces.append('    ')
rjustvect = numpy.ones(10,int);
rjustvect = ['8', '7', '5', '15', '15', '15', '6', '6', '6']


fname = options.grofile

infile = open(fname, 'r')
lines = infile.readlines()
infile.close()

outfile = open(fname, 'w')

linecount=1
for line in lines:
    if linecount<3 or linecount==len(lines):
	outfile.write('    ' + line)
    else:
	tokens = line.split()
	linestring=""
	for i,token in enumerate(tokens):
	    if i==3 or i==4 or i==5:
		token = str("%.8f" % round(float(token)*scaling_factor,8))
	    linestring = linestring +  str(token).rjust(int(rjustvect[i]))
	outfile.write(linestring + '\n')
    linecount=linecount+1
outfile.close()
