#/usr/bin/python
#PYTHON SCRIPT TO ADD FULL PRECISION TO A TINKER XYZ FILE GIVEN A GRO FILE
#param x - The name of the .xyz file
#param n - The scaling factor to apply to each xyz coordinate

import numpy # numerical array library
from optparse import OptionParser # for parsing command-line options
import pdb

#=============================================================================================
# READ IN USER INPUTS
#=============================================================================================
parser = OptionParser()
parser.add_option('-x', dest = 'xyzfile', help = 'Input xyz file', default='frame.xyz')
parser.add_option('-n', dest = 'num' , help = 'scaling factor', default='10')
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
rjustvect = ['6', '3', '15', '15', '15', '5', '6', '6', '6']


fname = options.xyzfile

infile = open(fname, 'r')
lines = infile.readlines()
infile.close()

outfile = open(fname, 'w')

linecount=1
for line in lines:
    if linecount==1:
	outfile.write(line)
    else:
	tokens = line.split()
	linestring=""
	for i,token in enumerate(tokens):
	    if i==2 or i==3 or i==4:
		token = round(float(token)*scaling_factor,8)
	    #linestring = linestring + spaces[i] + str(token)
	    linestring = linestring +  str(token).rjust(int(rjustvect[i]))
	outfile.write(linestring + '\n')
        #outfile.write('     ' + tokens[0] + '  ' + tokens[1] + '      ' + tokens[2] + '    ' + tokens[3] + '    ' + tokens[4] + '     ' + tokens[5] + '\n')
    linecount=linecount+1
    if linecount==11:
	spaces[0]='    '
    elif linecount==101:
	spaces[0]='   '
    elif linecount==1001:
	spaces[0]='  '
outfile.close()
