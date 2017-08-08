#/usr/bin/python
#PYTHON SCRIPT TO TAKE IN AN XVG FILE AND CREATE THE SAME DATA IN DIFFERENT UNITS
import numpy # numerical array library
import pymbar # multistate Bennett acceptance ratio
from pymbar import timeseries # timeseries analysis
from optparse import OptionParser # for parsing command-line options
import pdb

#=============================================================================================
# READ IN USER INPUTS
#=============================================================================================
parser = OptionParser()
parser.add_option('-f', dest = 'infile', help = 'Input xvg file')
parser.add_option('-s', dest = 'scalefactor', help = 'factor to scale all energies by', default=1.0)
parser.add_option('-i', dest = 'ignore', help = 'column to ignore from the conversion', default='5')
parser.add_option('-H','--hinge', dest = 'hinge', help = 'Optional hinge to new file', default='converted')
(options, args) = parser.parse_args()
scale = float(options.scalefactor)
hinge = '_' + str(options.hinge) + '.xvg'

ignorecols=[];
for i,token in enumerate(options.ignore.split()):
    ignorecols.append(token)

ignore_symbols = ['#', '@', '@TYPE', 'STEP', '=====================', 'TIME(ps)'];
fname = str(options.infile)
infile = open(fname,'r')
lines = infile.readlines()
infile.close()
print "loading " + str(fname)
fname = str(options.infile)[:len(options.infile)-4]+hinge
outfile = open(fname, 'w')

counter = 1
for line in lines:
    tokens = line.split()
    if tokens[0] in ignore_symbols:
	outfile.write(line)
	continue
    for i,token in enumerate(tokens):
	if i==0 or str(i) in ignorecols:
	    continue
	else:
    	    tokens[i] = str(scale*float(token))
    outfile.write('   '.join(tokens) + '\n')
    counter+=1
