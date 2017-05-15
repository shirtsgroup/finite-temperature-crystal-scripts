#/usr/bin/python
#PYTHON SCRIPT TO TAKE IN AN XVG FILE OF BOX VECTORS AT EACH STEP AND RETURN AN AVERAGE
#NOTE: The Box ZX number is taken as the absolute value due to flipping in Benzene 3
import numpy # numerical array library
from optparse import OptionParser # for parsing command-line options
import os
import pdb

#=============================================================================================
# READ IN USER INPUTS
#=============================================================================================
parser = OptionParser()
parser.add_option('-f', dest = 'infile', help = 'Input xvg file')
parser.add_option('-g', dest = 'defaultgro', help = 'Default box vector used to match the sign', default='none')
(options, args) = parser.parse_args()

ignore_symbols = ['#', '@', '@TYPE', 'STEP', '=====================', 'TIME(ps)'];
defaultvect=[]
if options.defaultgro != 'none':
    fname = str(options.defaultgro)
    infile = open(fname,'r')
    lines = infile.readlines()
    infile.close()
    defaultvect = lines[len(lines)-1].split()
    defaultvect = [float(j) for j in defaultvect]
    for i,X in enumerate(defaultvect):
	if X!=0:
	    defaultvect[i] = numpy.round(numpy.abs(X)/X,0)
    #Remove the extraneous terms
    if len(defaultvect)==9:
	defaultvect=[defaultvect[j] for j in [0,1,2,5,7,8]]
	print defaultvect   

fname = str(options.infile)
infile = open(fname,'r')
lines = infile.readlines()
infile.close()
print "loading " + str(fname)
vectors=numpy.zeros([6,100000])
tol=0.01
counter = 0
for line in lines:
    tokens = line.split()
    if tokens[0] in ignore_symbols:
	continue
    boxvect=[float(j) for j in tokens[1:]] 
    #If the value is close to zero, set it to exactly zero
    for i,X in enumerate(boxvect):
	if numpy.abs(X) < tol:
	    boxvect[i]=0.0000
    #If a default vector was given, ensure that the sign is the same
    if defaultvect!=[]:
	boxvect = [numpy.abs(j) for j in boxvect]
	for i in range(len(defaultvect)):
	    boxvect[i] = boxvect[i]*defaultvect[i] 
    vectors[:len(boxvect),counter]=boxvect
    counter+=1

averagevect=numpy.zeros(len(boxvect))
for i in range(len(boxvect)):
    averagevect[i] = numpy.average(vectors[i,:counter])
averagevect = [str(j) for j in averagevect]
averagevectstr='   '.join(averagevect)
os.system('echo "' +  averagevectstr + '" > avgboxvect.xvg')



