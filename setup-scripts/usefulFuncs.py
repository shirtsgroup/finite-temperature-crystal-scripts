#!/usr/bin/python
#
# Collection of useful functions for python 
#
# Copyright Eric C. Dybeck, University of Virginia, 2016
#

import numpy
import pdb

#=============================================================================================
# OUTPUT VARIABLES
#=============================================================================================
#f - Array of function values at each point f(a) -> f(b)
#dx - step size between x points
def simpson(f, dx):
    n=len(f)-1
    simp=f[0]+f[n]
    
    #Ensure that the length of f is even
    if n % 2 == 1:
	print "f must be even for simpsons rule!"
	sys.exit()

    for i in range(1,n,2):
	simp += 4*f[i]
    for i in range(2,n-1,2):
	simp += 2*f[i]

    return simp*dx/3.0

#x=numpy.arange(0,11)
#x=x/10.0
#y=[a*a*a*a*a for a in x]
#
#print x
#print y
#area=simpson(y,0.1)
#print area
