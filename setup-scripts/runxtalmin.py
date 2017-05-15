#!/home/nasc5274/anaconda2/bin/python

import numpy
import xtalmin
import resize_gro
import pdb

fname="MIN.gro"

#Vector="1.657047390000   2.208901680000   2.197299030000 -0.775322190000"

#Vector=numpy.array([float(j) for j in Vector.split()])

#energy = xtalmin.potenergy(Vector)

#print energy


xtalmin.xtalmin(fname)

