from __future__ import print_function
# !/usr/bin/python
#
# Interpolate an itp file
# All interpolation lines are headed with the keyword INTERPOLATE followed by the number of parameters that will be
#    interpolated.
# This line is then followed by the (0.0) and (1.0) values for each parameter.
#
# Copyright Michael R. Shirts, University of Virginia, 2014
#

import numpy as np
from optparse import OptionParser
import os
import sys

# =============================================================================================
# INPUT VARIABLES
# =============================================================================================
# ITP file to be interpolated
# Fraction of the interpolation to take place


# =============================================================================================
# OUTPUT VARIABLES
# =============================================================================================
# Final ITP file with interpolated parameters

parser = OptionParser()
parser.add_option('-f', '--file', dest='file', help='Initial itp file', default='benzene.itp')
parser.add_option('-d', '--frac', dest='frac', help='Fraction of interpolation', default=0.0)
(options, args) = parser.parse_args()

fname = options.file
frac = float(options.frac)

# =============================================================================================
# Ensure that all inputs are correct
# =============================================================================================
if frac < 0.0 or frac > 1.0:
    print("Fraction (" + str(frac) + ") must be between 0.0 and 1.0!")
    sys.exit()

# =============================================================================================
# Loop through all instances of 'INTERPOLATE' and interpolate the parameters
# =============================================================================================

infile = open(fname, 'r')
lines = infile.readlines()
infile.close()

print("INTERPOLATING...")
for line in lines:
    tokens = line.split()
    if 'INTERPOLATE' not in tokens:
        continue
    
    for i, token in enumerate(tokens):
        if token != 'INTERPOLATE':
            continue
        # Determine the number of parameters to interpolate
        numParams = int(tokens[i + 1])
        Params = np.zeros([numParams, 3], float)  # Row 1 is state A, Row 2 is state B, Row 3 is the interpolated state
        Params[:, 0] = [float(j) for j in tokens[i + 2: i + 2 + numParams]]  # Grab the State A Parameters
        Params[:, 1] = [float(j) for j in tokens[i + 2 + numParams: i + 2 + 2 * numParams]]  # Grab the State B Parameters
        Params[:, 2] = (1.0 - frac) * Params[:, 0] + frac * Params[:, 1]  # Interpolate to determine the new parameters
        Params[:, 2] = [round(j, 6) for j in Params[:, 2]]
        StrParams= [str(j) for j in Params[:, 2]]

        # Use sed to replace the line with the new parameters
        sedcommand = "sed -i \'0,/INTERPOLATE.*/{s/INTERPOLATE.*/" + '   '.join(StrParams[:]) + "/}\' " + fname
        os.system(sedcommand)
