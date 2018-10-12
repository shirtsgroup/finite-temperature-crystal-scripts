# /usr/bin/python

# PYTHON SCRIPT TO ADD FULL PRECISION TO A TINKER XYZ FILE GIVEN A GRO FILE
# param x - The name of the .xyz file
# param n - The scaling factor to apply to each xyz coordinate

import numpy as np
from optparse import OptionParser


def lattice_parameters_from_gro(grofile):
    gro_values = np.array(open(grofile, "r").readlines()[-1].strip("\n").split()).astype(float)
    crystal_tensor = np.zeros((3, 3))
    crystal_tensor[0, 0] = gro_values[0]
    crystal_tensor[1, 1] = gro_values[1]
    crystal_tensor[2, 2] = gro_values[2]
    if len(gro_values) > 3:
        crystal_tensor[0, 1] = gro_values[5]
        crystal_tensor[0, 2] = gro_values[7]
        crystal_tensor[1, 2] = gro_values[8]
    crystal_tensor = 10.*crystal_tensor

    # Computing lattice parameters
    a = np.linalg.norm(crystal_tensor[:, 0])
    b = np.linalg.norm(crystal_tensor[:, 1])
    c = np.linalg.norm(crystal_tensor[:, 2])

    gamma = np.arccos(np.dot(np.squeeze(np.asarray(crystal_tensor[:, 0])), np.squeeze(np.asarray(crystal_tensor[:, 1])))
                      / (a * b)) * 180. / np.pi
    alpha = np.arccos(np.dot(np.squeeze(np.asarray(crystal_tensor[:, 1])), np.squeeze(np.asarray(crystal_tensor[:, 2])))
                      / (b * c)) * 180. / np.pi
    beta = np.arccos(np.dot(np.squeeze(np.asarray(crystal_tensor[:, 2])), np.squeeze(np.asarray(crystal_tensor[:, 0])))
                     / (c * a)) * 180. / np.pi

    # Creating an array of lattice parameters
    lattice_parameters = "   " + str(a) + "   " + str(b) + "   " + str(c) + "   " + str(alpha) + "   " + str(beta) + \
						 "   " + str(gamma) + "\n"
    return lattice_parameters


# =============================================================================================
# READ IN USER INPUTS
# =============================================================================================
parser = OptionParser()
parser.add_option('-x', dest='xyzfile', help='Input xyz file', default='frame.xyz')
parser.add_option('-n', dest='num', help='scaling factor', default='10')
(options, args) = parser.parse_args()
scaling_factor = float(options.num)
# Create the space buffer
spaces = []
spaces.append('     ')
spaces.append('  ')
spaces.append('     ')
spaces.append('    ')
spaces.append('    ')
spaces.append('     ')
spaces.append('    ')
spaces.append('    ')
spaces.append('    ')
rjustvect = np.ones(10, int)
rjustvect = ['6', '3', '15', '15', '15', '5', '6', '6', '6']


fname = options.xyzfile

infile = open(fname, 'r')
lines = infile.readlines()
infile.close()

outfile = open(fname, 'w')

linecount = 1
for line in lines:
    if linecount == 1:
        outfile.write(line)
        if options.grofile != '':
            outfile.write(lattice_parameters_from_gro(options.grofile))
    else:
        tokens = line.split()
        linestring = ""
        for i, token in enumerate(tokens):
            if i == 2 or i == 3 or i == 4:
                token = round(float(token) * scaling_factor, 12)
            linestring = linestring + str(token).rjust(int(rjustvect[i]))
        outfile.write(linestring + '\n')
    linecount = linecount + 1
    if linecount == 11:
        spaces[0] = '    '
    elif linecount == 101:
        spaces[0] = '   '
    elif linecount == 1001:
        spaces[0] = '  '

outfile.close()
