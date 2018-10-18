#!/usr/bin/env python

from __future__ import print_function
import os
import sys
import yaml
import numpy as np

current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, current_path + '/analysis-scripts')
from dA_Lambda_MBAR import dA_Lambda_MBAR
from dA_Gamma_MBAR import dA_Gamma_MBAR
from dGvsT import dGvsT

def yaml_loader(file_path):
    # Loads in a ymal file
    with open(file_path, "r") as input_file:
        data = yaml.load(input_file)
    return data

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Sets up directories to run MD simulations for the free energy '
                                                 'differences between organic crystals')
    parser.add_argument('-i', '--input_file', dest='input_file', default='input.yaml',
                        help='Input file containing all parameters to set up the directories to run the MD simulations')

    args = parser.parse_args()
    inputs = yaml_loader(args.input_file)

    if inputs['gen_in']['independent'] == 'same':
        independent = inputs['gen_in']['number_of_molecules']
    else:
        independent = inputs['gen_in']['independent']

    dA_L, ddA_L = dA_Lambda_MBAR(MinL=inputs['PSCP_in']['min_lambda'], MaxL=inputs['PSCP_in']['max_lambda'],
                                 dL=inputs['PSCP_in']['lambda_spacing'], GAMMA=inputs['PSCP_in']['gamma'],
                                 exponent=inputs['PSCP_in']['lambda_exponent'], Molecule=inputs['gen_in']['molecule'],
                                 polymorphs=inputs['gen_in']['polymorph_num'],
                                 Molecules=inputs['gen_in']['number_of_molecules'], Independent=independent,
                                Temp=inputs['PSCP_in']['PSCP_temperature'], Pressure=inputs['gen_in']['pressure'],
                                 k=inputs['PSCP_in']['k_max'], potential=inputs['gen_in']['potential'],
                                 hinge=inputs['gen_in']['hinge'])

    dA_G, ddA_G = dA_Gamma_MBAR(MINGAMMA=inputs['PSCP_in']['min_gamma'], MAXGAMMA=inputs['PSCP_in']['max_gamma'],
                                GSPACING=inputs['PSCP_in']['gamma_spacing'], LAMBDA=inputs['PSCP_in']['lambda'],
                                exponent=inputs['PSCP_in']['gamma_exponent'],
                                polymorphs=inputs['gen_in']['polymorph_num'], Molecule=inputs['gen_in']['molecule'],
                                Molecules=inputs['gen_in']['number_of_molecules'], Independent=independent,
                                Temp=inputs['PSCP_in']['PSCP_temperature'],
                                Pressure=inputs['gen_in']['pressure'], k=inputs['PSCP_in']['k_max'],
                                potential=inputs['gen_in']['potential'], hinge=inputs['gen_in']['hinge'])

    refdG = (dA_L - dA_L[0]) + (dA_G - dA_G[0])
    #NSA: This is wrong!!!!
    refddG = ddA_L + ddA_G
    absolutedU = np.array(inputs['analysis_in']['potential_energy'].split()).astype(float)
    refdU = absolutedU - absolutedU[0]

    dGvsT(Temperatures=np.array(inputs['temp_in']['temperatures'].split()).astype(int), Pressure=inputs['gen_in']['pressure'],
          Molecules=inputs['gen_in']['number_of_molecules'], molecule=inputs['gen_in']['molecule'],
          Independent=independent, potential=inputs['gen_in']['potential'],
          simulation=inputs['temp_in']['simulation_package'], hinge=inputs['gen_in']['hinge'],
          Polymorphs=inputs['gen_in']['polymorph_num'].split(), refT=inputs['PSCP_in']['PSCP_temperature'],
          refdG=refdG, refddG=refddG, refdU=refdU, absolutedU=absolutedU)




