#!/usr/bin/env python

from __future__ import print_function
import os
import sys
import yaml
import numpy as np
import Run_setup as RS

current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, current_path + '/analysis-scripts')
from dA_Lambda_MBAR import dA_Lambda_MBAR
from dA_Gamma_MBAR import dA_Gamma_MBAR
from dGvsT import dGvsT
from dGvsT import old_systems_dictionary

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Runs analysis of a PSCP and ')
    parser.add_argument('-i', '--input_file', dest='input_file', default='input.yaml',
                        help='Input file containing all parameters to set up the directories to run the MD simulations')

    args = parser.parse_args()
    inputs = RS.yaml_loader(args.input_file)

    if inputs['gen_in']['independent'] == 'same':
        independent = inputs['gen_in']['number_of_molecules']
    else:
        independent = inputs['gen_in']['independent']

    if (inputs['PSCP_out']['refdG'] == None) or (inputs['PSCP_out']['refddG'] == None):
        if (inputs['PSCP_in']['run_restraints'] == True) and (inputs['PSCP_in']['run_interactions'] == True):
            # Computing the free energy difference from the PSCP
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

            inputs['PSCP_out']['refdG'] = (dA_L - dA_L[0]) + (dA_G - dA_G[0])
            inputs['PSCP_out']['refddG'] = np.sqrt(ddA_L**2 + ddA_G**2)
#            absolutedU = np.array(inputs['PSCP_out']['absolutedU']).astype(float)
#            refdU = absolutedU - absolutedU[0]
        else:
            # Pulling old data from Temperature Transformation paper 2017
            inputs['gen_in']['polymorph_num'], inputs['PSCP_in']['PSCP_temperature'], inputs['PSCP_out']['refdG'], \
            inputs['PSCP_out']['refddG'], refdU, absolutedU = \
                old_systems_dictionary(inputs['gen_in']['potential'], inputs['gen_in']['molecule'])
#        inputs['PSCP_out']['absolutedU'] = absolutedU

        # Writing out the input file with updated dG and ddG values
        with open(args.input_file, 'w') as yaml_file:
            yaml.dump(inputs, yaml_file, default_flow_style=False)
#    else:
#        absolutedU = np.array(inputs['PSCP_out']['absolutedU']).astype(float)
#        refdU = absolutedU - absolutedU[0]

    if inputs['temp_in']['run_temperature'] == True:
        dGvsT(Temperatures=np.array(inputs['temp_in']['temperatures'].split()).astype(float),
              Temperatures_unsampled=np.array(inputs['temp_in']['temperatures_unsampled'].split()).astype(float),
              Pressure=inputs['gen_in']['pressure'],
              Molecules=inputs['gen_in']['number_of_molecules'], molecule=inputs['gen_in']['molecule'],
              Independent=independent, potential=inputs['gen_in']['potential'],
              simulation=inputs['temp_in']['simulation_package'], hinge=inputs['gen_in']['hinge'],
              Polymorphs=inputs['gen_in']['polymorph_num'].split(), refT=inputs['PSCP_in']['PSCP_temperature'],
              refdG=inputs['PSCP_out']['refdG'], refddG=inputs['PSCP_out']['refddG'])
        #,refdU=refdU,
#              absolutedU=absolutedU)




