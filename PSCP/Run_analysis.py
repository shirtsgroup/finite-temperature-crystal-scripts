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
from dA_endpoint_MBAR import dA_endpoint_MBAR
from dA_MBAR import dA_MBAR
from dGvsT import dGvsT

if __name__ == '__main__':
    # Importing user specified input .yaml file
    import argparse
    parser = argparse.ArgumentParser(description='Runs analysis of a PSCP and ')
    parser.add_argument('-i', '--input_file', dest='input_file', default='input.yaml',
                        help='Input file containing all parameters to set up the directories to run the MD simulations')

    args = parser.parse_args()
    inputs = RS.yaml_loader(args.input_file)

    # Determining the number of independent molecules
    if inputs['gen_in']['independent'] == 'same':
        independent = inputs['gen_in']['number_of_molecules']
    else:
        independent = inputs['gen_in']['independent']

    # Run through all PSCP states to compute dA
    if inputs['PSCP_out']['dA'] == None:
        dA = np.zeros((len(inputs['gen_in']['polymorph_num'].split()), len(inputs['PSCP_in']['run_PSCP'])))
        dA[:, :] = np.nan
        ddA = np.zeros((len(inputs['gen_in']['polymorph_num'].split()), len(inputs['PSCP_in']['run_PSCP'])))
        save_dA = True
    else:
        dA = np.array(inputs['PSCP_out']['dA']).astype(float)
        ddA = np.array(inputs['PSCP_out']['ddA']).astype(float)
        save_dA = False

    interactions_count = 0
    restraints_count = 0

    for i, run in enumerate(inputs['PSCP_in']['run_PSCP']):
        # Getting the correct directory name
        if inputs['PSCP_in']['k'][i] == inputs['PSCP_in']['k'][i + 1]:
            if interactions_count == 0:
                extra_name = ''
            else:
                extra_name = '_' + str(interactions_count)
            directory_name = 'interactions' + extra_name
            interactions_count += 1
        else:
            if restraints_count == 0:
                extra_name = ''
            else:
                extra_name = '_' + str(restraints_count)
            directory_name = 'restraints' + extra_name
            restraints_count += 1

        # Running the analysis for this PSCP step
        if run and np.any(np.isnan(dA[:, i])):
            dA[:, i], ddA[:, i] = dA_MBAR(spacing=inputs['PSCP_in']['spacing'][i], exponent=inputs['PSCP_in']['exponent'][i],
                                          polymorphs=inputs['gen_in']['polymorph_num'],
                                          Molecules=inputs['gen_in']['number_of_molecules'], Independent=independent,
                                          Temp=inputs['PSCP_in']['PSCP_temperature'],
                                          bonds=inputs['PSCP_in']['run_bonded_interactions'],
                                          primary_directory=directory_name)

            inputs['PSCP_out']['dA'] = dA.tolist()
            inputs['PSCP_out']['ddA'] = ddA.tolist()
        
            with open(args.input_file, 'w') as yaml_file:
               yaml.dump(inputs, yaml_file, default_flow_style=False)

#    # Computing the free energy from PSCP
#    if inputs['PSCP_in']['run_restraints']:
#        # Computing the free energy to restrain the molecules
#        # Values are based on the restrained non-interacting state being the reference
#        dA_L, ddA_L = dA_Lambda_MBAR(MinL=inputs['PSCP_in']['min_lambda'], MaxL=inputs['PSCP_in']['max_lambda'],
#                                     dL=inputs['PSCP_in']['lambda_spacing'], GAMMA=inputs['PSCP_in']['gamma'],
#                                     exponent=inputs['PSCP_in']['lambda_exponent'],
#                                     polymorphs=inputs['gen_in']['polymorph_num'],
#                                     Molecules=inputs['gen_in']['number_of_molecules'], Independent=independent,
#                                     Temp=inputs['PSCP_in']['PSCP_temperature'], Pressure=inputs['gen_in']['pressure'],
#                                     potential=inputs['gen_in']['potential'],
#                                     hinge=inputs['gen_in']['hinge'])
#        inputs['PSCP_out']['dLambda'] = dA_L.tolist()
#        inputs['PSCP_out']['ddLambda'] = ddA_L.tolist()
#
#    if inputs['PSCP_in']['run_interactions']:
#        # Computing the free energy to turn off the interactions
#        dA_G, ddA_G = dA_Gamma_MBAR(MINGAMMA=inputs['PSCP_in']['min_gamma'], MAXGAMMA=inputs['PSCP_in']['max_gamma'],
#                                    GSPACING=inputs['PSCP_in']['gamma_spacing'], LAMBDA=inputs['PSCP_in']['lambda'],
#                                    exponent=inputs['PSCP_in']['gamma_exponent'],
#                                    polymorphs=inputs['gen_in']['polymorph_num'],
#                                    Molecules=inputs['gen_in']['number_of_molecules'], Independent=independent,
#                                    Temp=inputs['PSCP_in']['PSCP_temperature'],
#                                    Pressure=inputs['gen_in']['pressure'], k=inputs['PSCP_in']['k_max'],
#                                    potential=inputs['gen_in']['potential'], hinge=inputs['gen_in']['hinge'],
#                                    bonds=inputs['PSCP_in']['run_bonded_interactions'])
#        inputs['PSCP_out']['dGamma'] = dA_G.tolist()
#        inputs['PSCP_out']['ddGamma'] = ddA_G.tolist()
#
#        # Computing the free energy for the endpoints
#        dA_end, ddA_end = dA_endpoint_MBAR(polymorphs=inputs['gen_in']['polymorph_num'],
#                                           Molecules=inputs['gen_in']['number_of_molecules'], Independent=independent,
#                                           Temp=inputs['PSCP_in']['PSCP_temperature'])
#        inputs['PSCP_out']['dEnd'] = dA_end.tolist()
#        inputs['PSCP_out']['ddEnd'] = ddA_end.tolist()

    if (inputs['PSCP_out']['dG'] == None) or (np.any(np.isnan(inputs['PSCP_out']['dG']))):

        # Adding the free energy differences to the inputs to be saved
        dG = np.zeros(len(inputs['gen_in']['polymorph_num'].split()))
        ddG = np.zeros(len(inputs['gen_in']['polymorph_num'].split()))
        for i, ply in enumerate(inputs['gen_in']['polymorph_num'].split()):
            dG[i] = np.sum(dA[i] - dA[0])
            ddG[i] = np.sqrt(np.sum(ddA[i]**2))
        inputs['PSCP_out']['dG'] = dG.tolist()
        inputs['PSCP_out']['ddG'] = ddG.tolist()

        # Writing out the input file with updated dG and ddG values
        with open(args.input_file, 'w') as yaml_file:
            yaml.dump(inputs, yaml_file, default_flow_style=False)
    else:
        print("Using user specified values of dG!")
        print("   If this is wrong, please remove dG and ddG from the input file")

    # Determing the free energy across the entire temperature range
    if inputs['temp_in']['run_temperature'] == True:
        if inputs['temp_in']['temperatures_unsampled'] == None:
            inputs['temp_in']['temperatures_unsampled'] = []
        else:
            inputs['temp_in']['temperatures_unsampled'] = np.array(inputs['temp_in']['temperatures_unsampled'].split()).astype(float)
        dGvsT(Temperatures=np.array(inputs['temp_in']['temperatures'].split()).astype(float),
              Temperatures_unsampled=inputs['temp_in']['temperatures_unsampled'],
              Pressure=inputs['gen_in']['pressure'],
              Molecules=inputs['gen_in']['number_of_molecules'], molecule=inputs['gen_in']['molecule'],
              Independent=independent, potential=inputs['gen_in']['potential'],
              simulation=inputs['temp_in']['simulation_package'], hinge=inputs['gen_in']['hinge'],
              Polymorphs=inputs['gen_in']['polymorph_num'].split(), refT=inputs['PSCP_in']['PSCP_temperature'],
              refdG=inputs['PSCP_out']['dG'], refddG=inputs['PSCP_out']['ddG'], output_directory=inputs['gen_in']['output_directory'])

