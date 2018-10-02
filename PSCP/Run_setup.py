#!/usr/bin/env python

from __future__ import print_function
import os
import sys
import subprocess
import yaml

def yaml_loader(file_path):
    # Loads in a ymal file
    with open(file_path, "r") as input_file:
        data = yaml.load(input_file)
    return data

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Sets up directories to run MD simulations for the free energy differences between organic crystals')
    parser.add_argument('-i', '--input_file', dest='input_file', default='input.yaml', help='Input file containing all parameters to set up the directories to run the MD simulations')

    args = parser.parse_args()
    inputs = yaml_loader(args.input_file)

    path = os.path.realpath(__file__).strip('Run_setup.py')
    subprocess.call([path + 'setup-scripts/setup_Restraints -n "' + (inputs['polymorph_num']) + '"'
                                            + ' -M ' + str(inputs['molecule']) 
                                            + ' -U ' + str(inputs['min_lambda']) 
                                            + ' -R ' + str(inputs['max_lambda']) 
                                            + ' -L ' + str(inputs['lambda_spacing']) 
                                            + ' -f ' + str(inputs['lambda_exponent']) 
                                            + ' -G ' + str(inputs['gamma']) 
                                            + ' -T ' + str(inputs['PSCP_temperature']) 
                                            + ' -P ' + str(inputs['pressure']) 
                                            + ' -N ' + str(inputs['number_of_molecules']) 
                                            + ' -I ' + str(inputs['independent'])
                                            + ' -e ' + str(inputs['lambda_equil_steps']) 
                                            + ' -p ' + str(inputs['lambda_prod_steps']) 
                                            + ' -i ' + str(inputs['integrator']) 
                                            + ' -t ' + str(inputs['thermostat']) 
                                            + ' -a ' + str(inputs['cores']) 
                                            + ' -k ' + str(inputs['k_min']) 
                                            + ' -K ' + str(inputs['k_max']) 
                                            + ' -r ' + str(inputs['cutoff']) 
                                            + ' -u ' + str(inputs['potential'])
                                            + ' -h ' + str(inputs['hinge'])], shell=True)

    subprocess.call([path + 'setup-scripts/setup_Interactions -n "' + str(inputs['polymorph_num']) + '"'
                                              + ' -M ' + str(inputs['molecule']) 
                                              + ' -A ' + str(inputs['max_gamma']) 
                                              + ' -B ' + str(inputs['min_gamma']) 
                                              + ' -g ' + str(inputs['gamma_spacing']) 
                                              + ' -f ' + str(inputs['gamma_exponent']) 
                                              + ' -L ' + str(inputs['lambda']) 
                                              + ' -T ' + str(inputs['PSCP_temperature']) 
                                              + ' -P ' + str(inputs['pressure']) 
                                              + ' -N ' + str(inputs['number_of_molecules']) 
                                              + ' -I ' + str(inputs['independent'])
                                              + ' -e ' + str(inputs['gamma_equil_steps']) 
                                              + ' -p ' + str(inputs['gamma_prod_steps']) 
                                              + ' -i ' + str(inputs['integrator']) 
                                              + ' -t ' + str(inputs['thermostat']) 
                                              + ' -a ' + str(inputs['cores']) 
                                              + ' -K ' + str(inputs['k_max']) 
                                              + ' -r ' + str(inputs['cutoff']) 
                                              + ' -u ' + str(inputs['potential'])
                                              + ' -h ' + str(inputs['hinge'])], shell=True)
    
    subprocess.call([path + 'setup-scripts/setup_Temperature -T "' + inputs['temperatures'] + '" -C ' + str(inputs['charge']) 
                                           + ' -n "' + str(inputs['polymorph_num']) + '"'
                                           + ' -M ' + str(inputs['molecule']) 
                                           + ' -P ' + str(inputs['pressure']) 
                                           + ' -N ' + str(inputs['number_of_molecules']) 
                                           + ' -I ' + str(inputs['independent'])
                                           + ' -e ' + str(inputs['temp_equil_steps']) 
                                           + ' -p ' + str(inputs['temp_prod_steps']) 
                                           + ' -i ' + str(inputs['integrator']) 
                                           + ' -t ' + str(inputs['thermostat']) 
                                           + ' -b ' + str(inputs['barostat'])
                                           + ' -a ' + str(inputs['cores']) 
                                           + ' -k ' + str(inputs['k_max']) 
                                           + ' -r ' + str(inputs['cutoff']) 
                                           + ' -u ' + str(inputs['potential']) 
                                           + ' -z ' + str(inputs['simulation_package'])
                                           + ' -h ' + str(inputs['hinge']) 
                                           + ' -o ' + str(inputs['prodoutputs'])], shell=True)

    




