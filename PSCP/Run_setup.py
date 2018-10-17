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
    parser = argparse.ArgumentParser(description='Sets up directories to run MD simulations for the free energy '
                                                 'differences between organic crystals')
    parser.add_argument('-i', '--input_file', dest='input_file', default='input.yaml',
                        help='Input file containing all parameters to set up the directories to run the MD simulations')

    args = parser.parse_args()
    inputs = yaml_loader(args.input_file)

    path = os.path.realpath(__file__).strip('Run_setup.py')
    if inputs['PSCP_in']['run_restraints'] == True:
        subprocess.call([path + 'setup-scripts/setup_Restraints -n "' + (inputs['gen_in']['polymorph_num']) + '"'
                                            + ' -M ' + str(inputs['gen_in']['molecule'])
                                            + ' -U ' + str(inputs['PSCP_in']['min_lambda'])
                                            + ' -R ' + str(inputs['PSCP_in']['max_lambda'])
                                            + ' -L ' + str(inputs['PSCP_in']['lambda_spacing'])
                                            + ' -f ' + str(inputs['PSCP_in']['lambda_exponent'])
                                            + ' -G ' + str(inputs['PSCP_in']['gamma'])
                                            + ' -T ' + str(inputs['PSCP_in']['PSCP_temperature'])
                                            + ' -P ' + str(inputs['gen_in']['pressure'])
                                            + ' -N ' + str(inputs['gen_in']['number_of_molecules'])
                                            + ' -I ' + str(inputs['gen_in']['independent'])
                                            + ' -e ' + str(inputs['PSCP_in']['lambda_equil_steps'])
                                            + ' -p ' + str(inputs['PSCP_in']['lambda_prod_steps'])
                                            + ' -i ' + str(inputs['gen_in']['integrator'])
                                            + ' -t ' + str(inputs['gen_in']['thermostat'])
                                            + ' -a ' + str(inputs['gen_in']['cores'])
                                            + ' -k ' + str(inputs['PSCP_in']['k_min'])
                                            + ' -K ' + str(inputs['PSCP_in']['k_max'])
                                            + ' -r ' + str(inputs['gen_in']['cutoff'])
                                            + ' -u ' + str(inputs['gen_in']['potential'])
                                            + ' -h ' + str(inputs['gen_in']['hinge'])], shell=True)

    if inputs['PSCP_in']['run_interactions'] == True:
        subprocess.call([path + 'setup-scripts/setup_Interactions -n "' + str(inputs['gen_in']['polymorph_num']) + '"'
                                              + ' -M ' + str(inputs['gen_in']['molecule'])
                                              + ' -A ' + str(inputs['PSCP_in']['max_gamma'])
                                              + ' -B ' + str(inputs['PSCP_in']['min_gamma'])
                                              + ' -g ' + str(inputs['PSCP_in']['gamma_spacing'])
                                              + ' -f ' + str(inputs['PSCP_in']['gamma_exponent'])
                                              + ' -L ' + str(inputs['PSCP_in']['lambda'])
                                              + ' -T ' + str(inputs['PSCP_in']['PSCP_temperature'])
                                              + ' -P ' + str(inputs['gen_in']['pressure'])
                                              + ' -N ' + str(inputs['gen_in']['number_of_molecules'])
                                              + ' -I ' + str(inputs['gen_in']['independent'])
                                              + ' -e ' + str(inputs['PSCP_in']['gamma_equil_steps'])
                                              + ' -p ' + str(inputs['PSCP_in']['gamma_prod_steps'])
                                              + ' -i ' + str(inputs['gen_in']['integrator'])
                                              + ' -t ' + str(inputs['gen_in']['thermostat'])
                                              + ' -a ' + str(inputs['gen_in']['cores'])
                                              + ' -K ' + str(inputs['PSCP_in']['k_max'])
                                              + ' -r ' + str(inputs['gen_in']['cutoff'])
                                              + ' -u ' + str(inputs['gen_in']['potential'])
                                              + ' -h ' + str(inputs['gen_in']['hinge'])], shell=True)

    if inputs['temp_in']['run_temperature'] == True:
        subprocess.call([path + 'setup-scripts/setup_Temperature -T "' + inputs['temp_in']['temperatures'] + '" -C ' + str(inputs['temp_in']['charge'])
                                           + ' -n "' + str(inputs['gen_in']['polymorph_num']) + '"'
                                           + ' -M ' + str(inputs['gen_in']['molecule'])
                                           + ' -P ' + str(inputs['gen_in']['pressure'])
                                           + ' -N ' + str(inputs['gen_in']['number_of_molecules'])
                                           + ' -I ' + str(inputs['gen_in']['independent'])
                                           + ' -e ' + str(inputs['temp_in']['temp_equil_steps'])
                                           + ' -p ' + str(inputs['temp_in']['temp_prod_steps'])
                                           + ' -i ' + str(inputs['gen_in']['integrator'])
                                           + ' -t ' + str(inputs['gen_in']['thermostat'])
                                           + ' -b ' + str(inputs['temp_in']['barostat'])
                                           + ' -a ' + str(inputs['gen_in']['cores'])
                                           + ' -k ' + str(inputs['temp_in']['k_max'])
                                           + ' -r ' + str(inputs['gen_in']['cutoff'])
                                           + ' -u ' + str(inputs['gen_in']['potential'])
                                           + ' -z ' + str(inputs['temp_in']['simulation_package'])
                                           + ' -h ' + str(inputs['gen_in']['hinge'])
                                           + ' -o ' + str(inputs['temp_in']['prodoutputs'])], shell=True)

    




