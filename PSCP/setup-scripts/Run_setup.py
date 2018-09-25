#!/usr/bin/env python

from __future__ import print_function
import os
import sys
import subprocess
import yaml

def ymal_loader(file_path):
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

    subprocess.call(['setup_Restraints', '-n ' + inputs['polymorph_num'], '-M ' + inputs['molecule'], '-U ' + inputs['min_lambda'], 
        '-R ' + inputs['max_lambda'], '-L ' + inputs['lambda_spacing'], '-f ' + inputs['lambda_exponent'], '-G ' + inputs['gamma'], 
        '-T ' + inputs['PSCP_temperature'], '-P ' + inputs['pressure'], '-N ' + inputs['number_of_molecules'], '-I ' + inputs['independent']
        '-e ' + inputs['lambda_equil_steps'], '-p ' + inputs['lambda_prod_steps'], '-i ' + inputs['integrator'], '-t ' + inputs['thermostat'], 
        '-a ' + inputs['cores'], '-k ' + inputs['k_min'], '-K ' + inputs['k_max'], '-r ' + inputs['cutoff'], '-u ' + inputs['potential'],
        '-h ' + inputs['hinge']])

    subprocess.call(['setup_Interactions', '-n ' + inputs['polymorph_num'], '-M ' + inputs['molecule'], '-A ' + inputs['max_gamma'], 
        '-B ' + inputs['min_gamma'], '-g ' + inputs['gamma_spacing'], '-f ' + inputs['gamma_exponent'], '-L ' + inputs['lambda'], 
        '-T ' + inputs['PSCP_temperature'], '-P ' + inputs['pressure'], '-N ' + inputs['number_of_molecules'], '-I ' + inputs['independent']
        '-e ' + inputs['gamma_equil_steps'], '-p ' + inputs['gamma_prod_steps'], '-i ' + inputs['integrator'], '-t ' + inputs['thermostat'], 
        '-a ' + inputs['cores'], '-K ' + inputs['k_max'], '-r ' + inputs['cutoff'], '-u ' + inputs['potential'],
        '-h ' + inputs['hinge']])

    subprocess.call(['setup_Temperature', '-C ' inputs['charge'], '-n ' + inputs['polymorph_num'], '-M ' + inputs['molecule'], 
        '-T ' + inputs['PSCP_temperature'], '-P ' + inputs['pressure'], '-N ' + inputs['number_of_molecules'], '-I ' + inputs['independent']
        '-e ' + inputs['temp_equil_steps'], '-p ' + inputs['temp_prod_steps'], '-i ' + inputs['integrator'], '-t ' + inputs['thermostat'], '-b ' + inputs['barostat'],
        '-a ' + inputs['cores'], '-k ' + inputs['k_max'], '-r ' + inputs['cutoff'], '-u ' + inputs['potential'], '-z ' + inputs['simulation'],
        '-h ' + inputs['hinge'], '-o ' + inputs['prodoutputs']])

    




