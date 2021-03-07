#!/usr/bin/env python

from __future__ import print_function
import os
import sys
import subprocess
import yaml
#import mdtraj as md
import numpy as np
from scipy.optimize import fsolve
from scipy.special import erf
from pymbar.timeseries import detectEquilibration
import panedr as pdr
from scipy.special import erf

path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, path + '/setup-scripts')
import setup_directories as setup
from setup_boltz_QHA import setup_boltz_QHA

def setdefault(input_data, default_values):
    # Function to fill in the input_data if the default values are not set
    for k in default_values:
        if isinstance(default_values[k], dict):
            setdefault(input_data.setdefault(k, {}), default_values[k])
        else:
            input_data.setdefault(k, default_values[k])

def yaml_loader(file_path):
    # Loads in a ymal file
    with open(file_path, "r") as input_file:
        data = yaml.load(input_file, Loader=yaml.FullLoader)
    # Load in the default values
    with open(path + '/setup-scripts/default.yaml', "r") as default_file:
        default_input = yaml.load(default_file, Loader=yaml.FullLoader)

    # Setting the default v   alues if not specified
    setdefault(data, default_input)
    return data

def load_potenergy(fil):
    U = pdr.edr_to_df(fil)
    U = np.array(U['Potential'])
    [t0,_,_] = detectEquilibration(U) 
    return U[t0:]

def prob_diff(dT, P, T0, dmu, b_mu, dsig, b_sig):
    # Determines the difference between the desired probability overlap P and then probability overlap using dT
    dT = np.absolute(dT)
    mu_1 = dmu * T0 + b_mu
    mu_2 = dmu * (T0 + dT) + b_mu
    sig_1 = dsig * T0 + b_sig
    sig_2 = dsig * (T0 + dT) + b_sig
    return (P - compute_prob(mu_1, mu_2, sig_1, sig_2))**2
#    return (P - determine_exchange_overlap([T0+dT, T0], [mu_2, mu_1], [sig_2, sig_1]))**2

def compute_prob(mu_1, mu_2, sig_1, sig_2):
    c = (mu_2 * sig_1**2 - sig_2*(mu_1*sig_2 + sig_1*np.sqrt((mu_1 - mu_2)**2 + 2*(sig_1**2 - sig_2**2)*np.log10(sig_1/sig_2)))) / (sig_1**2 - sig_2**2)
    P = 1 - 0.5*erf((c- mu_1) / (np.sqrt(2) *sig_1)) + 0.5*erf((c- mu_2) / (np.sqrt(2) *sig_2))
    return P

def determine_exchange_overlap(T, mu, sigma):
    k_b = 0.0019872041 # Units in kcal/(mol*K)
    mu_12 = mu[0] - mu[1]
    sigma_12 = np.sqrt(sigma[0]**2 + sigma[1]**2)
    C = (1 / (k_b * T[0]) - 1 / (k_b * T[1]))
    P = 0.5 * (1 + erf(-mu_12 / (sigma_12 * np.sqrt(2)))) + 0.5 * np.exp(C * mu_12 + 0.5 * C**2 * sigma_12**2) * (1 + erf((mu_12 + C * sigma_12**2) / (sigma_12 * np.sqrt(2))))
    return P

def return_dT(dmu, b_mu, dsig, b_sig, P, T0):
    # Finds the dT that provides the desired proability overlap
    x = fsolve(prob_diff, 0.001, args=(P, T0, dmu, b_mu, dsig, b_sig))
    #mus = [dmu * (T0 + x) + b_mu, dmu * T0 + b_mu]
    #sigmas = [dsig * (T0 + x) + b_sig, dsig * T0 + b_sig]
    #Ts = [T0 + x, T0]
    #print(determine_exchange_overlap(Ts, mus, sigmas), compute_prob(mus[1], mus[0], sigmas[1], sigmas[0]))
    return x

def setup_ReplicaExchange_temperatures(inputs, check):
    # Extracing specific data from the inputs
    polymorph_num = inputs['gen_in']['polymorph_num'].split()
    temperatures = np.array(inputs['temp_in']['temperatures'].split()).astype(float)

    # Setting up arrays to store values from previous run
    average = np.zeros((len(polymorph_num), len(temperatures)))
    average[:, :] = np.nan
    st_dev = np.zeros((len(polymorph_num), len(temperatures)))
    st_dev[:, :] = np.nan

    # Pulling data from previous run and reporting if certain temperatures are missing
    run_replica_exchange = True
    for i in range(len(polymorph_num)):
        for j in range(len(temperatures)):
            if os.path.isfile(polymorph_num[i] + '/temperature/' + str(int(j)) + '/PROD.edr'):
                U = load_potenergy(polymorph_num[i] + '/temperature/' + str(int(j)) + '/PROD.edr')
                average[i, j] = np.mean(U)
                st_dev[i, j] = np.std(U)
            else:
                #run_replica_exchange = False
                print('WARNING: File /' + polymorph_num[i] + '/temperature/' + str(int(j)) + '/potenergy.xvg' + ' is missing')

    # Ending the run if previous temperature runs are not present
    if run_replica_exchange == False:
        print('Necessary files for Replica Exhange are missing, system exiting!')
        sys.exit()

    # Using the average and standard deviation of the previous runs to determine the temperature spacing fro replica exchange
# NSA: right now this is just taking the first polymorph to determine the temperature spacing
    temperatures = temperatures.astype(float)
    #dmu = np.zeros(len(polymorph_num[0]))
    #b_mu = np.zeros(len(polymorph_num[0]))
    #dsig = np.zeros(len(polymorph_num[0]))
    #b_sig = np.zeros(len(polymorph_num[0]))

    #print(len(polymorph_num))
    #for i in range(len(polymorph_num[0])):
    dmu, b_mu = np.polyfit(temperatures, average[0, :], 1)
    dsig, b_sig = np.polyfit(temperatures, st_dev[0, :], 1)
    
    temp = float(inputs['rep_exch_in']['T_min'])
    T_out = [temp]
    while temp < float(inputs['rep_exch_in']['T_max']):
        dt = np.around(return_dT(dmu, b_mu, dsig, b_sig, inputs['rep_exch_in']['prob_overlap'], temp), 4)
        temp += np.absolute(dt[0])
        if temp < float(inputs['rep_exch_in']['T_max']):
            T_out.append(temp)
        else:
            T_out.append(float(inputs['rep_exch_in']['T_max']))

    # Correcting for the number of nodes and changing the temperautre back to a string
    T_out, nodes = correct_for_nodes(T_out, inputs['rep_exch_in']['cores_per_replica'], inputs['rep_exch_in']['cores_per_node'])

    if check == True:
        print('Temperatures: ', str(T_out))
        print('Nodes: ', nodes)
        sys.exit()

    # Turning the desired temperatures into a string for the input file
    T = ""
    for i in T_out:
        T += str(np.around(i, 4)) + " "
    return T, nodes
    

def correct_for_nodes(T, cores_per_replica, cores_per_node):
    # Determining the number of nodes required to run REMD
    nodes = len(T) * cores_per_replica / cores_per_node

    # Determineing the number of extra temperatures that need to be added to have an integer value of nodes
    extra_T = int(np.ceil(nodes) * cores_per_node / cores_per_replica - len(T))

    # Adding Extra Temperatures
    if extra_T != 0:
        for i in range(extra_T):
            append = False
            while append == False:
                T_rand = np.around(np.random.random(1)[0] * np.max(T), 1)
                if (T_rand > np.min(T)) and (np.all(T_rand != T)):
                    append = True
                    T = np.append(T, T_rand)
    return np.sort(T), int(np.ceil(nodes))


if __name__ == '__main__':
    import argparse
    # Importing the input file
    parser = argparse.ArgumentParser(description='Sets up directories to run MD simulations for the free energy '
                                                 'differences between organic crystals')
    parser.add_argument('-i', '--input_file', dest='input_file', default='input.yaml',
                        help='Input file containing all parameters to set up the directories to run the MD simulations')
    parser.add_argument('--REP', action='store_true',
                        help='Setups the system for replica exchange if previous temperature have been run as reference')
    parser.add_argument('--check', action='store_true')

    args = parser.parse_args()

    # Loading the input file
    inputs = yaml_loader(args.input_file)

    if args.REP:
        inputs['temp_in']['temperatures'], RE_nodes = setup_ReplicaExchange_temperatures(inputs, args.check)
        subprocess.call(['mv', args.input_file, args.input_file.strip('.yaml') + '_prep.yaml'])
        with open(args.input_file, 'w') as yaml_file:
            yaml.dump(inputs, yaml_file, default_flow_style=False)

    # Setup a directory for re-setting up the simulations. 
    # Instead of looking for a general template path directory, a re_setup directory will be created with all necessary input files
    if (inputs['gen_in']['template_path'] != 're_setup'):
        subprocess.call(['mkdir', 're_setup'])
        for i in inputs['gen_in']['polymorph_num'].split():
            setup.re_setup(polymorph_num=i, #temperature=inputs['temp_in']['temperatures'].split()[0], 
                           #pressure=inputs['gen_in']['pressure'], 
                           molecule=inputs['gen_in']['molecule'],
                           number_of_molecules=inputs['gen_in']['number_of_molecules'], 
                           independent=inputs['gen_in']['independent'], 
                           potential=inputs['gen_in']['potential'],
                           simulation=inputs['temp_in']['simulation_package'], 
                           templatepath=inputs['gen_in']['template_path'])
        inputs['gen_in']['template_path'] = 're_setup'
        with open(args.input_file, 'w') as yaml_file:
            yaml.dump(inputs, yaml_file, default_flow_style=False)


    # Creating a directory for each polymorph
    for i in inputs['gen_in']['polymorph_num'].split():
        if args.REP:
            subprocess.call(['mv', i, i + '_prep'])

        if not os.path.isdir(i):
            subprocess.call(['mkdir', i])

    # Running setup for the temperature scan
    if inputs['temp_in']['run_temperature'] == True:
        run_production = False

        # Determining if there needs to be a production for individual temperatures
        # If replica exchange is occuring, the production run for the submit_cluster.slurm scripts is turned off
        if (inputs['temp_in']['temp_prod_steps'] > 0) and not args.REP:
            run_production = True

        # Setting up the directories for the temperature run
        setup.setup_temperature(inputs, run_production)

        if args.REP:
            # Creating a string with the list of directories 0,....N
            DIRS = ''
            for i in range(len(inputs['temp_in']['temperatures'].split())):
                DIRS += str(i) + ','
            DIRS = DIRS.strip(',')

            # Determining the number of total replicas
            process_num = len(inputs['temp_in']['temperatures'].split())

            # Determining the number of total exchanges to consider
            exchange_num = process_num ** 3

            # Setting up the replica exchange script to be run from each directory
            for i in inputs['gen_in']['polymorph_num'].split():
                setup.setup_replica_exchange(int(RE_nodes), DIRS, int(process_num), int(exchange_num),
                                             str(i) + '/temperature')
    # Running setups for the pseudosupercritical path
    setup.setup_PSCP(inputs)

    # Running setup for boltzmann weighting of QHA
    ## This is currently only setup for Gromacs MD and Tinker QHA compatability
    setup_boltz_QHA(inputs)




