#!/usr/bin/env python
import sys
import subprocess
import os
import numpy as np
from resize_gro_individual import run_resize_gro_individual
from interpolate_itp import interpolate_itp

path = os.path.realpath(__file__).strip('setup_directories.py')

def setup_temperature(inputs, run_production):
    # Python script to automatically set up simulations over a given range of temperatures
    # Original script written in bash by: Eric Dybeck on 03/31/2014
    # Converted to python by: Nate Abraham on 01/25/2019

    # Going through the setup for each polymorph
    for i in inputs['gen_in']['polymorph_num'].split():
        subprocess.call(['mkdir', i + '/temperature'])
        count = 0
        # Setting up each temperature for the polymorph
        for j in inputs['temp_in']['temperatures'].split():
            setup_molecule(polymorph_num=i, temperature=j, pressure=inputs['gen_in']['pressure'],
                           molecule=inputs['gen_in']['molecule'],
                           number_of_molecules=inputs['gen_in']['number_of_molecules'],
                           independent=inputs['gen_in']['independent'],
                           equil_steps=inputs['temp_in']['temp_equil_steps'],
                           prod_steps=inputs['temp_in']['temp_prod_steps'],
                           prodoutputs=inputs['temp_in']['prodoutputs'], integrator=inputs['gen_in']['integrator'],
                           thermostat=inputs['gen_in']['thermostat'], barostat=inputs['temp_in']['barostat'],
                           cores=inputs['gen_in']['cores'], k_max=inputs['temp_in']['k_max'],
                           cutoff=inputs['gen_in']['cutoff'], potential=inputs['gen_in']['potential'],
                           simulation=inputs['temp_in']['simulation_package'], jobpath=i + '/temperature/' + str(count),
                           templatepath=inputs['gen_in']['template_path'], anneal_temp=inputs['gen_in']['anneal_temp'],
                           anneal_steps=inputs['temp_in']['temp_anneal_steps'], run_production=run_production,
                           charge=inputs['temp_in']['charge'], hinge=inputs['gen_in']['hinge'],
                           submission_script=inputs['gen_in']['submission_script'])
            count += 1

def setup_restraints(inputs):
    # Going through the setup for each polymorph
    for i in inputs['gen_in']['polymorph_num'].split():
        subprocess.call(['mkdir', i + '/restraints'])  # making a directory for turning off the interactions
        lambd = inputs['PSCP_in']['min_lambda']
        # Setting up the directory for each lambda value between the minimum and maximum lambda
        while lambd <= inputs['PSCP_in']['max_lambda']:
            setup_molecule(polymorph_num=i, temperature=inputs['PSCP_in']['PSCP_temperature'],
                           pressure=inputs['gen_in']['pressure'], molecule=inputs['gen_in']['molecule'],
                           number_of_molecules=inputs['gen_in']['number_of_molecules'],
                           independent=inputs['gen_in']['independent'],
                           equil_steps=inputs['PSCP_in']['lambda_equil_steps'],
                           prod_steps=inputs['PSCP_in']['lambda_prod_steps'],
                           prodoutputs=inputs['temp_in']['prodoutputs'], integrator=inputs['gen_in']['integrator'],
                           thermostat=inputs['gen_in']['thermostat'], barostat=inputs['temp_in']['barostat'],
                           cores=inputs['gen_in']['cores'], k_min=inputs['PSCP_in']['k_min'],
                           k_max=inputs['PSCP_in']['k_max'], lambd=lambd,
                           min_lambda=inputs['PSCP_in']['min_lambda'], max_lambda=inputs['PSCP_in']['max_lambda'],
                           lambda_spacing=inputs['PSCP_in']['lambda_spacing'],
                           lambda_exponent=inputs['PSCP_in']['lambda_exponent'], gamma=inputs['PSCP_in']['gamma'],
                           cutoff=inputs['gen_in']['cutoff'], potential=inputs['gen_in']['potential'],
                           simulation=inputs['temp_in']['simulation_package'], ensemble='NVT',
                           jobpath=i + '/restraints/' + str(lambd), templatepath=inputs['gen_in']['template_path'],
                           anneal_temp=inputs['gen_in']['anneal_temp'], anneal_steps=0, run_production=True,
                           charge=inputs['temp_in']['charge'], hinge=inputs['gen_in']['hinge'],
                           submission_script=inputs['gen_in']['submission_script'])
            lambd += inputs['PSCP_in']['lambda_spacing']


def setup_interactions(inputs):
    # Going through the setup for each polymorph
    for i in inputs['gen_in']['polymorph_num'].split():
         subprocess.call(['mkdir', i + '/interactions'])  # making a directory for turning off the interactions
         gamma = inputs['PSCP_in']['min_gamma']
         # Setting up the directory for each gamma value between the minimum and maximum gamma
         while gamma <= inputs['PSCP_in']['max_gamma']:
             setup_molecule(polymorph_num=i, temperature=inputs['PSCP_in']['PSCP_temperature'],
                            pressure=inputs['gen_in']['pressure'], molecule=inputs['gen_in']['molecule'],
                            number_of_molecules=inputs['gen_in']['number_of_molecules'],
                            independent=inputs['gen_in']['independent'],
                            equil_steps=inputs['PSCP_in']['gamma_equil_steps'],
                            prod_steps=inputs['PSCP_in']['gamma_prod_steps'],
                            prodoutputs=inputs['temp_in']['prodoutputs'], integrator=inputs['gen_in']['integrator'],
                            thermostat=inputs['gen_in']['thermostat'], barostat=inputs['temp_in']['barostat'],
                            cores=inputs['gen_in']['cores'], k_min=inputs['PSCP_in']['k_min'],
                            k_max=inputs['PSCP_in']['k_max'], lambd=inputs['PSCP_in']['lambda'], gamma=gamma,
                            min_lambda=inputs['PSCP_in']['lambda'], max_lambda=inputs['PSCP_in']['lambda'],
                            min_gamma=inputs['PSCP_in']['min_gamma'], max_gamma=inputs['PSCP_in']['max_gamma'],
                            gamma_exponent=inputs['PSCP_in']['gamma_exponent'],
                            gamma_spacing=inputs['PSCP_in']['gamma_spacing'], cutoff=inputs['gen_in']['cutoff'],
                            potential=inputs['gen_in']['potential'], simulation=inputs['temp_in']['simulation_package'],
                            ensemble='NVT', jobpath=i + '/interactions/' + str(gamma),
                            templatepath=inputs['gen_in']['template_path'], anneal_temp=inputs['gen_in']['anneal_temp'],
                            anneal_steps=0, run_production=True, charge=inputs['temp_in']['charge'],
                            hinge=inputs['gen_in']['hinge'], submission_script=inputs['gen_in']['submission_script'], 
                            remove_bonded_interactions=inputs['PSCP_in']['run_bonded_interactions'])
             gamma += inputs['PSCP_in']['gamma_spacing']


def setup_replica_exchange(nodes, directories, process_number, exchange_number, jobpath):
    # Copying submission script into temperature direcotry
    subprocess.call(['cp', path + '/run_files/submit_cluster_REP.slurm', jobpath + '/'])

    # Replacing specific strings in the submit script to match the user input
    replace_string_in_text(jobpath + '/submit_cluster_REP.slurm', 'NODES', nodes)
#    replace_string_in_text(jobpath + '/submit_cluster_REP.slurm', 'DIRS', directories)
    replace_string_in_text(jobpath + '/submit_cluster_REP.slurm', 'NNPP', process_number)
#    replace_string_in_text(jobpath + '/submit_cluster_REP.slurm', 'NNEEXX', exchange_number)




def setup_mdp_lambdas(current_lambda, current_gamma, polymorph_num='all', min_lambda=0, max_lambda=100, 
                      lambda_spacing=-1, lambda_exponent=2, min_gamma=0, max_gamma=100,
                      gamma_spacing=-1, gamma_exponent=2, jobpath='DefaultPath', 
                      remove_bonded_interactions=False, equil_output_frequency=1000,
                      prod_output_frequency=1000):
    # Python script to automatically set up the mdp files to take on multiple lambda values
    # Original script written in bash by: Eric Dybeck on 09/12/2014
    # Converted to python by: Nate Abraham on 01/23/2019

    # Ensure that the parameters are properly entered
    # Lambda
    if (min_lambda < 0) or (max_lambda > 100) or (min_lambda > max_lambda):
        print('Minimum Lambda: ', min_lambda)
        print('Maximum Lambda: ', max_lambda)
        print('Is not a valid lambda range!')
        return

    if (lambda_spacing <= 0) or (lambda_spacing > 100):
        print('Invalid Lambda Spacing: ', lambda_spacing)
        return

    if (min_lambda == max_lambda) and (current_lambda != max_lambda):
        print('Minimum Lambda: ', min_lambda, ' Maximum Lambda: ', max_lambda, ' and Lambda: ', current_lambda,
              ' are not the same!')
        return

    # Gamma
    if (min_gamma < 0) or (max_gamma > 100) or (min_gamma > max_gamma):
        print('ERROR - ')
        print('Minimum Gamma: ', min_gamma)
        print('Maximum Gamma: ', max_gamma)
        print('Is not a valid lambda range!')
        return

    if gamma_spacing <= 0:
        print('ERROR - Invalid Gamma Spacing: ', gamma_spacing)
        return

    if (min_gamma == max_gamma) and (current_gamma != max_gamma):
        print('ERROR - Minimum Gamma: ', min_gamma, ' Maximum Gamma: ', max_gamma, ' and Gamma: ', current_gamma,
              ' are not the same!')
        return

    #JOBPATH
    if jobpath == 'DefaultPath':
        print('ERROR - Enter the job path!')
        return

    #If we have no harmonic restraints and full interactions, no need to proceed further
    if (min_lambda == max_lambda) and (max_lambda == 0) and (min_gamma == max_gamma) and (max_gamma == 100):
        print('ERROR - No lambda values added.')
        return

    #If we are adding harmonic restraints, we should not be changing gamma (and vice versa)
    if (min_lambda != max_lambda) and (min_gamma != max_gamma):
        print('ERROR - Harmonic restraints and Interactions changing simultaneously!')
        return

    # Change the free energy setting from 'no' to 'yes' and the output from 0 to nstxout
    replace_line_starting_with(jobpath + '/equilibration.mdp', 'free_energy', 'free_energy = yes')
    replace_line_starting_with(jobpath + '/production.mdp', 'free_energy', 'free_energy = yes')
    replace_line_starting_with(jobpath + '/equilibration.mdp', 'nstdhdl', 'nstdhdl = ' + str(equil_output_frequency))
    replace_line_starting_with(jobpath + '/production.mdp', 'nstdhdl', 'nstdhdl = ' + str(prod_output_frequency))

    # Setting arrays for 
    if min_lambda == max_lambda:
        # Setting up vectors for turning off interactions
        gammas = np.arange(min_gamma, max_gamma + 1, gamma_spacing)
        indicies = np.arange(0, (max_gamma - min_gamma) / gamma_spacing + 1, 1)
        lambda_points = np.ones(len(indicies))
        if gamma_exponent < 0:
            gamma_points = 1 - ((max_gamma - gammas) / max_gamma) ** abs(gamma_exponent)
        else:
            gamma_points = (gammas / max_gamma) ** abs(gamma_exponent)
        init_lambda = np.where(current_gamma == gammas)[0][0]

        # Setting interaction end points
#        replace_line_starting_with(jobpath + '/equilibration.mdp', 'couple-lambda0', ';couple-lambda0           = none') 
#        replace_line_starting_with(jobpath + '/production.mdp', 'couple-lambda0', ';couple-lambda0           = none')
#        replace_line_starting_with(jobpath + '/equilibration.mdp', 'couple-lambda1', ';couple-lambda1           = vdw-q') 
#        replace_line_starting_with(jobpath + '/production.mdp', 'couple-lambda1', ';couple-lambda1           = vdw-q')
#        replace_line_starting_with(jobpath + '/equilibration.mdp', 'couple-intramol', ';couple-intramol          = yes') 
#        replace_line_starting_with(jobpath + '/production.mdp', 'couple-intramol', ';couple-intramol          = yes')
            
    elif min_gamma == max_gamma:
        # Setting up vectors for restraining atoms
        lambdas = np.arange(min_lambda, max_lambda + 1, lambda_spacing)
        indicies = np.arange(0, (max_lambda - min_lambda) / lambda_spacing + 1, 1)
        gamma_points = np.zeros(len(indicies))
        if lambda_exponent < 0:
            lambda_points = 1 - ((max_lambda - lambdas) / max_lambda) ** abs(lambda_exponent)
        else:
            lambda_points = (lambdas / max_lambda) ** abs(lambda_exponent)
        init_lambda = np.where(current_lambda == lambdas)[0][0]

    # Setting the initial lambda state
    replace_line_starting_with(jobpath + '/equilibration.mdp', 'init_lambda_state', 'init_lambda_state        = ' + str(init_lambda))
    replace_line_starting_with(jobpath + '/production.mdp', 'init_lambda_state', 'init_lambda_state        = ' + str(init_lambda))

    # Write in the lambda states
    replace_line_starting_with(jobpath + '/equilibration.mdp', 'restraint_lambdas', 'restraint_lambdas   = ' + float_array_to_string(lambda_points))
    replace_line_starting_with(jobpath + '/equilibration.mdp', 'coul-lambdas', 'coul-lambdas   = ' + float_array_to_string(gamma_points))
    replace_line_starting_with(jobpath + '/equilibration.mdp', 'vdw-lambdas', 'vdw-lambdas   = ' + float_array_to_string(gamma_points))
    replace_line_starting_with(jobpath + '/production.mdp', 'restraint_lambdas', 'restraint_lambdas   = ' + float_array_to_string(lambda_points))
    replace_line_starting_with(jobpath + '/production.mdp', 'coul-lambdas', 'coul-lambdas   = ' + float_array_to_string(gamma_points))
    replace_line_starting_with(jobpath + '/production.mdp', 'vdw-lambdas', 'vdw-lambdas   = ' + float_array_to_string(gamma_points))

    if remove_bonded_interactions == True:
        replace_line_starting_with(jobpath + '/equilibration.mdp', ';bonded-lambdas', 'bonded-lambdas   = ' + float_array_to_string(gamma_points))
        replace_line_starting_with(jobpath + '/production.mdp', ';bonded-lambdas', 'bonded-lambdas   = ' + float_array_to_string(gamma_points))
    
   
def float_array_to_string(array_values):
    string = ''
    for i in array_values:
        string += str(np.around(i, 6)) + ' '
    return string

################################################################################
#### Specific functions for setup_molecule
################################################################################

def replace_string_in_text(fil, string, replacement_string):
    # Replaces the 'string' with 'replacement_string' in the provided file 'fil'
    s = open(fil).read()
    s = s.replace(str(string), str(replacement_string))
    f = open(fil, 'w')
    f.write(s)
    f.close()

def replace_line_starting_with(fil, string, line):
    # Replacing the line starting with 'string' with a new 'line' in a given file 'fil'
    with open(fil) as f:
        out = ''
        for l in f:
            if len(l.split()) > 0:
                if l.split()[0] == str(string):
                    l = str(line) + '\n'
            out += l
    with open(fil, 'w') as F:
        F.write(out)

def grofil_number_of_atoms(fil):
    # Determining the number of atoms in the gro file
    with open(fil) as f:
        number_of_atoms = int(f.read().split('\n')[1].split()[0])
    return number_of_atoms

def append_files(file_1, file_2):
    f1 = open(file_1)
    f2 = open(file_2)
    f = open('hold', 'w')
    f.write(f1.read())
    f.write(f2.read())
    subprocess.call(['mv', 'hold', file_1])

def setup_molecule(polymorph_num='p1', temperature=[], pressure=1, molecule='', number_of_molecules=0,
                   independent='same', equil_steps=100000, prod_steps=40000000, prodoutputs=20000, integrator='sd',
                   thermostat='nose-hoover', barostat='Parrinello-Rahman', cores=1, k_min=0, k_max=1000000, lambd=0,
                   min_lambda=0, max_lambda=0, lambda_spacing=100, lambda_exponent=2, gamma=100, min_gamma=100,
                   max_gamma=100, gamma_exponent=2, gamma_spacing=100, cutoff=8, potential='oplsaa',
                   simulation='gromacs', ensemble='NPT', jobpath='./', templatepath='', anneal_temp=400,
                   anneal_steps=10000, run_production=True, volume=-1, charge=0.1150, hinge='DefaultHinge', delta=0,
                   SigmaH=100, SigmaC=100, drude_k=100, submission_script='submit_cluster.slurm', remove_bonded_interactions=False):
    # =============================================================================================
    # ENSURE THAT INPUTS HAVE BEEN PROPERLY ENTERED
    # =============================================================================================
    # TEMPERATURE
    if temperature == -1:
        print("Invalid Temperature: ", temperature)
        sys.exit()

    # PRESSURE
    if pressure < 0:
        print("Invalid Pressure: ", pressure)
        sys.exit()

    ##POLYMORPH NUMBER
    if polymorph_num == 'gas':
        integrator = "sd"
        cutoff = "20"
        thermostat = "andersen"

    # INDEPENDENT MOLECULES
    if (independent == 'same') or (independent == number_of_molecules):
        independent = number_of_molecules
        independenthinge = str(number_of_molecules)
    else:
        independenthinge = str(independent) + 'ind'

    # SYSTEM
#    xyzfiles =$(ls ${TEMPLATEPATH} / *.xyz | grep "${MOLECULE}_" | grep "_${polymorph_num}_" | grep "_${MOLECULES}")
#    elif (simulation == 'tinker') and (xyzfiles == ''):
#        print("There are no available files in the runfiles directory for the combination: ")
#        print("Molecule: " + molecule)
#        print("Polymorph: " + polymorph_num)
#        print("Number: ", str(independent))
#        sys.exit()

    # THERMOSTAT
    if thermostat not in ["berendsen", "v-rescale", "andersen", "nose-hoover", "bussi"]:
        print("Unrecognized thermostat: " + thermostat)
        print("Supported thermostats: berendsen v-rescale andersen nose-hoover bussi")
        sys.exit()

    if "$barostat" == "parrinello-rahman":
        thermostat = "nose-hoover"

    # INTEGRATOR
    if integrator not in ["md", "md-vv", "sd"]:
        print("Unrecognized integrator: " + integrator)
        sys.exit()
    elif polymorph_num == "gas":
        integrator = "md-vv"

    if thermostat == "andersen":
        integrator = "md-vv"

    # BAROSTAT
    if barostat not in ["berendsen", "Parrinello-Rahman", "MTTK", "montecarlo"]:
        print("Unrecognized barostat: " + barostat)
        print("Supported barostats: berendsen Parrinello-Rahman MTTK montecarlo")
        sys.exit()

    # DRUDE SPRING CONSTANT
    if drude_k > 999:
        print("Drude spring constant too strong: ", drude_k)
        sys.exit()
    elif drude_k < 1:
        print("Drude spring constant too weak: ", drude_k)
        sys.exit()

    # LAMBDA POINT
    if (lambd < 0) or (lambd > 100):
        print("Invalid Lambda point: ", lambd)
        sys.exit()

    if (min_lambda < 0) or (max_lambda > 100) or (min_lambda > max_lambda):
        print("Minimum Lambda: ", min_lambda)
        print("Maximum Lambda: ", max_lambda)
        print("Is not a valid lambda range!")
        sys.exit()

    if lambda_spacing <= 0:
        print("Invalid Lambda Spacing: ", lambda_spacing)
        sys.exit()

#    if (lambda_exponent < 0) or (lambda_exponent > 4):
#        print("Invalid Lambda Exponent: ", lambda_exponent)
#        sys.exit()

    # GAMMA POINT
    if gamma_spacing < 0:
        print("Invalid Gambda Spacing: ", gamma_spacing)
        sys.exit()

#    if (gamma_exponent < 0) or (gamma_exponent > 4):
#        print("Invalid Gamma Exponent: ", gamma_exponent)
#        sys.exit()

    # SIGMAC
    if SigmaC < 0:
        print("Invalid SigmaC value: ", SigmaC)
        sys.exit()

    # SIGMAH
    if SigmaH < 0:
        print("Invalid SigmaC value: ", SigmaH)
        sys.exit()

    # CUTOFF RADIUS
    if cutoff < 0:
        print("Invalid Cutoff Radius: ", cutoff)
        sys.exit()

    # POTENTIAL
    potentiallist = ["oplsaa", "amber", "smirnoff", 
                     "gromos", "designedg", "oplsaatodesignedg", "designeda", "oplsaatodesigneda",
                     "amoeba09", "DMA", "PCA", "amoeba09todesa", "amoeba09restraint", "amoeba09interactions",
                     "amoeba09multinp", "amoeba09mononp", "amoeba09monoopls", "amoeba09opls", "day", "drude", "oplsaal"]
    valid = False
    for pot in potentiallist:
        if potential == pot:
            valid = True

    if valid == False:
        print("Unsupported potential: " + potential)
        print("Supported potentials: " + potentiallist)
        sys.exit()

    # SIMULATION PACKAGE
    if simulation not in ['gromacs', 'tinker']:
        print("Invalid Simulation Package: " + simulation)
        print("Supported Simulations: gromacs tinker")
        sys.exit()

    # ENSEMBLE
    if ensemble not in ['NVE', 'NVT', 'NPT']:
        print("Invalid Thermodynamic Ensemble: ",ensemble)
        print("Supported Ensembles: NVE, NVT, NPT")
        sys.exit()

    # =============================================================================================
    # FORMAT INPUTS FOR THE NAME OF THE JOB
    # =============================================================================================
    def number_to_string(value):
        if value < 10:
            return_string = '00' + str(value)
        elif value < 100:
            return_string = '0' + str(value)
        else:
            return_string = str(value)
        return return_string

    # Format the temperature name
    tempname = number_to_string(int(float(temperature)))

    # Format the number of molecules
    if number_of_molecules == independent:
        molnum = str(number_of_molecules)
    else:
        molnum = str(number_of_molecules) + '_' + str(independent) + 'ind'

    if pressure == '-1':
        pname = ''
    else:
        pname = '_' + number_to_string(pressure) + 'P'

    print(potential)

    # Format the potential
    if potential == 'oplsaa':
        potname = 'OPLS'
    elif potential == 'amber':
        potname = 'AMBER'
    elif potential == 'smirnoff':
        potname = 'SMIRNOFF'
    elif potential == 'gromos':
        potname = 'GROM'
    elif potential == 'designedg':
        potname = 'DESG'
    elif potential == 'oplsaatodesignedg':
        potname = 'OPLSDESG'
    elif potential == 'designeda':
        potname = 'DESA'
    elif potential == 'oplsaatodesigneda':
        potname = 'OPLSDESA'
    elif potential == 'DMA':
        potname = 'DMA'
    elif potential == 'PCA':
        potname = 'PCA'
    elif potential in ['amoeba09', 'amoeba09todesa', 'amoeba09restraint', 'amoeba09interactions']:
        potname = 'AMO'
    elif potential == 'amoeba09multinp':
        potname = 'MULTI'
    elif potential == 'amoeba09mononp':
        potname = 'MONONP'
    elif potential == 'amoeba09monoopls':
        potname = 'MONOOPLS'
    elif potential == 'amoeba09opls':
        potname = 'AMOOPLS'
    elif potential == 'day':
        potname = 'DAY'
    elif potential == 'drude':
        potname = 'DRUDE'
    elif potential == 'oplsaal':
        potname = 'OPLSL'

    # Format the simulation
    if simulation == 'gromacs':
        simname = 'GRO'
    elif simulation == 'tinker':
        simname = 'TIN'

    # Format the hinge if specified
    if hinge != 'DefaultHinge':
        hinge = '_' + hinge

    # Name the job
    gname = "pre_EQ"
    tname = "topology"

    # make the directory if it does not already exist
    print('Making Directory: ', jobpath,'...')
    subprocess.call(['mkdir', jobpath])

    # For PSCP setting if an NPT equilibration should be run prior to NPT
    NPT_equil = False
    if simulation == 'gromacs':
        # COPY OVER THE INITIAL GRO FILE
        print('Copying .gro file...')
        if polymorph_num == 'gas':
            subprocess.call(['cp', templatepath + '/' + molecule + '_gas_1.gro', jobpath + '/' + gname + '.gro'])
        else:
            if os.path.isfile(templatepath + '/' + molecule + '_' + polymorph_num + '_' + molnum + '_' + tempname + 'K_'
                              + pname + 'bar_' + potname + '.gro'):
                grofile = templatepath + '/' + molecule + '_' + polymorph_num + '_' + molnum + '_' + tempname + 'K_' + \
                          pname + 'bar_' + potname
            elif os.path.isfile(templatepath + '/' + molecule + '_' + polymorph_num + '_' + molnum + '_' + tempname +
                                'K_1bar_' + potname + '.gro'):
                 grofile = templatepath + '/' + molecule + '_' + polymorph_num + '_' + molnum + '_' + tempname + \
                          'K_1bar_' + potname
                 if ensemble == 'NVT':
                     NPT_equil = True
            elif os.path.isfile(templatepath + '/' + molecule + '_' + polymorph_num + '_' + molnum + '_REMD_' + potname
                                + '.gro'):
                 grofile = templatepath + '/' + molecule + '_' + polymorph_num + '_' + molnum + '_REMD_' + potname
                 if ensemble == 'NVT':
                     NPT_equil = True
            elif os.path.isfile(templatepath + '/' + molecule + '_' + polymorph_num + '_' + molnum + '_000K_' + potname
                                + '.gro'):
                 grofile = templatepath + '/' + molecule + '_' + polymorph_num + '_' + molnum + '_000K_' + potname
                 if ensemble == 'NVT':
                     NPT_equil = True
            elif os.path.isfile(templatepath + '/' + molecule + '_' + polymorph_num + '_' + molnum + '.gro'):
                 grofile = templatepath + '/' + molecule + '_' + polymorph_num + '_' + molnum
                 if ensemble == 'NVT':
                     NPT_equil = True
            else:
                print("There are no available files in the runfiles directory for the combination: ")
                print("Runfiles Directory: " + templatepath)
                print("Molecule: " + molecule)
                print("Polymorph: " + polymorph_num)
                print("Number: ", str(number_of_molecules))
                print("Independent: ", str(independent))
                sys.exit()

            # Copying over pre-equilibrated structure from templatepaths
            subprocess.call(['cp', grofile + '.gro', jobpath + '/' + gname + '.gro'])

            # Scale the box vectors if necessary
            # Determine the current volume of the unit cell
            if volume != -1:
                print('Resizing to ', volume * 0.01)
                run_resize_gro_individual(molecule, jobpath + '/' + gname + '.gro', volume)
        print('Using initial structure: ' + grofile)
        # COPY OVER THE RESTRAINT GRO FILE
        print('Copying restraint file...')
        if polymorph_num == 'gas':
            subprocess.call(['cp', templatepath + '/' + molecule + '_gas_1.gro', jobpath + '/restraint.gro'])
        else:
            if os.path.isfile(grofile + '_restraint.gro'):
                subprocess.call(['cp', grofile + '_restraint.gro', jobpath + '/restraint.gro'])
            else:
                subprocess.call(['cp', grofile + '.gro', jobpath + '/restraint.gro'])

        # Copy the template equilibration and production mdp file in to the new directory
        print('Copying .mdp files...')
        if polymorph_num == 'gas':
            subprocess.call(['cp', path + 'run_files/equilibration_gas.mdp', jobpath + '/equilibration.mdp'])
            subprocess.call(['cp', path + 'run_files/production_gas.mdp', jobpath + '/production.mdp'])
            subprocess.call(['cp', path + 'run_files/minimization.mdp', jobpath + '/minimization.mdp'])
            subprocess.call(['cp', path + 'run_files/relaxation.mdp', jobpath + '/relaxation.mdp'])
            subprocess.call(['cp', path + 'run_files/anneal.mdp', jobpath + '/anneal.mdp'])
        else:
            subprocess.call(['cp', path + 'run_files/equilibration.mdp', jobpath + '/equilibration.mdp'])
            subprocess.call(['cp', path + 'run_files/production.mdp', jobpath + '/production.mdp'])
            subprocess.call(['cp', path + 'run_files/anneal.mdp', jobpath + '/anneal.mdp'])
            subprocess.call(['cp', path + 'run_files/minimization.mdp', jobpath + '/minimization.mdp'])
            subprocess.call(['cp', path + 'run_files/relaxation.mdp', jobpath + '/relaxation.mdp'])
            if NPT_equil == True:
                subprocess.call(['cp', path + 'run_files/equilibration.mdp', jobpath + '/npt_equilibration.mdp'])
        replace_string_in_text(jobpath + '/minimization.mdp', 'MOLMOLMOLMOL', molecule)
        replace_string_in_text(jobpath + '/relaxation.mdp', 'MOLMOLMOLMOL', molecule)
        replace_string_in_text(jobpath + '/equilibration.mdp', 'MOLMOLMOLMOL', molecule)
        replace_string_in_text(jobpath + '/production.mdp', 'MOLMOLMOLMOL', molecule)
        replace_string_in_text(jobpath + '/anneal.mdp', 'MOLMOLMOLMOL', molecule)
        if NPT_equil == True:
            replace_string_in_text(jobpath + '/npt_equilibration.mdp', 'MOLMOLMOLMOL', molecule)

        print('Editing .mdp files...')
        # TEMPERATURE COUPLING
        if ensemble in ['NVT', 'NPT']:
            # Changing the thermostat in the .mdp files
            replace_line_starting_with(jobpath + '/equilibration.mdp', 'tcoupl', 'tcoupl = ' + thermostat)
            replace_line_starting_with(jobpath + '/production.mdp', 'tcoupl', 'tcoupl = ' + thermostat)
            replace_line_starting_with(jobpath + '/anneal.mdp', 'tcoupl', 'tcoupl = ' + thermostat)

            # Changing the reference temperature in the .mdp file
            replace_line_starting_with(jobpath + '/equilibration.mdp', 'ref_t', 'ref_t = ' + str(temperature))
            replace_line_starting_with(jobpath + '/production.mdp', 'ref_t', 'ref_t = ' + str(temperature))
            replace_line_starting_with(jobpath + '/anneal.mdp', 'ref_t', 'ref_t = ' + str(temperature))

            # Changing the gen_temp in the .mdp file
            replace_line_starting_with(jobpath + '/equilibration.mdp', 'gen_temp', 'gen_temp = ' + str(temperature))
            replace_line_starting_with(jobpath + '/production.mdp', 'gen_temp', 'gen_temp = ' + str(temperature))
            replace_line_starting_with(jobpath + '/anneal.mdp', 'gen_temp', 'gen_temp = ' + str(temperature))

            # Updating annealing setting in the anneal.mdp file
            replace_string_in_text(jobpath + '/anneal.mdp', 'END_ANNEAL_TIME', anneal_steps / 2000 - anneal_steps
                                   / 20000)
            replace_string_in_text(jobpath + '/anneal.mdp', 'STARTTEMP', anneal_temp)
            replace_string_in_text(jobpath + '/anneal.mdp', 'ENDTEMP', str(temperature))
            if NPT_equil == True:
                replace_line_starting_with(jobpath + '/npt_equilibration.mdp', 'tcoupl', 'tcoupl = ' + thermostat)
                replace_line_starting_with(jobpath + '/npt_equilibration.mdp', 'ref_t', 'ref_t = ' + str(temperature))
                replace_line_starting_with(jobpath + '/npt_equilibration.mdp', 'gen_temp', 'gen_temp = ' + str(temperature))
                replace_line_starting_with(jobpath + '/npt_equilibration.mdp', 'pcoupl', 'pcoupl = berendsen')
                replace_line_starting_with(jobpath + '/npt_equilibration.mdp', 'tau_p', 'tau_p = 1.0')

        # PRESSURE COUPLING
        if (ensemble == 'NPT'):
            # Updating the pressure coupling in the .mdp file
            replace_line_starting_with(jobpath + '/equilibration.mdp', 'pcoupl', 'pcoupl = berendsen')
            replace_line_starting_with(jobpath + '/production.mdp', 'pcoupl', 'pcoupl = ' + barostat)
            replace_line_starting_with(jobpath + '/anneal.mdp', 'pcoupl', 'pcoupl = berendsen')

            replace_line_starting_with(jobpath + '/equilibration.mdp', 'tau_p', 'tau_p = 1.0')
            replace_line_starting_with(jobpath + '/production.mdp', 'tau_p', 'tau_p = 10.0')
            replace_line_starting_with(jobpath + '/anneal.mdp', 'tau_p', 'tau_p = 1.0')

            # Making sure that we are sampling full anisotropy in the crystal lattice
            replace_line_starting_with(jobpath + '/equilibration.mdp', 'compressibility',
                                       'compressibility = 4.5e-5 4.5e-5 4.5e-5 4.5e-5 4.5e-5 4.5e-5')
            replace_line_starting_with(jobpath + '/production.mdp', 'compressibility',
                                       'compressibility = 4.5e-5 4.5e-5 4.5e-5 4.5e-5 4.5e-5 4.5e-5')
            replace_line_starting_with(jobpath + '/anneal.mdp', 'compressibility',
                                       'compressibility = 4.5e-5 4.5e-5 4.5e-5 4.5e-5 4.5e-5 4.5e-5')

        # CUTOFF RADIUS
        if cutoff != 9:
            coulombswitch = cutoff * 0.1 - 0.02
            rcoulomb = cutoff * 0.1
            vdwswitch = cutoff * 0.1 - 0.05
            rvdw = cutoff * 0.1

            # Updating the coulombswitch
            replace_line_starting_with(jobpath + '/equilibration.mdp', 'rcoulomb-switch', 'rcoulomb-switch = ' +
                                       str(coulombswitch))
            replace_line_starting_with(jobpath + '/production.mdp', 'rcoulomb-switch', 'rcoulomb-switch = ' +
                                       str(coulombswitch))
            replace_line_starting_with(jobpath + '/minimization.mdp', 'rcoulomb-switch', 'rcoulomb-switch = ' +
                                       str(coulombswitch))
            replace_line_starting_with(jobpath + '/relaxation.mdp', 'rcoulomb-switch', 'rcoulomb-switch = ' +
                                       str(coulombswitch))
            replace_line_starting_with(jobpath + '/anneal.mdp', 'rcoulomb-switch', 'rcoulomb-switch = ' +
                                       str(coulombswitch))

            replace_line_starting_with(jobpath + '/equilibration.mdp', 'rcoulomb', 'rcoulomb = ' + str(rcoulomb))
            replace_line_starting_with(jobpath + '/production.mdp', 'rcoulomb', 'rcoulomb = ' + str(rcoulomb))
            replace_line_starting_with(jobpath + '/minimization.mdp', 'rcoulomb', 'rcoulomb = ' + str(rcoulomb))
            replace_line_starting_with(jobpath + '/relaxation.mdp', 'rcoulomb', 'rcoulomb = ' + str(rcoulomb))
            replace_line_starting_with(jobpath + '/anneal.mdp', 'rcoulomb', 'rcoulomb = ' + str(rcoulomb))

            replace_line_starting_with(jobpath + '/equilibration.mdp', 'rvdw-switch', 'rvdw-switch = ' + str(vdwswitch))
            replace_line_starting_with(jobpath + '/production.mdp', 'rvdw-switch', 'rvdw-switch = ' + str(vdwswitch))
            replace_line_starting_with(jobpath + '/minimization.mdp', 'rvdw-switch', 'rvdw-switch = ' + str(vdwswitch))
            replace_line_starting_with(jobpath + '/relaxation.mdp', 'rvdw-switch', 'rvdw-switch = ' + str(vdwswitch))
            replace_line_starting_with(jobpath + '/anneal.mdp', 'rvdw-switch', 'rvdw-switch = ' + str(vdwswitch))

            replace_line_starting_with(jobpath + '/equilibration.mdp', 'rvdw', 'rvdw = ' + str(rvdw))
            replace_line_starting_with(jobpath + '/production.mdp', 'rvdw', 'rvdw = ' + str(rvdw))
            replace_line_starting_with(jobpath + '/minimization.mdp', 'rvdw', 'rvdw = ' + str(rvdw))
            replace_line_starting_with(jobpath + '/relaxation.mdp', 'rvdw', 'rvdw = ' + str(rvdw))
            replace_line_starting_with(jobpath + '/anneal.mdp', 'rvdw', 'rvdw = ' + str(rvdw))

            if NPT_equil == True:
                replace_line_starting_with(jobpath + '/npt_equilibration.mdp', 'rcoulomb-switch', 'rcoulomb-switch = ' + str(coulombswitch))
                replace_line_starting_with(jobpath + '/npt_equilibration.mdp', 'rcoulomb', 'rcoulomb = ' + str(rcoulomb))
                replace_line_starting_with(jobpath + '/npt_equilibration.mdp', 'rvdw-switch', 'rvdw-switch = ' + str(vdwswitch))
                replace_line_starting_with(jobpath + '/npt_equilibration.mdp', 'rvdw', 'rvdw = ' + str(rvdw))

            if cutoff > 10:
                replace_line_starting_with(jobpath + '/equilibration.mdp', 'rlist', 'rlist = ' + str(rvdw))
                replace_line_starting_with(jobpath + '/production.mdp', 'rlist', 'rlist = ' + str(rvdw))
                replace_line_starting_with(jobpath + '/minimization.mdp', 'rlist', 'rlist = ' + str(rvdw))
                replace_line_starting_with(jobpath + '/relaxation.mdp', 'rlist', 'rlist = ' + str(rvdw))
                replace_line_starting_with(jobpath + '/anneal.mdp', 'rlist', 'rlist = ' + str(rvdw))

                if NPT_equil == True:
                    replace_line_starting_with(jobpath + '/npt_equilibration.mdp', 'rlist', 'rlist = ' + str(rvdw))

        # TIMESTEPS
        replace_line_starting_with(jobpath + '/anneal.mdp', 'nsteps', 'nsteps = ' + str(anneal_steps))
        replace_line_starting_with(jobpath + '/equilibration.mdp', 'nsteps', 'nsteps = ' + str(equil_steps))
        replace_line_starting_with(jobpath + '/production.mdp', 'nsteps', 'nsteps = ' + str(prod_steps))
        if NPT_equil == True:
            replace_line_starting_with(jobpath + '/npt_equilibration.mdp', 'nsteps', 'nsteps = ' + str(equil_steps))

        # GENERATE VELOCITIES
        if anneal_steps > 0:
            replace_line_starting_with(jobpath + '/equilibration.mdp', 'gen_vel', 'gen_vel = no')
        if (anneal_steps > 0) or (equil_steps > 0):
            replace_line_starting_with(jobpath + '/production.mdp', 'gen_vel', 'gen_vel = no')

        # OUTPUT FREQUENCY
        replace_line_starting_with(jobpath + '/equilibration.mdp', 'nstlog', 'nstlog = ' + str(prodoutputs))
        replace_line_starting_with(jobpath + '/equilibration.mdp', 'nstenergy', 'nstenergy = ' +
                                   str(prodoutputs))
        replace_line_starting_with(jobpath + '/equilibration.mdp', 'nstxout', 'nstxout = ' +
                                   str(prodoutputs))

        replace_line_starting_with(jobpath + '/production.mdp', 'nstlog', 'nstlog = ' + str(prodoutputs))
        replace_line_starting_with(jobpath + '/production.mdp', 'nstenergy', 'nstenergy = ' +
                                   str(prodoutputs))
        replace_line_starting_with(jobpath + '/production.mdp', 'nstxout', 'nstxout = ' +
                                   str(prodoutputs))

        replace_line_starting_with(jobpath + '/anneal.mdp', 'nstlog', 'nstlog = ' + str(prodoutputs))
        replace_line_starting_with(jobpath + '/anneal.mdp', 'nstenergy', 'nstenergy = ' +
                                   str(prodoutputs))
        replace_line_starting_with(jobpath + '/anneal.mdp', 'nstxout', 'nstxout = ' +
                                   str(prodoutputs))
        if NPT_equil == True:
            replace_line_starting_with(jobpath + '/npt_equilibration.mdp', 'nstlog', 'nstlog = ' + str(prodoutputs))
            replace_line_starting_with(jobpath + '/npt_equilibration.mdp', 'nstenergy', 'nstenergy = ' +
                                       str(prodoutputs))
            replace_line_starting_with(jobpath + '/npt_equilibration.mdp', 'nstxout', 'nstxout = ' +
                                       str(prodoutputs))

        # INTEGRATOR
        replace_line_starting_with(jobpath + '/equilibration.mdp', 'integrator', 'integrator = ' + integrator)
        replace_line_starting_with(jobpath + '/production.mdp', 'integrator', 'integrator = ' + integrator)
        replace_line_starting_with(jobpath + '/anneal.mdp', 'integrator', 'integrator = ' + integrator)
        if NPT_equil == True:
            replace_line_starting_with(jobpath + '/npt_equilibration.mdp', 'integrator', 'integrator = ' + integrator)

        # FREE ENERGY PARAMETERS
        setup_mdp_lambdas(lambd, gamma, polymorph_num=polymorph_num, 
                          min_lambda=min_lambda, max_lambda=max_lambda, 
                          lambda_spacing=lambda_spacing, 
                          lambda_exponent=lambda_exponent, min_gamma=min_gamma, 
                          max_gamma=max_gamma, gamma_spacing=gamma_spacing, 
                          gamma_exponent=gamma_exponent, jobpath=jobpath, equil_output_frequency=prodoutputs,
                          prod_output_frequency=prodoutputs, remove_bonded_interactions=remove_bonded_interactions)

        #subprocess.call(['setup_mdpLambdas', '-L', str(lambd), '-W', str(min_lambda), '-S', str(max_lambda),
        #                 '-s', str(lambda_spacing), '-A', str(max_gamma), '-B', str(min_gamma), '-G', str(gamma),
        #                 '-g', str(gamma_spacing), '-f', str(lambda_exponent), '-F', str(gamma_exponent),
        #                 '-d', jobpath])

        # Copy over the molecule itp file and make the necessary modifications to the bond lengths, charges, and sigma values
        print('Copying itp file...')
        subprocess.call(['cp', templatepath + '/' + molecule + '_' + potential + '.itp', jobpath + '/molecule.itp'])
        subprocess.call(['cp', templatepath + '/' + molecule + '_' + potential + '.itp', './'])

        # Conduct any interpolations in the itp file
        if potential in ['oplsaatofakeg', 'oplsaatofakea']:
            interpolate_itp(jobpath + '/molecule.itp', delta)

        subprocess.call(['cp', templatepath + '/parameters.txt', jobpath + '/parameters.txt'])
        if potential in ['oplsaa', 'oplsaal', 'oplsaafakegd']:
            replace_string_in_text(jobpath + '/parameters.txt', 'CCCCC', '0.140')
            replace_string_in_text(jobpath + '/parameters.txt', 'HHHHH', '0.108')
        elif potential in ['gromos', 'oplsaafakeg', 'oplsaafakegb']:
            replace_string_in_text(jobpath + '/parameters.txt', 'CCCCC', '0.139')
            replace_string_in_text(jobpath + '/parameters.txt', 'HHHHH', '0.109')
        elif potential == 'drude':
            replace_string_in_text(jobpath + '/parameters.txt', 'CCCCC', '0.1375')
            replace_string_in_text(jobpath + '/parameters.txt', 'HHHHH', '0.108')
            replace_string_in_text(jobpath + '/parameters.txt', 'drude_k', np.around(drude_k * 4184, 8))
        elif potential == 'oplsaafakea':
            replace_string_in_text(jobpath + '/parameters.txt', 'CCCCC', '0.1382')
            replace_string_in_text(jobpath + '/parameters.txt', 'HHHHH', '0.1079')

        if charge != '':
            replace_string_in_text(jobpath + '/molecule.itp', '-0.115    12.011', '-' + str(np.around(charge, 3)) +
                                   '    12.011')
            replace_string_in_text(jobpath + '/molecule.itp', '0.115     1.008', str(np.around(charge, 3)) +
                                   '     1.008')
        replace_string_in_text(jobpath + '/parameters.txt', 'CCCHARGE', str(np.around(charge, 3)))

        # CREATE THE POSITION RESTRAINT ITP FILE
        c = subprocess.Popen(['echo', '0'], stdout=subprocess.PIPE)
        output = subprocess.check_output(['gmx', 'genrestr', '-f', jobpath + '/' + gname + '.gro',
                                                 '-o', jobpath + '/posre.itp',
                                                 '-fc', str(k_max), str(k_max), str(k_max), '-quiet'], stdin=c.stdout)
        c.wait()

        # Now lop off all but the first apermol + 4 lines
        atoms = grofil_number_of_atoms(jobpath + '/' + gname + '.gro')
        # echo "atoms: $atoms"
        apermol = int(atoms / number_of_molecules)
        with open(jobpath + '/posre.itp') as f:
            out = ''
            f = f.read().split('\n')
            for l in range(apermol + 4):
                if l > 3:
                    hold_posre = f[l].split()
                    out += '  ' + hold_posre[0] + '  ' + hold_posre[1] + '  0  0  0  ' + hold_posre[2] + '  ' + hold_posre[3] + '  ' + hold_posre[4] + '\n'
                else:
                    out += f[l] + '\n'

        with open(jobpath + '/restr.itp', 'w') as F:
            F.write(out)

        # COPY OVER THE INDEX FILE(for the force-averaging code)
        if number_of_molecules == independent:
            subprocess.call(['cp', templatepath + '/index.ndx', jobpath + '/index.ndx'])
        else:
            subprocess.call(['cp', templatepath + '/' + molecule + '_' + polymorph_num + '_' + molnum + '.ndx', jobpath
                             + '/index.ndx'])
            #replace_line_starting_with(jobpath + '/equilibration.mdp', 'symmetry-averaging', 'symmetry-averaging = yes')
            #replace_line_starting_with(jobpath + '/production.mdp', 'symmetry-averaging', 'symmetry-averaging = yes')
            #replace_line_starting_with(jobpath + '/minimization.mdp', 'symmetry-averaging', 'symmetry-averaging = yes')
            #replace_line_starting_with(jobpath + '/relaxation.mdp', 'symmetry-averaging', 'symmetry-averaging = yes')
            append_files(jobpath + '/equilibration.mdp', templatepath + '/' + molecule + '_' + str(independent) + 'ind_symmetry_groups.txt')
            append_files(jobpath + '/production.mdp', templatepath + '/' + molecule + '_' + str(independent) + 'ind_symmetry_groups.txt')
            append_files(jobpath + '/minimization.mdp', templatepath + '/' + molecule + '_' + str(independent) + 'ind_symmetry_groups.txt')
            append_files(jobpath + '/relaxation.mdp', templatepath + '/' + molecule + '_' + str(independent) + 'ind_symmetry_groups.txt')

        # Now edit the position restraint file to indicate a transformation over lambda space (handle all possible spacing cases)
        for i in range(3, 11):
            spaces = ''
            for j in range(i):
                spaces += ' '
            replace_string_in_text(jobpath + '/restr.itp', str(k_max) + spaces + str(k_max) + spaces + str(k_max),
                                   str(k_min) + ' ' + str(k_min) + ' ' + str(k_min) + ' ' + str(k_max) + ' ' +
                                   str(k_max) + ' ' + str(k_max))

        #subprocess.call(['mv', jobpath + '/restr.itp', jobpath + '/posre.itp'])
        subprocess.call(['cp', jobpath + '/restr.itp', jobpath + '/posre.itp'])

        # COPY OVER THE TOPOLOGY FILE
        print('Copying topology file...')
        subprocess.call(['cp', templatepath + '/topology.top', jobpath + '/' + tname + '.top'])

        # Edit the topology file based on the potential and system inputs
        replace_string_in_text(jobpath + '/topology.top', 'MOLMOLMOL', molecule)
        replace_string_in_text(jobpath + '/topology.top', 'NUMNUMNUM', number_of_molecules)

        if potential == 'gromos':
            replace_string_in_text(jobpath + '/' + tname + '.top', 'oplsaa.ff', 'gromos54a7.ff')
        elif potential in ['day', 'amber', 'smirnoff']:
            replace_string_in_text(jobpath + '/' + tname + '.top', '#include "oplsaa.ff/forcefield.itp"', '')
        elif potential == 'designedg':
            replace_string_in_text(jobpath + '/' + tname + '.top', 'oplsaa.ff', molecule + '_designedg.ff')
        elif potential == molecule + '_oplsaatodesignedg':
            replace_string_in_text(jobpath + '/' + tname + '.top', 'oplsaa.ff', molecule + '_oplsaatodesignedg')
        elif potential == 'designeda':
            replace_string_in_text(jobpath + '/' + tname + '.top', 'oplsaa.ff', 'designeda.ff')
        elif potential == 'oplsaatodesigneda':
            replace_string_in_text(jobpath + '/' + tname + '.top', 'oplsaa.ff', molecule + '_oplsaatodesigneda.ff')
        elif potential == 'designedd':
            replace_string_in_text(jobpath + '/' + tname + '.top', 'oplsaa.ff', molecule + '_designedd.ff')
        elif potential == 'drude':
            replace_string_in_text(jobpath + '/' + tname + '.top', 'oplsaa.ff', 'drude.ff')

        # Copy over local and cluster submission scripts
        print('Copying local and cluster submission scripts...')
        #subprocess.call(['cp', templatepath + '/submit_minimization_local.sh', jobpath + '/'])
        #subprocess.call(['cp', templatepath + '/submit_minimization.sh', jobpath + '/'])
        #subprocess.call(['cp', templatepath + '/submit_relaxation.sh', jobpath + '/'])
        print(path)
        subprocess.call(['cp', path + 'run_files/submit_local.sh', jobpath + '/'])
        subprocess.call(['cp', path + 'run_files/' + submission_script, jobpath + '/submit_cluster.slurm'])

        # if the number of Annealing steps is 0, skip the equilibration
        annealing = True
        if anneal_steps == 0:
            annealing = False
            print('Skipping Annealing...')

        # if the number of equilibration steps is 0, skip the equilibration
        equilibration = True
        if equil_steps == 0:
            equilibration = False
            print('Skipping Equilibration...')

        # if the number of production steps is 0, skip the equilibration
        production = True
        if (prod_steps == 0) or (run_production == False):
            production = False
            print('Skipping Production...')

        # If there is no force-averaging, remove the index file command
        indexing = True
        if number_of_molecules == independent:
            indexing = False
            #replace_string_in_text(jobpath + '/submit_minimization_local.sh', '-n index.ndx', '')
            #replace_string_in_text(jobpath + '/submit_minimization.sh', '-n index.ndx', '')

        # If we are not using the drude oscillator potential, remove the drude settings from the mdp file
        if potential != 'drude':
            for i in ['equilibration.mdp', 'production.mdp', 'anneal.mdp']:
                replace_line_starting_with(jobpath + '/' + i, 'drude', '')
                replace_line_starting_with(jobpath + '/' + i, 'drude-mode', '')
                replace_line_starting_with(jobpath + '/' + i, 'drude-hardwall', '')
                replace_line_starting_with(jobpath + '/' + i, 'drude-hyper', '')
                replace_line_starting_with(jobpath + '/' + i, 'drude-r', '')

        # IF THIS IS THE PSCP, DONT REWEIGHT, ENERGY MINIMIZE, OR PICK THE AVERAGE CONFIGURATION
        reweight = True
        if hinge in ['_L', '_G']:
            reweight = False

        # SET THE APPROPRIATE POST-SIMULATION REWEIGHTING
        if potential == 'oplsaa':
            potential_to_pass = 'oplsaa'
        elif potential in ['designeda', 'oplsaatodesigneda']:
            potential_to_pass = 'designeda'
        elif potential == 'designedd':
            potential_to_pass = 'designedd'
        elif potential == 'drude':
            "WARNING: DRUDE NO LONGER SET UP"
            sys.exit()
        elif potential in ['gromos', 'designedg']:
            potential_to_pass = 'gromos54a7 oplsaa designedg'
        else:
            replace_line_starting_with(jobpath + '/submit_cluster.slurm', 'convertjobtinker', '')

        replace_string_in_text(jobpath + '/submit_local.sh', 'rrrr', str(reweight))
        replace_string_in_text(jobpath + '/submit_local.sh', 'iiii', str(indexing))
#        replace_string_in_text(jobpath + '/submit_local.sh', 'uuuu', potential_to_pass)
        replace_string_in_text(jobpath + '/submit_local.sh', 'AAAA', str(annealing))
        replace_string_in_text(jobpath + '/submit_local.sh', 'EEEE', str(equilibration))
        replace_string_in_text(jobpath + '/submit_local.sh', 'PPPP', str(production))

        replace_string_in_text(jobpath + '/submit_cluster.slurm', 'rrrr', str(reweight))
        replace_string_in_text(jobpath + '/submit_cluster.slurm', 'iiii', str(indexing))
#        replace_string_in_text(jobpath + '/submit_cluster.slurm', 'uuuu', potential_to_pass)
        replace_string_in_text(jobpath + '/submit_cluster.slurm', 'aaaa', str(cores))
        replace_string_in_text(jobpath + '/submit_cluster.slurm', 'AAAA', str(annealing))
        replace_string_in_text(jobpath + '/submit_cluster.slurm', 'EEEE', str(equilibration))
        replace_string_in_text(jobpath + '/submit_cluster.slurm', 'PPPP', str(production))

    print('Create the job status file')
#echo "UNSUBMITTED" > ${JOBPATH} / jobstatus.txt
    print('Creating .initial directory')
    subprocess.call(['mkdir', jobpath + '/.initial'])
#for file in $(ls ${JOBPATH}); do
#    cp ${JOBPATH} /$file ${JOBPATH} /.initial /$file
#done
    print('Done!')


def re_setup(polymorph_num='p1', molecule='', number_of_molecules=0, 
             independent='same', potential='oplsaa', simulation='gromacs', 
             templatepath=''):
    # =============================================================================================
    # ENSURE THAT INPUTS HAVE BEEN PROPERLY ENTERED
    # =============================================================================================
    # POTENTIAL
    potentiallist = ["oplsaa", "amber", "smirnoff",
                     "gromos", "designedg", "oplsaatodesignedg", "designeda", "oplsaatodesigneda",
                     "amoeba09", "DMA", "PCA", "amoeba09todesa", "amoeba09restraint", "amoeba09interactions",
                     "amoeba09multinp", "amoeba09mononp", "amoeba09monoopls", "amoeba09opls", "day", "drude", "oplsaal"]

    valid = False
    for pot in potentiallist:
        if potential == pot:
            valid = True

    if valid == False:
        print("Unsupported potential: " + potential)
        print("Supported potentials: " + potentiallist)
        sys.exit()

    # SIMULATION PACKAGE
    if simulation not in ['gromacs', 'tinker']:
        print("Invalid Simulation Package: " + simulation)
        print("Supported Simulations: gromacs tinker")
        sys.exit()

    # =============================================================================================
    # FORMAT INPUTS FOR THE NAME OF THE JOB
    # =============================================================================================
    def number_to_string(value):
        if value < 10:
            return_string = '00' + str(value)
        elif value < 100:
            return_string = '0' + str(value)
        else:
            return_string = str(value)
        return return_string

    # INDEPENDENT MOLECULES
    if (independent == 'same') or (independent == number_of_molecules):
        independent = number_of_molecules
        independenthinge = str(number_of_molecules)
    else:
        independenthinge = str(independent) + 'ind'

    # Format the number of molecules
    if number_of_molecules == independent:
        molnum = str(number_of_molecules)
    else:
        molnum = str(number_of_molecules) + '_' + str(independent) + 'ind'

    # Format the potential
    if potential == 'oplsaa':
        potname = 'OPLS'
    elif potential == 'amber':
        potname = 'AMBER'
    elif potential == 'smirnoff':
        potname = 'SMIRNOFF'
    elif potential == 'gromos':
        potname = 'GROM'
    elif potential == 'designedg':
        potname = 'DESG'
    elif potential == 'oplsaatodesignedg':
        potname = 'OPLSDESG'
    elif potential == 'designeda':
        potname = 'DESA'
    elif potential == 'oplsaatodesigneda':
        potname = 'OPLSDESA'
    elif potential == 'DMA':
        potname = 'DMA'
    elif potential == 'PCA':
        potname = 'PCA'
    elif potential in ['amoeba09', 'amoeba09todesa', 'amoeba09restraint', 'amoeba09interactions']:
        potname = 'AMO'
    elif potential == 'amoeba09multinp':
        potname = 'MULTI'
    elif potential == 'amoeba09mononp':
        potname = 'MONONP'
    elif potential == 'amoeba09monoopls':
        potname = 'MONOOPLS'
    elif potential == 'amoeba09opls':
        potname = 'AMOOPLS'
    elif potential == 'day':
        potname = 'DAY'
    elif potential == 'drude':
        potname = 'DRUDE'
    elif potential == 'oplsaal':
        potname = 'OPLSL'

    # Format the simulation
    if simulation == 'gromacs':
        simname = 'GRO'
    elif simulation == 'tinker':
        simname = 'TIN'

    # Name the job
    gname = "pre_EQ"
    tname = "topology"

    if simulation == 'gromacs':
        # COPY OVER THE INITIAL GRO FILE
        print('Copying .gro file...')
        if polymorph_num == 'gas':
            subprocess.call(['cp', templatepath + '/' + molecule + '_gas_1.gro', 're_setup/'])
        else:
            if os.path.isfile(templatepath + '/' + molecule + '_' + polymorph_num + '_' + molnum + '_000K_' + potname
                                + '.gro'):
                grofile = templatepath + '/' + molecule + '_' + polymorph_num + '_' + molnum + '_000K_' + potname
            elif os.path.isfile(templatepath + '/' + molecule + '_' + polymorph_num + '_' + molnum + '.gro'):
                grofile = templatepath + '/' + molecule + '_' + polymorph_num + '_' + molnum
            else:
                print(templatepath + '/' + molecule + '_' + polymorph_num + '_' + molnum + '_000K_' + potname,
                      templatepath + '/' + molecule + '_' + polymorph_num + '_' + molnum)
                print("There are no available files in the runfiles directory for the combination: ")
                print("Runfiles Directory: " + templatepath)
                print("Molecule: " + molecule)
                print("Polymorph: " + polymorph_num)
                print("Number: ", str(number_of_molecules))
                print("Independent: ", str(independent))
                sys.exit()

            # Copying over pre-equilibrated structure from templatepaths
            subprocess.call(['cp', grofile + '.gro', 're_setup/'])

        print('Using initial structure: ' + grofile)
        # COPY OVER THE RESTRAINT GRO FILE
        print('Copying restraint file...')
        if polymorph_num == 'gas':
            subprocess.call(['cp', templatepath + '/' + molecule + '_gas_1.gro', 're_setup/'])
        else:
            if os.path.isfile(grofile + '_restraint.gro'):
                subprocess.call(['cp', grofile + '_restraint.gro', 're_setup/'])
            else:
                subprocess.call(['cp', grofile + '.gro', 're_setup/'])

        # Copy over the molecule itp file
        print('Copying itp file...')
        subprocess.call(['cp', templatepath + '/' + molecule + '_' + potential + '.itp', 're_setup/'])
        subprocess.call(['cp', templatepath + '/parameters.txt', 're_setup/'])

        # COPY OVER THE INDEX FILE(for the force-averaging code)
        if number_of_molecules == independent:
            subprocess.call(['cp', templatepath + '/index.ndx', 're_setup/'])
        else:
            subprocess.call(['cp', templatepath + '/' + molecule + '_' + polymorph_num + '_' + molnum + '.ndx', 're_setup/'])

        # COPY OVER THE TOPOLOGY FILE
        print('Copying topology file...')
        subprocess.call(['cp', templatepath + '/topology.top', 're_setup/'])


