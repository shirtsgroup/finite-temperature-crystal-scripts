#!/usr/bin/env python
import sys
import subprocess
import os
import numpy as np
from resize_gro_individual import run_resize_gro_individual
from interpolate_itp import interpolate_itp

path = os.path.realpath(__file__).strip('Run_LatticeDynamics.py')

def setup_interactions(polymorph_num=['p1', 'p2' 'p3'], molecule="benzene", max_gamma=100, min_gamma=0, gamma_spacing=10,
                       gamma_exponent=2, lambd=100, PSCP_temperature=200, pressure=1, number_of_molecules=72,
                       independent="same", gamma_equil_steps=100000, gamma_prod_steps=1000000, integrator="sd",
                       thermostat="nose-hoover", cores=7, k_max=1000, cutoff="8", potential="oplsaa",
                       hinge="DefaultHinge"):
    #NSA: The checks Eric initially had, should they go here or integrated elsewhere? (Are they more general then just this?)
    # Going through the setup for each polymorph
    for i in polymorph_num:
         gamma = min_gamma
         subprocess.call(['mkdir', i + '/interactions'])  # making a directory for turning off the interactions
         
         # Setting up the directory for each gamma value between the minimum and maximum gamma
         while gamma <= max_gamma:
             setup_molecule(i, PSCP_temperature, pressure, molecule, number_of_molecules, independent, 
                            gamma_equil_steps, gamma_prod_steps, integrator, thermostat, lambd, lambd, lambd, 
                            gamma_exponent, cores, 'NVT', k_max, max_gamma, min_gamma, gamma_spacing, gamma, cutoff, 
                            potential, hinge, i + '/interactions/' + str(gamma))
             gamma += gamma_spacing

def setup_restraints(polymorph_num=['p1', 'p2', 'p3'], molecule="benzene", min_lambda=0, max_lambda=100, 
                     lambda_spacing=5, lambda_exponent=4, gamma=100, PSCP_temperature=200, pressure=1, 
                     number_of_molecules=72, independent="same", lambda_equil_steps=1000000, lambda_prod_steps=10000000,
                     integrator="sd", thermostat="nose-hoover", cores=7, k_min=0, k_max=1000, cutoff="8", 
                     potential="oplsaa", hinge="DefaultHinge"):
# NSA: The checks Eric initially had, should they go here or integrated elsewhere? (Are they more general then just this?)
    # Going through the setup for each polymorph
    for i in polymorph_num:
        lambd = min_lambda
        subprocess.call(['mkdir', i + '/restraints'])  # making a directory for turning off the interactions

        # Setting up the directory for each lambda value between the minimum and maximum lambda
        while lambd <= max_lambda:
            setup_molecule()
            lambd += lambda_spacing


def setup_mdp_lambdas(current_lambda, current_gamma, polymorph_num='all', min_lambda=0, max_lambda=100, lambda_spacing=-1, lambda_exponent=2, min_gamma=0, max_gamma=100,
                      gamma_spacing=-1, gamma_exponent=2, jobpath='DefaultPath'):
    # Python script to automatically set up the mdp files to take on multiple lambda values
    # Original script written in bash by: Eric Dybeck on 09/12/2014
    # Converted to python by: Nate Abraham on 01/23/2019

    # Ensure that the parameters are properly entered
    # Lambda
    if (min_lambda < 0) or (max_lambda > 100) or (min_lambda > max_lambda):
        print('Minimum Lambda: ', min_lambda)
        print('Maximum Lambda: ', max_lambda)
        print('Is not a valid lambda range!')
        sys.exit()

    if lambda_spacing <= 0:
        print('Invalid Lambda Spacing: ', lambda_spacing)
        sys.exit()

    if (min_lambda == max_lambda) and (current_lambda!= max_lambda):
        print('Minimum Lambda: ', min_lambda, ' Maximum Lambda: ', max_lambda, ' and Lambda: ', current_lambda,
              ' are not the same!')
        sys.exit()

    if not 1 <= lambda_exponent <= 4:
        print('Invalid Lambda Exponent: ', lambda_exponent)
        sys.exit()

    # Gamma
    if (min_gamma < 0) or (max_gamma > 100) or (min_gamma > max_gamma):
        print('Minimum Gamma: ', min_gamma)
        print('Maximum Gamma: ', max_gamma)
        print('Is not a valid lambda range!')
        sys.exit()

    if gamma_spacing <= 0:
        print('Invalid Gamma Spacing: ', gamma_spacing)
        sys.exit()

    if (min_gamma == max_gamma) and (current_gamma != max_gamma):
        print('Minimum Gamma: ', min_gamma, ' Maximum Gamma: ', max_gamma, ' and Gamma: ', current_gamma,
              ' are not the same!')
        sys.exit()

    if not 1 <= gamma_exponent <= 4:
        print('Invalid Gamma Exponent: ', gamma_exponent)
        sys.exit()

    #JOBPATH
    if jobpath == 'DefaultPath':
        print('Enter the job path!')
        sys.exit()

    #If we have no harmonic restraints and full interactions, no need to proceed further
    if (min_lambda == max_lambda) and (max_lambda == 0) and (min_gamma == max_gamma) and (max_gamma == 100):
        print('No lambda values added.')
        sys.exit()

    #If we are adding harmonic restraints, we should not be changing gamma (and vice versa)
    if (min_lambda != max_lambda) and (min_gamma != max_gamma):
        print('Harmonic restraints and Interactions changing simultaneously!!')
        sys.exit(

    # Change the free energy setting from 'no' to 'yes' and the output from 0 to nstxout
    replace_line_starting_with(jobpath + '/equilibration.mdp', 'free_energy', 'free_energy = yes')
    replace_line_starting_with(jobpath + '/production.mdp', 'free_energy', 'free_energy = yes')
    #log_equil =$(less ${JOBPATH} / equilibration.mdp | grep -m 1 'nstenergy' | awk '{print $3}')
    #log_prod =$(less ${JOBPATH} / production.mdp | grep -m 1 'nstenergy' | awk '{print $3}')
    replace_line_starting_with(jobpath + '/equilibration.mdp', 'nstdhdl', 'nstdhdl = ' + str(log_equil))
    replace_line_starting_with(jobpath + '/production.mdp', 'nstdhdl', 'nstdhdl = ' + str(log_equil))

    lambda_indicies = ''
    lambda_vector = ''
    gamma_indicies = ''
    gamma_vector = ''

    # If harmonic restraints are being added, loop over all lambda points and set up the lambda vector in the mdp files
    if min_lambda != max_lambda:
        raw_lambda = min_lambda
        lambda_hold = 0.
        i = 0
        while raw_lambda < max_lambda:
            if i < 10:
                lambda_indicies = lambda_indicies + '    ' + str(i) + '    '
            else:
                lambda_indicies = lambda_indicies + '    ' + str(i) + '   '

            lambda_1 = np.around(raw_lambda ** lambda_exponent / (max_lambda ** (lambda_exponent -1))), 6)
            lambda_2 = np.around(lambda_1 / 100., 6)
            gamma_1 = np.around(max_gamma ** 2/ 100, 6)
            gamma_2 = np.around(gamma_1 / 100., 6)
            lambda_vector =








"""
       Lambda_vect="${Lambda_vect}${lambda} "
       Gamma_vect="${Gamma_vect}${gamma} "

       if[$LAMBDA == $RawLambda]; then
       sed -i "/init_lambda_state/c init_lambda_state        = $i    ;Which state are we sampling from?" ${JOBPATH} / equilibration.mdp
       sed -i "/init_lambda_state/c init_lambda_state        = $i    ;Which state are we sampling from?" ${JOBPATH} / production.mdp
       fi

       let "RawLambda=$RawLambda+$LSPACING"
       let "i=i+1"
       done

       # Catch the final temperature off-by-one exception
       Lambda_indicies="${Lambda_indicies}   ${i}    "
       lambda =$(echo "x=$MAXLAMBDA*0.010000; if(x==0) print \"0.00000\"; if(x<1 && x>0) print 0; x" | bc)
       Gamma=$(echo "scale=6; ($MAXGAMMA^2) / (100)" | bc)
       gamma=$(echo "x=$Gamma*0.010000; if(x==0) print \"0.00000\"; if(x<1 && x>0) print 0; x" | bc)
       Lambda_vect="${Lambda_vect}${lambda} "
       Gamma_vect="${Gamma_vect}${gamma} "

       if[$LAMBDA == $MAXLAMBDA]; then
       sed -i "/init_lambda_state/c init_lambda_state        = $i    ;Which state are we sampling from?" ${JOBPATH} / equilibration.mdp
       sed -i "/init_lambda_state/c init_lambda_state        = $i    ;Which state are we sampling from?" ${JOBPATH} / production.mdp
       fi

       # Set the coupling specifications for turning the interactions off
       if["$GAMMA" == "0"]; then
       sed -i "/couple-lambda0/c couple-lambda0           = none" ${JOBPATH} / equilibration.mdp
       sed -i "/couple-lambda1/c couple-lambda1           = vdw-q" ${JOBPATH} / equilibration.mdp
       sed -i "/couple-intramol/c couple-intramol          = yes" ${JOBPATH} / equilibration.mdp
       sed -i "/couple-lambda0/c couple-lambda0           = none" ${JOBPATH} / production.mdp
       sed -i "/couple-lambda1/c couple-lambda1           = vdw-q" ${JOBPATH} / production.mdp
       sed -i "/couple-intramol/c couple-intramol          = yes" ${JOBPATH} / production.mdp
       fi
       Lambda_indicies=";Lambda Indicies         = ${Lambda_indicies}"
       Lambda_vect="restraint_lambdas        = ${Lambda_vect}"
       Gamma_vect1="coul-lambdas             = ${Gamma_vect}"
       Gamma_vect2="vdw-lambdas              = ${Gamma_vect}"
       Gamma_vect3="bonded-lambdas           = "

       else
       # If the interactions are changing, loop over all gamma points and set up the gamma vector in the mdp files
       RawGamma=$MINGAMMA
       Lambda_indicies=""
       Lambda_vect=""
       Gamma_indicies=""
       Gamma_vect=""
       # Gamma_vect_bond=""
       gamma=0.0
       i=0
       while["$RawGamma" -lt "$MAXGAMMA"]; do

       if["$i" -lt "10"]; then
       Lambda_indicies="${Lambda_indicies}    ${i}    "
       else
       Lambda_indicies="${Lambda_indicies}    ${i}   "
       fi

       Gamma=$(echo "scale=6; ($RawGamma^$EXPONENT) / ($MAXGAMMA^($EXPONENT-1))" | bc)
       lambda =$(echo "x=$MAXLAMBDA*0.010000; if(x==0) print \"0.00000\"; if(x<1 && x>0) print 0; x" | bc)
       gamma=$(echo "x=$Gamma*0.010000; if(x==0) print \"0.00000\"; if(x<1 && x>0) print 0; x" | bc)
       # gamma_bond=$(echo "x=1.000000-$gamma; if(x==0) print \"0.00000\"; if(x<1 && x>0) print 0; x" | bc)
       # gamma="1.000000"
       Lambda_vect="${Lambda_vect}${lambda} "
       Gamma_vect="${Gamma_vect}${gamma} "
       Gamma_vect_bond="${Gamma_vect_bond}${gamma_bond} "

       if[$GAMMA == $RawGamma] & &[$MINLAMBDA == $MAXLAMBDA]; then
       sed -i "/init_lambda_state/c init_lambda_state        = $i    ;Which state are we sampling from?" ${JOBPATH} / equilibration.mdp
       sed -i "/init_lambda_state/c init_lambda_state        = $i    ;Which state are we sampling from?" ${JOBPATH} / production.mdp
       fi

       let "RawGamma=$RawGamma+$GSPACING"
       let "i=i+1"
       done


       # Catch the final gamma off-by-one exception
       Lambda_indicies="${Lambda_indicies}    ${i}   "
       lambda =$(echo "x=$MAXLAMBDA*0.010000; if(x==0) print \"0.00000\"; if(x<1 && x>0) print 0; x" | bc)
       gamma=$(echo "x=$MAXGAMMA*0.010000; if(x==0) print \"0.00000\"; if(x<1 && x>0) print 0; x" | bc)
       # gamma_bond=$(echo "x=1.000000-$gamma; if(x==0) print \"0.00000\"; if(x<1 && x>0) print 0; x" | bc)
       # gamma="1.000000"
       Lambda_vect="${Lambda_vect}${lambda} "
       Gamma_vect="${Gamma_vect}${gamma} "
       Gamma_vect_bond="${Gamma_vect_bond}${gamma_bond} "

       if[$GAMMA == $MAXGAMMA] & &[$MINLAMBDA == $MAXLAMBDA]; then
       sed -i "/init_lambda_state/c init_lambda_state        = $i    ;Which state are we sampling from?" ${JOBPATH} / equilibration.mdp
       sed -i "/init_lambda_state/c init_lambda_state        = $i    ;Which state are we sampling from?" ${JOBPATH} / production.mdp
       fi

       Lambda_indicies=";Lambda Indicies         = ${Lambda_indicies}"
       Lambda_vect="restraint_lambdas        = ${Lambda_vect}"
       Gamma_vect1="coul-lambdas             = ${Gamma_vect}"
       Gamma_vect2="vdw-lambdas              = ${Gamma_vect}"
       Gamma_vect3="bonded-lambdas           = ${Gamma_vect_bond}"

       # Set the coupling specifications for turning the interactions off
       sed -i "/couple-lambda0/c couple-lambda0           = none" ${JOBPATH} / equilibration.mdp
       sed -i "/couple-lambda1/c couple-lambda1           = vdw-q" ${JOBPATH} / equilibration.mdp
       sed -i "/couple-intramol/c couple-intramol          = yes" ${JOBPATH} / equilibration.mdp
       sed -i "/couple-lambda0/c couple-lambda0           = none" ${JOBPATH} / production.mdp
       sed -i "/couple-lambda1/c couple-lambda1           = vdw-q" ${JOBPATH} / production.mdp
       sed -i "/couple-intramol/c couple-intramol          = yes" ${JOBPATH} / production.mdp
       fi

       # Now replace the free-energy section with the new strings
       sed -i "/;Lambda Indicies/c ${Lambda_indicies}" ${JOBPATH} / equilibration.mdp
       sed -i "/restraint_lambdas/c ${Lambda_vect}" ${JOBPATH} / equilibration.mdp
       sed -i "/coul-lambdas/c ${Gamma_vect1}" ${JOBPATH} / equilibration.mdp
       sed -i "/vdw-lambdas/c ${Gamma_vect2}" ${JOBPATH} / equilibration.mdp
       sed -i "/bonded-lambdas/c ${Gamma_vect3}" ${JOBPATH} / equilibration.mdp
       sed -i "/;Lambda Indicies/c ${Lambda_indicies}" ${JOBPATH} / production.mdp
       sed -i "/restraint_lambdas/c ${Lambda_vect}" ${JOBPATH} / production.mdp
       sed -i "/coul-lambdas/c ${Gamma_vect1}" ${JOBPATH} / production.mdp
       sed -i "/vdw-lambdas/c ${Gamma_vect2}" ${JOBPATH} / production.mdp
       sed -i "/bonded-lambdas/c ${Gamma_vect3}" ${JOBPATH} / production.mdp
"""
################################################################################
#### Specific functions for setup_molecule
################################################################################

def replace_string_in_text(fil, string, replacement_string):
    # Replaces the 'string' with 'replacement_string' in the provided file 'fil'
    s = open(fil).read()
    s = s.replace(string, replacement_string)
    f = open(fil, 'w')
    f.write(s)
    f.close()

def replace_line_starting_with(fil, string, line):
    # Replacing the line starting with 'string' with a new 'line' in a given file 'fil'
    with open(fil) as f:
        out = ''
        for l in f:
            if len(l.split()) > 0:
                if l.split()[0] == string:
                    l = line
            out += l
    with open(fil, 'w') as F:
        F.write(out)

def grofil_number_of_atoms(fil):
    # Determining the number of atoms in the gro file
    with open(fil) as f:
        number_of_atoms = int(f[1].split()[0])
    return number_of_atoms

def append_files(file_1, file_2):
    f1 = open(file_1)
    f2 = open(file_2)
    f = open('hold', 'w')
    f.write(f1.read())
    f.write(f2.read())
    subprocess.call(['mv', 'hold', file_1])


def setup_molecule(polymorph_num, temperature, pressure, molecule, number_of_molecules, independent, equil_steps,
                   eqoutputs, prod_steps, prodoutputs, integrator, thermostat, barostat, cores, k_min, k_max, lambd, min_lambda,
                   max_lambda, lambda_spacing, lambda_exponent, gamma, min_gamma, max_gamma, gamma_exponent, gamma_spacing, cutoff, potential,
                   simulation, ensemble, jobpath, templatepath, anneal_temp, anneal_steps, run_production, volume,
                   charge, hinge, delta, SigmaH, SigmaC, drude_k):

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

    if (exponent < 0) or (exponent > 4):
        print("Invalid Exponent: ", exponent)
        sys.exit()

    # GAMMA POINT
    if gamma_spacing < 0:
        print("Invalid Gambda Spacing: ", gamma_spacing)
        sys.exit()

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
    potentiallist = ["oplsaa", "gromos", "designedg", "oplsaatodesignedg", "designeda", "oplsaatodesigneda",
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
    tempname = number_to_string(temperature)

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
    print('Making Directory: ${JOBPATH} ...')
    subprocess.call(['mkdir', jobpath])

    # OUTPUT FREQUENCY
    equil_output_frequency = equil_steps / 100  # The 100 was eqoutput
    equil_trr_output_frequency = equil_steps / 100  # The 100 was emtroutputs
    prod_output_frequency = prod_steps / prodoutputs
    prod_trr_output_frequency = prod_steps / prodoutputs
    anneal_output_frequency = anneal_steps / 10000  # The 10000 could be a user specified variable
    anneal_trr_output_frequency = anneal_steps / 10000
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
            elif os.path.isfile(templatepath + '/' + molecule + '_' + polymorph_num + '_' + molnum + '_000K_' + potname
                                + '.gro'):
                grofile = templatepath + '/' + molecule + '_' + polymorph_num + '_' + molnum + '_000K_' + potname
            elif os.path.isfile(templatepath + '/' + molecule + '_' + polymorph_num + '_' + molnum + '.gro'):
                grofile = templatepath + '/' + molecule + '_' + polymorph_num + '_' + molnum
            else:
                print("There are no available files in the runfiles directory for the combination: ")
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
            subprocess.call(['cp', templatepath + '/equilibration_gas.mdp', jobpath + '/equilibration.mdp'])
            subprocess.call(['cp', templatepath + '/production_gas.mdp', jobpath + '/production.mdp'])
            subprocess.call(['cp', templatepath + '/minimization.mdp', jobpath + '/minimization.mdp'])
            subprocess.call(['cp', templatepath + '/relaxation.mdp', jobpath + '/relaxation.mdp'])
            subprocess.call(['cp', templatepath + '/anneal_gas.mdp', jobpath + '/anneal.mdp'])
        else:
            subprocess.call(['cp', templatepath + '/equilibration.mdp', jobpath + '/equilibration.mdp'])
            subprocess.call(['cp', templatepath + '/production.mdp', jobpath + '/production.mdp'])
            subprocess.call(['cp', templatepath + '/minimization.mdp', jobpath + '/minimization.mdp'])
            subprocess.call(['cp', templatepath + '/relaxation.mdp', jobpath + '/relaxation.mdp'])
            subprocess.call(['cp', templatepath + '/anneal.mdp', jobpath + '/anneal.mdp'])
        replace_string_in_text(jobpath + '/minimization.mdp', 'MOLMOLMOLMOL', molecule)
        replace_string_in_text(jobpath + '/relaxation.mdp', 'MOLMOLMOLMOL', molecule)
        replace_string_in_text(jobpath + '/equilibration.mdp', 'MOLMOLMOLMOL', molecule)
        replace_string_in_text(jobpath + '/production.mdp', 'MOLMOLMOLMOL', molecule)
        replace_string_in_text(jobpath + '/anneal.mdp', 'MOLMOLMOLMOL', molecule)

        print('Editing .mdp files...')
        # TEMPERATURE COUPLING
        if ensemble in ['NVT', 'NPT']:
            # Changing the thermostat in the .mdp files
            replace_line_starting_with(jobpath + '/equilibration.mdp', 'tcoupl', 'tcoupl = ' + thermostat)
            replace_line_starting_with(jobpath + '/prodcution.mdp', 'tcoupl', 'tcoupl = ' + thermostat)
            replace_line_starting_with(jobpath + '/anneal.mdp', 'tcoupl', 'tcoupl = ' + thermostat)

            # Changing the reference temperature in the .mdp file
            replace_line_starting_with(jobpath + '/equilibration.mdp', 'ref_t', 'ref_t = ' + temperature)
            replace_line_starting_with(jobpath + '/production.mdp', 'ref_t', 'ref_t = ' + temperature)
            replace_line_starting_with(jobpath + '/anneal.mdp', 'ref_t', 'ref_t = ' + temperature)

            # Changing the gen_temp in the .mdp file
            replace_line_starting_with(jobpath + '/equilibration.mdp', 'gen_temp', 'gen_temp = ' + temperature)
            replace_line_starting_with(jobpath + '/production.mdp', 'gen_temp', 'gen_temp = ' + temperature)
            replace_line_starting_with(jobpath + '/anneal.mdp', 'gen_temp', 'gen_temp = ' + temperature)

            # Updating annealing setting in the anneal.mdp file
            replace_string_in_text(jobpath + '/anneal.mdp', 'END_ANNEAL_TIME', anneal_steps / 2000 - anneal_steps
                                   / 20000)
            replace_string_in_text(jobpath + '/anneal.mdp', 'STARTTEMP', anneal_temp)
            replace_string_in_text(jobpath + '/anneal.mdp', 'ENDTEMP', temperature)

        # PRESSURE COUPLING
        if ensemble == 'NPT':
            # Updating the pressure coupling in the .mdp file
            replace_line_starting_with(jobpath + '/equilibration.mdp', 'pcoupl', 'pcoupl = berendsen')
            replace_line_starting_with(jobpath + '/production.mdp', 'pcoupl', 'pcoupl = ' + barostat)
            replace_line_starting_with(jobpath + '/anneal.mdp', 'pcoupl', 'pcoupl = berendsen')

            replace_line_starting_with(jobpath + 'equilibration.mdp', 'tau_p', 'tau_p = 1.0')
            replace_line_starting_with(jobpath + 'production.mdp', 'tau_p', 'tau_p = 10.0')
            replace_line_starting_with(jobpath + 'anneal.mdp', 'tau_p', 'tau_p = 1.0')

            # Making sure that we are sampling full anisotropy in the crystal lattice
            replace_line_starting_with(jobpath + 'equilibration.mdp', 'compressibility',
                                       'compressibility = 4.5e-5 4.5e-5 4.5e-5 4.5e-5 4.5e-5 4.5e-5')
            replace_line_starting_with(jobpath + 'production.mdp', 'compressibility',
                                       'compressibility = 4.5e-5 4.5e-5 4.5e-5 4.5e-5 4.5e-5 4.5e-5')
            replace_line_starting_with(jobpath + 'anneal.mdp', 'compressibility',
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
            if cutoff > 10:
                replace_line_starting_with(jobpath + '/equilibration.mdp', 'rlist', 'rlist = ' + str(rvdw))
                replace_line_starting_with(jobpath + '/production.mdp', 'rlist', 'rlist = ' + str(rvdw))
                replace_line_starting_with(jobpath + '/minimization.mdp', 'rlist', 'rlist = ' + str(rvdw))
                replace_line_starting_with(jobpath + '/relaxation.mdp', 'rlist', 'rlist = ' + str(rvdw))
                replace_line_starting_with(jobpath + '/anneal.mdp', 'rlist', 'rlist = ' + str(rvdw))

        # TIMESTEPS
        replace_line_starting_with(jobpath + '/equilibration.mdp', 'nsteps', 'nsteps = ' + str(equil_steps))
        replace_line_starting_with(jobpath + '/production.mdp', 'nsteps', 'nsteps = ' + str(prod_steps))
        replace_line_starting_with(jobpath + '/anneal.mdp', 'nsteps', 'nsteps = ' + str(anneal_steps))

        # OUTPUT FREQUENCY
        replace_line_starting_with(jobpath + '/equilibraiton.mdp', 'nstlog', 'nstlog = ' + str(equil_output_frequency))
        replace_line_starting_with(jobpath + '/equilibraiton.mdp', 'nstenergy', 'nstenergy = ' +
                                   str(equil_output_frequency))
        replace_line_starting_with(jobpath + '/equilibraiton.mdp', 'nstxout', 'nstxout = ' +
                                   str(equil_trr_output_frequency))
        replace_line_starting_with(jobpath + '/equilibraiton.mdp', 'nstxout-compressed', 'nstxout-compressed = ' +
                                   str(equil_trr_output_frequency))

        replace_line_starting_with(jobpath + '/production.mdp', 'nstlog', 'nstlog = ' + str(prod_output_frequency))
        replace_line_starting_with(jobpath + '/production.mdp', 'nstenergy', 'nstenergy = ' +
                                   str(prod_output_frequency))
        replace_line_starting_with(jobpath + '/production.mdp', 'nstxout', 'nstxout = ' +
                                   str(prod_trr_output_frequency))
        replace_line_starting_with(jobpath + '/production.mdp', 'nstxout-compressed', 'nstxout-compressed = ' +
                                   str(prod_trr_output_frequency))

        replace_line_starting_with(jobpath + '/anneal.mdp', 'nstlog', 'nstlog = ' + str(prod_output_frequency))
        replace_line_starting_with(jobpath + '/anneal.mdp', 'nstenergy', 'nstenergy = ' +
                                   str(prod_output_frequency))
        replace_line_starting_with(jobpath + '/anneal.mdp', 'nstxout', 'nstxout = ' +
                                   str(prod_trr_output_frequency))
        replace_line_starting_with(jobpath + '/anneal.mdp', 'nstxout-compressed', 'nstxout-compressed = ' +
                                   str(prod_trr_output_frequency))

        # INTEGRATOR
        replace_line_starting_with(jobpath + '/equilibraton.mdp', 'integrator', 'integrator = ' + integrator)
        replace_line_starting_with(jobpath + '/production.mdp', 'integrator', 'integrator = ' + integrator)
        replace_line_starting_with(jobpath + '/anneal.mdp', 'integrator', 'integrator = ' + integrator)

        # FREE ENERGY PARAMETERS
#setup_mdpLambdas - L $Lambda - W $MINLAMBDA - S $MAXLAMBDA - s $LSPACING - A $MAXGAMMA - B $MINGAMMA - G $GAMMA - g $GSPACING - f $EXPONENT - d $JOBPATH


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
            replace_string_in_text(jobpath + 'molecule.itp', '-0.115    12.011', '-' + str(np.around(charge, 3)) +
                                   '    12.011')
            replace_string_in_text(jobpath + 'molecule.itp', '0.115     1.008', str(np.around(charge, 3)) +
                                   '     1.008')
        replace_string_in_text(jobpath + '/parameters.txt', 'CCCHARGE', str(np.around(charge, 3)))

        # CREATE THE POSITION RESTRAINT ITP FILE
        c = subprocess.Popen(['gmx', 'genrestr', '-f', jobpath + '/' + gname + '.gro'
                                                 '-o', jobpath + '/posre.itp',
                                                 '-fc', str(k_max) + ' ' + str(k_max) + ' ' + str(k_max),
                                                 '-quiet'])
        c.communicate(b'0 0\n')
        c.wait()

        # Now lop off all but the first apermol + 4 lines
        atoms = grofil_number_of_atoms(jobpath + '/' + gname + '.gro')
        # echo "atoms: $atoms"
        apermol = atoms / number_of_molecules
        with open(jobpath + '/posre.itp') as f:
            out = ''
            for l in range(apermol + 4):
                out += f[l]
        with open(jobpath + '/restr.itp', 'w') as F:
            F.write(out)

        # COPY OVER THE INDEX FILE(for the force-averaging code)
        if number_of_molecules == independent:
            subprocess.call(['cp', templatepath + '/index.ndx', jobpath + '/index.ndx'])
        else:
            subprocess.call(['cp', templatepath + '/' + molecule + '_' + polymorph_num + '_' + molnum + '.ndx', jobpath
                             + '/index.ndx'])
            replace_line_starting_with(jobpath + '/equilibration.mdp', 'symmetry-averaging', 'symmetry-averaging = yes')
            replace_line_starting_with(jobpath + '/production.mdp', 'symmetry-averaging', 'symmetry-averaging = yes')
            replace_line_starting_with(jobpath + '/minimization.mdp', 'symmetry-averaging', 'symmetry-averaging = yes')
            replace_line_starting_with(jobpath + '/relaxation.mdp', 'symmetry-averaging', 'symmetry-averaging = yes')
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

        subprocess.call(['mv', jobpath + '/restr.itp', jobpath + '/posre.itp'])

        # COPY OVER THE TOPOLOGY FILE
        print('Copying topology file...')
        subprocess.call(['cp', templatepath + '/topology.top', jobpath + '/' + tname + '.top'])

        # Edit the topology file based on the potential and system inputs
        replace_string_in_text(jobpath + '/topology.top', 'MOLMOLMOL', molecule)
        replace_string_in_text(jobpath + '/topology.top', 'NUMNUMNUM', number_of_molecules)

        if potential == 'gromos':
            replace_string_in_text(jobpath + '/' + tname + '.top', 'oplsaa.ff', 'gromos54a7.ff')
        elif potential == 'day':
            replace_string_in_text(jobpath + '/' + tname + '.top', 'oplsaa.ff', '')
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
        subprocess.call(['cp', path + 'submit_local.sh', jobpath + '/'])
        subprocess.call(['cp', path + 'submit_cluster.sh', jobpath + '/'])

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
            replace_string_in_text(jobpath + '/submit_minimization_local.sh', '-n index.ndx', '')
            replace_string_in_text(jobpath + '/submit_minimization.sh', '-n index.ndx', '')

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
        replace_string_in_text(jobpath + '/submit_local.sh', 'uuuu', potential_to_pass)
        replace_string_in_text(jobpath + '/submit_local.sh', 'AAAA', str(annealing))
        replace_string_in_text(jobpath + '/submit_local.sh', 'EEEE', str(equilibration))
        replace_string_in_text(jobpath + '/submit_local.sh', 'PPPP', str(production))

        replace_string_in_text(jobpath + '/submit_cluster.slurm', 'rrrr', str(reweight))
        replace_string_in_text(jobpath + '/submit_cluster.slurm', 'iiii', str(indexing))
        replace_string_in_text(jobpath + '/submit_cluster.slurm', 'uuuu', potential_to_pass)
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

"""
NSA: SEE LINE 114 before implimenting this.
    elif simulation == 'tinker':
        # Copy over the key and xyz file
        if polymorph_num == 'gas':
            subprocess.call(['cp', templatepath + '/' + molecule + '_gas_1.key', jobpath + '/keyfile.key'])
            subprocess.call(['cp', templatepath + '/' + molecule + '_gas_1.xyz', jobpath + '/molecule.xyz'])
            xyzfile = molecule + '_gas_1'
        else:
            if os.path.isfile(templatepath + '/' + molecule + '_' + polymorph_num + '_' + str(number_of_molecules) +
                              '_' + tempname + 'K_' + str(pressure) + 'bar_' + potname + '.xyz'):
                xyzfile = molecule + '_' + polymorph_num + '_' + str(number_of_molecules) + '_' + tempname + 'K_' + \
                          str(pressure) + 'bar_' + potname
            elif os.path.isfile(templatepath + '/' + molecule + '_' + polymorph_num + '_' + str(number_of_molecules) +
                                '_' + tempname + 'K_1bar_' + potname + '.xyz'):
                xyzfile = molecule + '_' + polymorph_num + '_' + str(number_of_molecules) + '_' + tempname + 'K_1bar_' \
                          + potname + '.xyz'
            elif os.path.isfile(templatepath + '/' + molecule + '_' + polymorph_num + '_' + str(number_of_molecules) +
                                '_000K_' + potname + '.xyz'):
                xyzfile = molecule + '_' + polymorph_num + '_' + str(number_of_molecules) + '_000K_' + potname + '.xyz'
            else:
                xyzfile = molecule + '_' + polymorph_num + '_' + str(number_of_molecules)
            subprocess.call(['cp', templatepath + '/' + xyzfile + '.xyz', jobpath + '/molecule.xyz'])
            subprocess.call(['cp', templatepath + '/' + xyzfile + '.key', jobpath + '/keyfile.key'])

        print('Using initial structure: ', xyzfile)

        # Modify the key and parameter file to match the potential
        if potential in ['oplsaa', 'day']:
            replace_string_in_text(jobpath + '/keyfile.key', 'amoeba09.prm', 'oplsaa.prm')
            # ${SCRIPT_LOC} / addatomtypesindividual - p oplsaa - x ${JOBPATH} / molecule.xyz - M $MOLECULE
        elif potential == 'oplsaafake':
            replace_string_in_text(jobpath + '/keyfile.key', 'amoeba09.prm', 'oplsaafake.prm')
            #${SCRIPT_LOC} / addatomtypesindividual - p oplsaa - x ${JOBPATH} / molecule.xyz - M $MOLECULE
        elif potential == 'amoeba09todesa':
            subprocess.call(['cp', templatepath + '/amoeba09todesa.prm', jobpath + '/param.prm'])
            #python ${SCRIPT_LOC} / interpolate_itp.py - f ${JOBPATH} / param.prm - d $DELTA
            #sed - i "s/PARAMETERS.*/PARAMETERS param.prm/g" ${JOBPATH} / keyfile.key
            #${SCRIPT_LOC} / addatomtypesindividual - p amoeba09 - x ${JOBPATH} / molecule.xyz - M $MOLECULE
        else:
            pass
            #sed - i "s/\/params.*/\/params\/${POTENTIAL}.prm/g" ${JOBPATH} / keyfile.key
            #${SCRIPT_LOC} / addatomtypesindividual - p amoeba09 - x ${JOBPATH} / molecule.xyz - M $MOLECULE

        # Copy over the job specs file
        subprocess.call(['cp', templatepath + '/job_specs.txt', jobpath + '/'])
        replace_string_in_text(jobpath + '/job_specs.txt', 'SSSS', str(prod_steps))
        replace_string_in_text(jobpath + '/job_specs.txt', 'dtdt', str(0.5))


cp ${TEMPLATEPATH} / job_specs.txt ${JOBPATH} / job_specs.txt
sed - i
"s/SSSS/$prod_steps/g" ${JOBPATH} / job_specs.txt
sed - i
"s/dtdt/0.5/g" ${JOBPATH} / job_specs.txt
output_frac =$(echo "scale=8;x=1.0/${prod_output_frequency}; if(x<1) print 0; x" | bc)
output_dt =$(echo "scale=4;x=${prod_output_frequency}*0.5*0.001; if(x<1) print 0; x" | bc)
sed - i
"s/XXXX/${output_dt}/g" ${JOBPATH} / job_specs.txt
sed - i
"s/TTTT/$TEMP/g" ${JOBPATH} / job_specs.txt
if [ $ENSEMBLE == "NPT"]; then
sed - i
"s/EEEE/4   /g" ${JOBPATH} / job_specs.txt
sed - i
"s/(NVT)/(NPT)/g" ${JOBPATH} / job_specs.txt
sed - i
"s/PPPP/$PRESSURE/g" ${JOBPATH} / job_specs.txt
else
sed - i
"s/EEEE/2   /g" ${JOBPATH} / job_specs.txt
sed - i
"s/PPPP//g" ${JOBPATH} / job_specs.txt
fi

# Add on the relevant supercell expansions to make the gromacs supercell
less ${TEMPLATEPATH} /${MOLNAME}
_${polymorph_num}
_${MOLECULES}
_expansion.txt >> ${JOBPATH} / job_specs.txt

# Copy over local and cluster submission scripts
echo
"Copying local and cluster submission scripts..."
cp ${TEMPLATEPATH} / submit_localitc_tinker.sh ${JOBPATH} / submit_localitc.sh
cp ${TEMPLATEPATH} / submit_cluster_tinker.slurm ${JOBPATH} / submit_cluster.slurm

# Modify the local and cluster submission scripts
sed - i
"s/SSSS/$prod_steps/g" ${JOBPATH} / submit_cluster.slurm
sed - i
"s/SSSS/$prod_steps/g" ${JOBPATH} / submit_localitc.sh
sed - i
"s/dtdt/0.5/g" ${JOBPATH} / submit_cluster.slurm
sed - i
"s/dtdt/0.5/g" ${JOBPATH} / submit_localitc.sh
output_dt =$(echo "scale=4;x=${prod_output_frequency}*0.5*0.001; if(x<1) print 0; x" | bc)
sed - i
"s/XXXX/${output_dt}/g" ${JOBPATH} / submit_cluster.slurm
sed - i
"s/XXXX/${output_dt}/g" ${JOBPATH} / submit_localitc.sh
sed - i
"s/TTTT/$TEMP/g" ${JOBPATH} / submit_cluster.slurm
sed - i
"s/TTTT/$TEMP/g" ${JOBPATH} / submit_localitc.sh
if [ $ENSEMBLE == "NPT"]; then
sed - i
"s/EEEE/4   /g" ${JOBPATH} / submit_cluster.slurm
sed - i
"s/EEEE/4   /g" ${JOBPATH} / submit_localitc.sh
sed - i
"s/PPPP/$PRESSURE/g" ${JOBPATH} / submit_cluster.slurm
sed - i
"s/PPPP/$PRESSURE/g" ${JOBPATH} / submit_localitc.sh
else
sed - i
"s/EEEE/2   /g" ${JOBPATH} / submit_cluster.slurm
sed - i
"s/EEEE/2   /g" ${JOBPATH} / submit_localitc.sh
sed - i
"s/PPPP//g" ${JOBPATH} / submit_cluster.slurm
sed - i
"s/PPPP//g" ${JOBPATH} / submit_localitc.sh
fi

# Handle any special exceptions for the interpolated potential
if ["$POTENTIAL" == "amoeba09todesa"]; then
sed - i
"s/reweightjobtinker.*/reweightjobtinker -s tinker -u amoeba09todesa -d 10/g" ${JOBPATH} / submit_cluster.slurm
sed - i
"s/reweightjobtinker.*/reweightjobtinker -s tinker -u amoeba09todesa -d 10/g" ${JOBPATH} / submit_localitc.sh
elif ["$POTENTIAL" == "amoeba09restraint"];
then
sed - i
"s/reweightjobtinker.*/reweightjobtinker -s tinker -u amoeba09 -d 10 -e 1/g" ${JOBPATH} / submit_cluster.slurm
sed - i
"s/reweightjobtinker.*/reweightjobtinker -s tinker -u amoeba09 -d 10 -e 1/g" ${JOBPATH} / submit_localitc.sh
elif ["$POTENTIAL" == "amoeba09interactions"];
then
sed - i
"s/reweightjobtinker.*/reweightjobtinker -s tinker -u amoeba09interactions -d 10 -e 1/g" ${
                                                                                          JOBPATH} / submit_cluster.slurm
sed - i
"s/reweightjobtinker.*/reweightjobtinker -s tinker -u amoeba09interactions -d 10 -e 1/g" ${JOBPATH} / submit_localitc.sh
elif ["$POTENTIAL" == "oplsaa"];
then
sed - i
"s/reweightjobgromacs.*/reweightjobgromacs -s tinker -u oplsaa/g" ${JOBPATH} / submit_cluster.slurm
sed - i
"s/reweightjobgromacs.*/reweightjobgromacs -s tinker -u oplsaa/g" ${JOBPATH} / submit_localitc.sh
sed - i
"s/reweightjobtinker.*/reweightjobtinker -s tinker -u oplsaa/g" ${JOBPATH} / submit_cluster.slurm
sed - i
"s/reweightjobtinker.*/reweightjobtinker -s tinker -u oplsaa/g" ${JOBPATH} / submit_localitc.sh
fi

# Set the number of threads
sed - i
"s/-nt 1/-nt ${cores}/g" ${JOBPATH} / submit_localitc.sh
sed - i
"s/-nt 1/-nt ${cores}/g" ${JOBPATH} / submit_cluster.slurm
sed - i
"s/ntasks=1/ntasks=${cores}/g" ${JOBPATH} / submit_cluster.slurm

fi

"""




