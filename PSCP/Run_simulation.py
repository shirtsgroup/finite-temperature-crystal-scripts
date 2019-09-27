#!/usr/bin/env python

# Script created by Nathan S. Abraham
# Created on March 30, 2019
# Shirts Group; Dept. of Chemical & Biological Engineering; University of Colorado Boulder
# There was an equivalent script created that was written in bash, now converted to python

import sys
import os
import subprocess
import numpy as np
import argparse
path = os.path.dirname(os.path.abspath(__file__))



parser = argparse.ArgumentParser(description='')
parser.add_argument('-c', '--num_cores', dest='num_cores', default=1, help='Number of cores to run simulation on')
parser.add_argument('-r', '--reweight', dest='reweight', default=True, help='Run reweighting of energy, True or False')
parser.add_argument('-e', '--equilibration', dest='equilibration', default=True, help='Should equilibration be run, True or False')
parser.add_argument('-i', '--indexing', dest='indexing', default=False, help='Adds on the -n to a gromacs simulation for indexing')
parser.add_argument('-u', '--potential', dest='potential', default='oplsaa', help='Potential to used in the simulation')
parser.add_argument('-p', '--production', dest='production', default=True, help='Should a production simulation be run, True or False')
parser.add_argument('-a', '--anneal', dest='anneal', default=False, help='Should annealing be run, True or False')
args = parser.parse_args()

def check_bool(val):
    if val == 'True':
        return True
    elif val == 'False':
        return False
    else:
        print("Option not valid")
        sys.exit()

# setting up flags for re-running simulations
npt_run = True
anneal_run = True
eq_run = True
if os.path.isfile('PROD.trr'):
    npt_run = False
    anneal_run = False
    eq_run = False
elif os.path.isfile('EQ.trr'):
    npt_run = False
    anneal_run = False
elif os.path.isfile('ANNEAL.trr'):
    npt_run = False

# Run a NPT simulation if interactions are being turned off and not starting from the right T and P
if os.path.isfile('npt_equilibration.mdp') and npt_run:
    # Running annealing of crystal
    subprocess.call(['mv', 'pre_EQ.gro', 'pre_NPT.gro'])

    if check_bool(args.indexing):
        # For force averaging
        subprocess.call('gmx grompp -f npt_equilibration.mdp -c pre_NPT.gro -r restraint.gro -p topology.top -o NPT_equil.tpr -n index.ndx -maxwarn 10', shell=True)
    else:
        subprocess.call('gmx grompp -f npt_equilibration.mdp -c pre_NPT.gro -r restraint.gro -p topology.top -o NPT_equil.tpr -maxwarn 10', shell=True)

    # Running the annealing simulation
    subprocess.call("gmx mdrun -nt " + str(args.num_cores) + " -v -deffnm NPT_equil -dhdl dhdl_NPT_equil", shell=True)

    # Extracing the final frame from the checkpoint file
    subprocess.call("echo '0' | gmx trjconv -f NPT_equil.cpt -s NPT_equil.tpr -o pre_EQ.gro -pbc whole -ndec 12", shell=True)

# Run temperature annealing
if check_bool(args.anneal) and anneal_run:
    # Running annealing of crystal
    subprocess.call(['cp', 'pre_EQ.gro', 'pre_ANNEAL.gro'])

    if check_bool(args.indexing):
        # For force averaging
        subprocess.call('gmx grompp -f anneal.mdp -c pre_ANNEAL.gro -r restraint.gro -p topology.top -o ANNEAL.tpr -n index.ndx -maxwarn 10', shell=True)
    else:
        subprocess.call('gmx grompp -f anneal.mdp -c pre_ANNEAL.gro -r restraint.gro -p topology.top -o ANNEAL.tpr -maxwarn 10', shell=True)

    # Running the annealing simulation
    subprocess.call("gmx mdrun -nt " + str(args.num_cores) + " -v -deffnm ANNEAL -dhdl dhdl_ANNEAL", shell=True)

    # Extracing the final frame from the checkpoint file
    subprocess.call("echo '0' | gmx trjconv -f ANNEAL.cpt -s ANNEAL.tpr -o ANNEAL.gro -pbc whole -ndec 12 -vel yes", shell=True)

# Run equilibration
if check_bool(args.equilibration) and eq_run:
    if check_bool(args.anneal):
        # If annealing was run, copy over file
        subprocess.call(['cp', 'ANNEAL.gro', 'pre_EQ.gro'])

    if check_bool(args.indexing):
        # For force averaging
        subprocess.call('gmx grompp -f equilibration.mdp -c pre_EQ.gro -r restraint.gro -p topology.top -o EQ.tpr -n index.ndx -maxwarn 10', shell=True)
    else:
        subprocess.call('gmx grompp -f equilibration.mdp -c pre_EQ.gro -r restraint.gro -p topology.top -o EQ.tpr -maxwarn 10', shell=True)
 
    # Running the annealing simulation
    subprocess.call("gmx mdrun -nt " + str(args.num_cores) + " -v -deffnm EQ -dhdl dhdl_EQ", shell=True)

    # Extracing the final frame from the checkpoint file
    subprocess.call("echo '0' | gmx trjconv -f EQ.cpt -s EQ.tpr -o EQ.gro -pbc whole -ndec 12 -vel yes", shell=True)

# Running a production run
if  check_bool(args.production):
    # Determine the starting structure
    if (check_bool(args.anneal) == True) and (check_bool(args.equilibration) == False):
        starting_structure = 'ANNEAL.gro'
    elif (args.equilibration == False):
        starting_structure = 'EQ.gro'
    else:
        starting_structure = 'pre_EQ.gro'
    
    # Setting up production run
    subprocess.call("gmx grompp -f production.mdp -c " + starting_structure + " -r restraint.gro -p topology.top -o PROD.tpr  -maxwarn 10", shell=True)

    # Runnning production run
    subprocess.call("gmx mdrun -nt " + str(args.num_cores) + " -v -deffnm PROD -dhdl dhdl_PROD -cpi PROD.cpt", shell=True)

    # Extracing the final frame from the checkpoint file
    subprocess.call("echo '0' | gmx trjconv -f PROD.cpt -s PROD.tpr -o PROD.gro -pbc whole -ndec 12 -vel yes", shell=True)

#    # Calculate the configuration energies and ensemble averages
#    subprocess.call(path + '/run-scripts/crunchjobenergy')
#
#    # Relax the molecule and energy minimize the intramolecular interactions to create the restraint files
#    if check_bool(args.indexing) == True:
#        index = 'true'
#    else:
#        index = 'false'
#
#    subprocess.call(path + '/run-scripts/relax_molecule -i' + index)
#
#    if check_bool(args.reweight) == True:
#        # Reweight the job into the other potentials if there are no restraints
#        subprocess.call([path + '/run-scripts/reweightjobgromacs', '-s', 'gromacs', '-u', args.potential])
#        if args.potential == 'designeda':
#            subprocess.call([path + '/run-scripts/reweightjobtinker', '-s', 'gromacs', '-u', 'amoeba09todesa', '-d', '10'])



