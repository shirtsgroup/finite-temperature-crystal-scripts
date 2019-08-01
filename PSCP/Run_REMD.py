#!/usr/bin/env python

# Script created by Nathan S. Abraham
# Created on March 28, 2019
# Shirts Group; Dept. of Chemical & Biological Engineering; University of Colorado Boulder
# Current versions of gromacs cannot re-run a REMD simulation from a checkpoint file if the simulation times out


# PRIOR SETUP: This script assumes replica are set-up by directories named {0,1,2,....N}
# ISSUES: This script has yet to be implimented to extend simulations. The easiest way to impliment this is to change 
#          nsteps in the last created .mdp file by Delta(nsteps) that you want to increase for the total simulations

import sys
import os
import subprocess
import numpy as np
import itertools as it

def run_REMD_gromacs(num_replica, output_string='PROD', mdp='production.mdp', initial_structure='EQ.gro', 
                     restraints=None, topology='topology.top', replex=50, nex=True):
    # Initializing some parameters
    if restraints == None:
        restraint_string = ['', '']
    else:
        restraint_string = ['-r', restraints]

    if nex == True:
        nex_string = ['-nex', str(num_replica ** 3)]
    else:
        nex_string = ['', '']

    # Generate the tpr file for all replicas
    replica_dirs = '{'
    for i in range(num_replica):
        dir_string = str(i) + '/'
        subprocess.call(['gmx', 'grompp', '-f', dir_string + mdp, '-c', dir_string + initial_structure, 
                         restraint_string[0], dir_string + restraint_string[1], '-p', dir_string + topology, 
                         '-o', dir_string + output_string + '.tpr', '-maxwarn', '10'])
        subprocess.call(['rm', 'mdout.mdp'])
        
        if i == 0:
            replica_dirs += str(i)
        else:
            replica_dirs += ',' + str(i)
    replica_dirs += '}'  

    # Running replica exchange
    subprocess.call('mpirun -np ' + str(num_replica) + ' gmx_mpi mdrun -v -multidir ' + replica_dirs + 
                     ' -replex ' +  str(replex) + ' ' +  nex_string[0] + ' ' + nex_string[1] + ' -deffnm ' + output_string + ' -cpt 0.5', shell=True)


def new_mdp(last_t, mdp, new_mdp):
    # Opening the old mdp file
    with open(mdp) as f:
        mdp_data = np.array(list(it.zip_longest(*[lines.split() for lines in f], fillvalue=' '))).T

    # determinign the dt value in the mdp file
    stepsize = float(mdp_data[np.where(mdp_data[:, 0] == 'dt')[0][0], 2])

    # Determineing the previous value of nsteps in mdp
    nsteps_old = int(mdp_data[np.where(mdp_data[:, 0] == 'nsteps')[0][0], 2])

    # Updating the nsteps in new_mdp to keep the total simulation length constant
    mdp_data[np.where(mdp_data[:, 0] == 'nsteps')[0][0], 2] = int(nsteps_old - last_t / stepsize)

    # Turning off gen_vel in the new_mdp since we are starting from a configuration at the last checkpoint file
    mdp_data[np.where(mdp_data[:, 0] == 'gen_vel')[0][0], 2] = 'no'

    # Re-writing the new_mdp file    
    string_mdp = ''
    for i in range(len(mdp_data[:, 0])):
        for j in range(len(mdp_data[i, :])):
            string_mdp += mdp_data[i, j] + ' ' 
        string_mdp = string_mdp + '\n'

    with open(new_mdp, 'w') as file_out:
        file_out.write(string_mdp)

def get_checkpoint_time(cpt_file):
    # runs gmx_mpi check on the input cpt_file to determine where the new simulation must start from
    hold = subprocess.check_output('gmx_mpi check -f ' + cpt_file, shell=True, stderr=subprocess.STDOUT).decode("utf-8").split('\n')
    for i in hold:
        if len(i.split()) > 1:
            if i.split()[0] == 'Last':
                last_t = float(i.split()[4])
    return(last_t)


def run_REMD_gromacs_multi(num_replica, output_string='PROD', mdp='production.mdp', initial_structure='EQ.gro',
                           restraints=None, topology='topology.top', replex=50, nex=True):
    # Determine where the number of simulations run prior
    files = os.listdir('0/')
    check = False
    count = 0
    while check == False:
        if output_string + '_' + str(count) + '.cpt' not in files:
            check = True
        else:
            count += 1

    # Running the first simulation if it has yet to be run
    if count == 0:
        # Making a copy of the mdp file for this run
        for i in range(num_replica):
            subprocess.call(['cp', str(i) + '/' + mdp, str(i) + '/' + output_string + '_0' + '.mdp'])

        # Running the simulation
        run_REMD_gromacs(num_replica, output_string=output_string + '_0', mdp=output_string + '_0' + '.mdp', initial_structure=initial_structure,
                         restraints=restraints, topology=topology, replex=replex, nex=nex)

        # If finished in time, moving the first run edr and trr files to the desired output
        for i in range(num_replica):
            subprocess.call(['cp', str(i) + '/' + output_string + '_0.trr', str(i) + '/' + output_string + '.trr'])
            subprocess.call(['cp', str(i) + '/' + output_string + '_0.edr', str(i) + '/' + output_string + '.edr'])

    else:
        # Getting the file names
        last_file = output_string + '_' + str(count - 1)
        new_file = output_string + '_' + str(count)

        # Check to see where the last checkpoint file left off
        last_t = get_checkpoint_time(' 0/' + last_file + '.cpt')

        for i in range(num_replica):
            # Re-writing the mdp file
            new_mdp(last_t, str(i) + '/' + last_file + '.mdp', str(i) + '/' + new_file + '.mdp')

            # Creating a new gro file with the trajectories form the last checkpoint
            c = subprocess.Popen(['echo', '0'], stdout=subprocess.PIPE)
            output = subprocess.check_output(['gmx', 'trjconv', 
                                              '-f', str(i) + '/' + last_file + '.cpt', 
                                              '-s', str(i) + '/' + last_file + '.tpr', 
                                              '-o', str(i) + '/' + new_file + '.gro', 
                                              '-vel'], stdin=c.stdout)
            c.wait()

        # Running the simulation
        run_REMD_gromacs(num_replica, output_string=new_file, mdp=new_file + '.mdp', initial_structure=new_file + '.gro',
                         restraints=restraints, topology=topology, replex=replex, nex=nex) 

        # Changing the start time of the trr and edr files
        file_list = [output_string + '_0']
        for j in range(1, count + 1):
            file_list.append(output_string + '_' + str(j))
            start_t = get_checkpoint_time(' 0/' + output_string + '_' + str(j - 1) + '.cpt')
            for i in range(num_replica):
                dirr = str(i) + '/' + output_string + '_' + str(j) 
                subprocess.call("gmx trjconv -f " + dirr + ".trr -t0 " + str(start_t) + " -o " + dirr + ".trr", shell=True)
                c = subprocess.Popen(['echo', str(start_t)], stdout=subprocess.PIPE)
                output = subprocess.check_output(['gmx', 'eneconv',
                                                  '-f', dirr + '.edr',
                                                  '-o', dirr + '.edr',
                                                  '-settime', 'yes'], stdin=c.stdout)
                c.wait()

        # Concatinating the edr and trr files together
        for i in range(num_replica):
            trr_list = ''
            edr_list = ''
            for j in file_list:
                trr_list += str(i) + '/' + j + '.trr '
                edr_list += str(i) + '/' + j + '.edr '
            dirr = str(i) + '/' + output_string 
            subprocess.call("gmx trjcat -f " + trr_list + "-o " + dirr + ".trr -cat", shell=True)
            subprocess.call("gmx eneconv -f " + edr_list + "-o " + dirr + ".edr", shell=True)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-n', '--num_replica', dest='num_replica', help="Number of replica in the simulation. WARNING: This code assumes the replica are in directories that go from {0,1,2,...N} ex. -n N+1")
    parser.add_argument('-o', '--output_string', dest='output_string', default='PROD', help='The string for the output trr and edr files ex. PROD.trr and PROD.edr')
    parser.add_argument('-m', '--mdp', dest='mdp', default='production.mdp', help='GROMACS .mdp file in each of the replcia directories')
    #parser.add_argument('-i', '--initial_structure', dest='initial_structure', default='EQ.gro', help='GROMACS .gro file to start from in each of the replica directories')
    parser.add_argument('-r', '--restraints', dest='restraints', help='GROMACS .gro file for the restraints to be applied to the simulation')
    parser.add_argument('-t', '--topology', dest='topology', default='topology.top', help='GROMACS topology file in each replica')
    parser.add_argument('-x', '--replex', dest='replex', default='50', help='Interval for replica exchange')
    parser.add_argument('-E', '--nex', dest='nex', default=True, help='True will set the number of random exchanges to (number of replica)^3, False will perform nearest neighbor exchange')

    args = parser.parse_args()
    if os.path.isfile('EQ.gro'):
        initial_structure = 'EQ.gro'
    elif os.path.isfile('ANNEAL.gro'):
        initial_structure = 'ANNEAL.gro'
    else:
        initial_structure = 'pre_EQ.gro'

    run_REMD_gromacs_multi(int(args.num_replica), output_string=args.output_string, mdp=args.mdp, initial_structure=initial_structure,
                     restraints=args.restraints, topology=args.topology, replex=int(args.replex), nex=bool(args.nex))

