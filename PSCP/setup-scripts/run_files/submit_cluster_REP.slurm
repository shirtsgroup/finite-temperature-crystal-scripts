#! /bin/bash

#SBATCH -N NODES
#SBATCH --ntasks-per-node 28
#SBATCH --time=48:00:00
#SBATCH --output=logfile
#SBATCH --partition=RM

module load mpi/intel_mpi
module load gromacs/2018_cpu

Run_REMD.py -n NNPP -r restraint.gro -i EQ.gro


