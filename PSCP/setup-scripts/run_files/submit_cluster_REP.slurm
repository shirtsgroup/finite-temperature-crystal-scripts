#! /bin/bash

#SBATCH -N NODES
#SBATCH --ntasks-per-node 28
#SBATCH --time=48:00:00
#SBATCH --output=logfile
#SBATCH --partition=RM

Run_REMD.py -n NNPP -r restraint.gro
find . -iname \#* -type f -delete

