#! /bin/bash

#SBATCH -N NODES
#SBATCH --ntasks-per-node 28
#SBATCH --time=48:00:00
#SBATCH --output=logfile
#SBATCH --partition=RM

module load mpi/intel_mpi
module load gromacs/2018_cpu

(for dir in {DIRS}; do cd $dir; gmx grompp -f production.mdp -c EQ.gro -r restraint.gro -p topology.top -o PROD.tpr  -maxwarn 10; cd ..; done)

mpirun -np NNPP  gmx_mpi mdrun -v -multidir {DIRS} -replex 50 -nex NNEEXX -deffnm PROD



