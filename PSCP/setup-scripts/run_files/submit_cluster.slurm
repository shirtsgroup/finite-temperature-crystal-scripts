#! /bin/bash

#SBATCH --time=12:00:00
#SBATCH -N 1
#SBATCH --output=logfile
#SBATCH --ntasks=aaaa
#SBATCH --partition=RM-shared

Run_simulation.py -c aaaa -r rrrr -i iiii -u uuuu -e EEEE -p PPPP -a AAAA

