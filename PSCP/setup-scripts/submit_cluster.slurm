#! /bin/bash
#SBATCH --time=48:00:00
#SBATCH --output=logfile
#SBATCH --ntasks=aaaa
#SBATCH --partition=RM-shared

Run_simulation -a aaaa -r rrrr -i iiii -u uuuu -E EEEE -P PPPP -A AAAA

