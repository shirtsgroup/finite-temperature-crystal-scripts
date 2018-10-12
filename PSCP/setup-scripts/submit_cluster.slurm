#! /bin/bash
#SBATCH --time=48:00:00
#SBATCH --output=logfile
#SBATCH --ntasks=aaaa
#SBATCH --partition=RM-shared

Run_simulation -a aaaa -r rrrr -E eeee -i iiii -u uuuu

