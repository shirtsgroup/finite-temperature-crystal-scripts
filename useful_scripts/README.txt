SCRIPTS DIRECTORY
=========================================

This directory contains the scripts necessary to create and manipulate jobs. The scripts fall in to these categories

JOB CREATION
=========================================
setup_molecule - 	This is the main job creation script that creates one MD jobs. This is called by every other job creation script

setup_Restraints -	This creates the jobs necessary to add harmonic restraints to a system according to the PSCP

setup_Interactions-	This creates the jobs necessary to turn off the interactions of a harmonically restrained system according to the PSCP

setup_Charges -		This creates the jobs to scan over a set of charges with an NPT simulation

setup_Interpolation-	This creates the jobs necessary to transform from one hamiltonian to another via linear interpolation

setup_Check -		This creates the jobs necessary to run the check-ensemble subroutine



JOB DIRECTORY MANIPULATION
=========================================
All scripts of this type have the flags -h "<include text>" -i "<ignore text>" and will act on only jobs in the current directory whose names containthe include text and do not contain the ignore text. 

listjobs -		List all jobs

movejobsitc -		Move all jobs to the scratch directory of Rivanna

movejobs<element>	Move all jobs to the local job directory of a desktop

submitjobs-		Submit all jobs to the Rivanna cluster

checkjobs-		Check the current status of all jobs. Either 'Running' 'Finished' 'Failed' or 'Unsubmitted'

killjobs-		Cancel the simulation of all jobs

resetjobs-		Reset all jobs to their pre-submitted status

resubmitjobs-		Resubmit all jobs that have failed. This will sequentially call reset and submit. This does not affect jobs that have the status of 'Finished' or 'Running'


POST-SIMULATION MANIPULATION
=========================================
convertjobtinker -	Reweights the current directory into the AMOEBA
potential

reweightjobgromacs -	Reweights the current directory into a potential
implemented within GROMACS

crunchallenergy - 	Extract information from a finished jobs including
pressure, volume, and potential energy.
