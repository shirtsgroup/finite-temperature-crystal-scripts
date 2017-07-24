#!/bin/bash

mkdir "icfile"

isocompress1="icfile/isocompressp1.txt"
press1="icfile/pressurep1.txt"
isocompress2="icfile/isocompressp2.txt"
press2="icfile/pressurep2.txt"
isocompress3="icfile/isocompressp3.txt"
press3="icfile/pressurep3.txt"



for direct in $(ls -l | egrep '^d' | awk '{print $9}'); do
	

	

	if [ "$(echo $direct | grep "p1")" == $direct ]; then
		( less $direct/benzene_production.mdp | grep "ref_p" | awk '{print $3}' )  >> $press1
		echo -e "Temperature \n Volume \n 0" | gmx energy -fluct_props -nmol 72 -f $direct/benzene_PROD.edr | grep "Kappa" | awk '{print $5}' >> $isocompress1

	elif [ "$(echo $direct | grep "p2")" == $direct ]; then
		( less $direct/benzene_production.mdp | grep "ref_p" | awk '{print $3}' )  >> $press2
		echo -e "Temperature \n Volume \n 0" | gmx energy -fluct_props -nmol 72 -f $direct/benzene_PROD.edr | grep "Kappa" | awk '{print $5}' >> $isocompress2
	else
		( less $direct/benzene_production.mdp | grep "ref_p" | awk '{print $3}' )  >> $press3
		echo -e "Temperature \n Volume \n 0" | gmx energy -fluct_props -nmol 72 -f $direct/benzene_PROD.edr | grep "Kappa" | awk '{print $5}' >> $isocompress3
	fi
done




