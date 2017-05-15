#! /bin/bash
#SBATCH --partition=serial
#SBATCH --time=1-00:00:00
#SBATCH --output=logfile
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ecd4bd@virginia.edu

# submit a single SLURM job for a free energy calc
#module load gromacs/5.0.2-sse
export GRO_LOC="/h3/n1/shirtsgroup/gromacs_forceavg3/install/bin"
export SCRIPT_LOC="/nv/blue/nps5kd/scripts"
export THREADINFO="-nt 1 "

export OMP_NUM_THREADS=1

sleep 1

#Change the job status to 'SUBMITTED'
echo "SUBMITTED" > jobstatus.txt

#EQUILIBRATION
${GRO_LOC}/grompp_d -f LJEquil.mdp -c fccinitial.gro -r fccinitial.gro -p binaryLJ.top -o LJ_EM.tpr -maxwarn 10
${GRO_LOC}/mdrun_d ${THREADINFO} -v -deffnm LJ_EM

sleep 1

#EXTRACT THE FINAL FRAME FROM THE OUTPUT XTC FILE AND USE IT AS THE FINAL EQUILIBRATION FRAME
echo "Total-Energy" | ${GRO_LOC}/g_energy_d -f LJ_EM.edr -o time.xvg
finaltime=$(tail -1 time.xvg | awk '{print $1}')
rm LJ_EM.gro
echo 0 | ${GRO_LOC}/trjconv_d -f LJ_EM.xtc -s LJ_EM.tpr -o LJ_EM.gro -pbc whole -ndec 12 -dump $finaltime

sleep 1

#PRODUCTION
${GRO_LOC}/grompp_d -f LJprod.mdp -c LJ_EM.gro -r fccinitial.gro -p binaryLJ.top -o LJ_PROD.tpr -maxwarn 10
${GRO_LOC}/mdrun_d ${THREADINFO} -v -deffnm LJ_PROD

sleep 1

#EXTRACT THE FINAL FRAME FROM THE OUTPUT TRR FILE AND USE IT AS THE FINAL PRODUCTION FRAME
echo "Total-Energy" | ${GRO_LOC}/g_energy_d -f LJ_PROD.edr -o time.xvg
finaltime=$(tail -1 time.xvg | awk '{print $1}')
rm LJ_PROD.gro
echo 0 | ${GRO_LOC}/trjconv_d -f LJ_PROD.trr -s LJ_PROD.tpr -o LJ_PROD.gro -pbc whole -ndec 12 -dump $finaltime
echo 0 | ${GRO_LOC}/trjconv_d -f LJ_PROD.trr -s LJ_PROD.tpr -o LJ_PROD.pdb -pbc nojuump -ndec 8

sleep 1

#CALCULATE THE CONFIGURATION ENERGIES AND ENSEMBLE AVERAGES
${SCRIPT_LOC}/crunchjobenergy


sleep 1


#DELETE THE LOGFILE IF THE JOB FINISHED CORRECTLY
if [ -f benzene_PROD.log ] && [ "$(tail -1 benzene_PROD.log | awk '{print $1}')" == "Finished" ]; then
    rm logfile
fi



echo "FINISHED" > jobstatus.txt
# print end time
echo
echo "Job Ended at `date`"
echo "###################################################################"
