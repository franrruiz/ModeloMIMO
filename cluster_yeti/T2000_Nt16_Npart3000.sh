#!/bin/sh
#T2000_Nt16_Npart3000.sh
#Torque script to run Matlab program

#Torque directives
#PBS -N T2000_Nt16_Npart3000
#PBS -W group_list=yetidsi
#PBS -l nodes=1:ppn=16,walltime=72:00:00,mem=16gb
#PBS -M fr2392@columbia.edu
#PBS -m abe
#PBS -V

#set output and error directories (SSCC example here)
#PBS -o localhost:/vega/dsi/users/fr2392/simMIMO/
#PBS -e localhost:/vega/dsi/users/fr2392/simMIMO/

cd /vega/dsi/users/fr2392/simMIMO
export OMP_NUM_THREADS=16

T=2000
Nt=16
Nr=12
M=2
L=5
Niter=30000
lHead=0
onOffModel=0
Nparticles=3000
blockNtSize=6
flagParallel=1
itCluster=1
simId=3

#Command to execute Matlab code
matlab -nosplash -nodisplay -nodesktop -r "simClusterWise_yeti($T,$Nt,$Nr,$M,$L,$Niter,$lHead,$onOffModel,$Nparticles,$blockNtSize,$flagParallel,$itCluster,$simId)" > sal_T2000_Nt16_Npart3000.txt

#Command below is to execute Matlab code for Job Array (Example 4) so that each part writes own output
#matlab -nosplash -nodisplay -nodesktop -r "simPoissGLM($LAMBDA)" > matoutfile.$PBS_ARRAYID

#End of script