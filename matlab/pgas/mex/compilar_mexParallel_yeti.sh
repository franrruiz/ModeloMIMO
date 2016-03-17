#!/bin/sh
#compilar_mexParallel_yeti.sh
#Torque script to run Matlab program

#Torque directives
#PBS -N compilar_mexParallel_yeti
#PBS -W group_list=yetidsi
#PBS -l nodes=1,walltime=00:10:00,mem=50mb
#PBS -M fr2392@columbia.edu
#PBS -m abe
#PBS -V

#set output and error directories (SSCC example here)
#PBS -o localhost:/vega/dsi/users/fr2392/
#PBS -e localhost:/vega/dsi/users/fr2392/

#Command to execute Matlab code
matlab -nosplash -nodisplay -nodesktop -r compilar_mexParallel_yeti > matoutfile

#Command below is to execute Matlab code for Job Array (Example 4) so that each part writes own output
#matlab -nosplash -nodisplay -nodesktop -r "simPoissGLM($LAMBDA)" > matoutfile.$PBS_ARRAYID

#End of script