#!/bin/bash
COMP="/opt/intel/composerxe/bin/compilervars.sh"
ENV="/mydata/opt/matlab/matlab_env.sh" 

export LD_PRELOAD=/opt/intel/composerxe/mkl/lib/intel64/libmkl_rt.so:/opt/intel/composerxe/lib/intel64/libiomp5.so
export OMP_NUM_THREADS=8

if [ -f $COMP ];
then
	 source $COMP intel64
fi

if [ -f $ENV ];
then
	 source $ENV 
fi
echo "Lanzando experimento:" $EXP 

# ABSOLUTE_PATH should look like this --> /export/usuarios01/idUser/.../FOLDER/ 
# Where FOLDER must contain your SARIA.sh script

# "LICENSE"
#/mydata/opt/matlab/R2012b/bin/matlab -c /ABSOLUTE_PATH_TO_THE_LICENSE/license.dat -nodesktop < /ABSOLUTE_PATH_TO_YOUR_SCRIPTS_FOLDER/$EXP
# NO LICENSE
/mydata/opt/matlab/matlab79/bin/matlab -nojvm -nodesktop < /export/usuarios01/franrruiz87/simMIMO/$EXP
