### Job script for OpenMPI_Intel ###

#
#! /bin/bash
# The walltime limit is not set, yet! E.g., you may run your code for infinite time by commenting out the line below
#$ -l h_rt=2880:00:00

##--- queue request
##--- debug.q, single.q, short_24.q, long_24.q, short_28.q, long_28.q,
#$ -q short_24.q
#### ,where 'mpi' is the name of a queue. Please read the fllows :                                      ####
#### mpi : For only mpi and for using nodes exclusively. But you can use only core numbers which are    ####
####       multiples of 24.  Max. Core number : 576 = 24 times 24, that is you can use up to 24 nodes   ####

#$ -l excl=true

#$ -v OMP_NUM_THREADS=1     #### If your code is not hybrid, you should always take 1.                  ####
                            #### If your code is hybrid, instead of '1', input your thread number       ####

#$ -pe orte 600                       #### Instead of 'number', input total core number which you want to use ####
#$ -N your_simulation_name
#$ -M your_email
#$ -m abe
# E.g. #$ -o /scratch/blabla/test.out
#$ -o /scratch/weihua/test_out         #### Input your output file name including the path    #####
#$ -e /scratch/weihua/test_err         #### Input error file name including the path          #####
#$ -S /bin/bash

echo "Preparing:"

cd /scratch/weihua                     #### Go to your working place in /scratch/yourname     #####

echo "Checking:"
pwd
hostname
date

echo "Environment:"

export OMP_NUM_THREADS=1        #### If your code is not hybrid, you should always take 1.                  ####
                                #### If your code is hybrid, instead of '1', input your thread number       ####
         #### Actually, this number should be the same as the number of "OMP_NUM_THREADS" of the above line ####

export PATH=/opt/openmpi.intel17/1.10.3/bin:$PATH
export LD_LIBRARY_PATH=/opt/openmpi.intel17/1.10.3/lib:/opt/intel/compilers_and_libraries_2017.1.132/linux/compiler/lib/intel64:$LD_LIBRARY_PATH


cat $PE_HOSTFILE | awk '{print $1, "slots="$2}' > hostfile

echo "Starting:"

# MPI run
mpirun -hostfile hostfile -n 600 /scratch/weihua/CircuitMap_SRN_LO.exe
####  If your code is not hybrid, instead of 'number', input your total core number which you want to use  ####
####  If your code is hybrid, instead of 'number', input "total core number/number of threads"             ####
####             Furthermore, you need to add the following option : "--bind-to core --map-by node:PE=number1" ####
####             Instead of 'number1', input your thread number                                                ####


#### Instead of '/scratch/hscheon/executable file name', input your executable file name including the path ####


echo "Stopping:"
rm -f hostfile
date

echo "Done."


