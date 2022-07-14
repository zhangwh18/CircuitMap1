# For pure serial jobs, @YNT=1
#!/bin/sh

#$ -N weihua

#$ -q single.q                            
#### single : For NOT high memory jobs (below 14GB memory per node). It uses just a single node.        ####
####          But your computation node may share with other user's jobs.                               ####    
####        

#$ -cwd  
#$ -pe openmpi 1        #### Instead of ' 12 ', input core number which you want to use, but there is a limitation  ####
                        #### in both matlab and mathematica : Matlab : 12, Mathematica : 8                          ####

# Run job through bash shell

#$ -S /bin/bash


set -e

source /etc/profile.d/modules.sh

MODULE_TO_LOAD=(python/gcc-4.8.5/3.5.6)
for m in ${MODULE_TO_LOAD[*]}; do
   module load $m
done

export OMP_NUM_THREADS=1   #### Instead of ' 12 ', input core number which you want to use, but there is a limitation  ####
                            #### in both matlab and mathematica : Matlab : 12, Mathematica : 8                          ####


for i in $(seq 1 28)
do

      gcc /home/weihua/CircuitMap_SRN_LO_initial_M1.c -o /home/weihua/circuit.out -lm  &
      /home/weihua/circuit.out i 

done

wait

      
####  Instead of 'command', input your corresponding command, such as matlab, math, and python  ####
####  : Input matlab (math) {python2.7, python3.5} for using Matlab (Mathematica) {Python2.7, python3.5}, respectively ####
    
####  Instead of '/scratch/hscheon/source_file_name', input your source file name including the path ####

####  Instead of '/scratch/hscheon/output_file_name.out', input your output file name including the path ####
