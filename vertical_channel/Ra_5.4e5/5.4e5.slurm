#!/bin/bash

#SBATCH -p normal
#SBATCH -J VC_5.4e5	      # Job name
#SBATCH -o job.%j.out         # Name of stdout output file (%j expands to jobId)
#SBATCH -e job.%j.err	      # Name of error file (%j expands to jobID)
#SBATCH -N 2                  # Total number of nodes requested
#SBATCH -n 136            # Total number of mpi tasks requested (1 node : 24 CPU)
#SBATCH -t 48:00:00           # Run time (hh:mm:ss), 12  hours
#SBATCH --mail-type=ALL       # mail alert
#SBATCH --mail-user=yunjaeh@stanford.edu 


echo "Running [your app] with: $cmd on $SLURM_JOB_NODELIST in directory "`pwd`
date

# Working directory
# /scratch/07279/tg865790/CharLES_scalar/vertical_channel/Ra_5.4e5
cd $SCRATCH
cd ./CharLES_scalar/vertical_channel/Ra_5.4e5

mpiexec -n 136 ../charles.exe -i charles_coarse.in > log.charles_coarse
mpiexec -n 136 ../charles.exe -i charles.in > log.charles

# mpiexec -n 68 $HOME/nextgen/bin/surfer.exe -i 01_surfer.in > log.surfer
# mpiexec -n 680 $HOME/nextgen/bin/charles.exe -i 05_charles.in > log.charles5
# mpiexec -n 68 $HOME/nextgen/bin/charles.exe -i 06_charles.in > log.charles6


