#!/bin/bash
####
#$ -N GOES-FP
#$ -q ARROMA
##$ -l h=argon-lc-g21-20
#$ -pe smp 6
#$ -cwd
#$ -M yi-wang-4@uiowa.edu
#$ -m a
#$ -o ./$JOB_ID.out
#$ -e ./$JOB_ID.err
####

module load netcdf-fortran/4.4.4_parallel_studio-2017.4

# dkh debug unlimit core size and force core dump
ulimit
#ulimit -s unlimited
export KMP_STACKSIZE=512M

cd work_directory

time ./GeosFpDriver0.x < date_filename
