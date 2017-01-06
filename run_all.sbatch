#!/bin/sh
#SBATCH --partition=sandyb
#SBATCH --array=1-192
#SBATCH --time=4:00:00
#SBATCH --output=slurm_out/output-%A_%a.txt

# cd into the nth line of jobs.txt, where n is $SLURM_ARRAY_TASK_ID
JOB_DIR=`head -n $SLURM_ARRAY_TASK_ID jobs.txt | tail -n 1`
echo $JOB_DIR
cd $JOB_DIR
../../run_job.R
