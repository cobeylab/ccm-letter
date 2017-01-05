#!/bin/sh
#SBATCH --array=1:192

# cd into the nth line of jobs.txt, where n is $SLURM_ARRAY_TASK_ID
cd `head -n $SLURM_ARRAY_TASK_ID jobs.txt | tail -n 1`
../../run_job.R
