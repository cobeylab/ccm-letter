This code is organized to run in parallel on a SLURM cluster or other batch system as many
independent jobs, one for each combination of model assumptions/causal direction, but can
also be run on a single machine.

These instructions assume Unix (e.g. Mac OS X or Linux) with R installed so that the
Rscript executable is in the current search path.


INSTRUCTIONS
------------
1. (Local and Cluster)

Generate individual job directories in jobs/us=?-rz=?-usf=?-use=?-ulf=?-ev=?-fic=?

Guide to directory names:

us: use splines (0 or 1)
rz: remove zeros from flu *and* corresponding (non-zero) data points in environmental
    variable (if rz=0, zeros are still removed from predicted variable) (0 or 1)
usf: use surrogates for flu (0 or 1)
use: use surrogates for environmental variable (0 or 1)
ulf: log-transform flu (0 or 1)
ev: first letter of environmental variable
    (A)bsolute humidity, (R)elative humidity, (T)emperature, or (P)recipitation
fic: flu is cause (0 or 1)

cd <this-directory>
./create_jobs.R


2. (Cluster)

Modify run_all.sbatch for your SLURM cluster, or create an analogous file
for your cluster system.


2. (Local)

Modify N_CORES in run_job.R to reflect the number of cores you want jobs to use on the
local machine.


3. (Cluster)

Submit the cluster jobs, e.g.:

sbatch run_all.sbatch


3. (Local)

Either run individual models manually, via, e.g.:

cd jobs/us=0-rz=0-usf=0-use=1-ulf=0-ev=A-fic=0
../../run_job.R
cd ../..

or run them all sequentially via:

./run_all.sh

Each job produces an output file output_db.sqlite


4. (Local and Cluster)

Combine output files from all jobs into a single SQLite database:

./gather.R

which produces the file output_db_all.sqlite
