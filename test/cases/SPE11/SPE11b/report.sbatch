#!/bin/bash

# execute with : sbatch report.sbatch conda_env job_id
# it creates the paraview, extract and report directories into output-SPE11_b
# conda_env is ComPASS-conda_env and must contain Paraview 5.12

#SBATCH -J pv_extract      # Job name
#SBATCH -p cpu             # Partition to submit to
#SBATCH --nodes=1          # Number of nodes
#SBATCH --ntasks-per-node=4         # maybe everything is sequential
#SBATCH -t 48:00:00        # Walltime
# Redirect output and error files to the result directory
#SBATCH -o ./Results/SPE11_b/report_%j.out  # STDOUT
#SBATCH -e ./Results/SPE11_b/report_%j.err  # STDERR
#SBATCH -d afterany:497934   # jobid of the *postprocess* job of scompass execution


# Stop the script if any command fails
set -e

# Verify if JOB_ID is provided as an argument
if [ -z "$2" ]; then
    echo "Error: No job ID provided. Usage: sbatch extract.sbatch <JOB_ID>"
    exit 1
fi

conda_evt=$1  # must contain Paraview 5.12
JOB_ID=$2

SCRATCH=/scratch/$USER/$SLURM_JOB_ID
SCRATCH_RESULT=$SCRATCH/output-SPE11_b
RESULT=$SLURM_SUBMIT_DIR/Results/SPE11_b/$JOB_ID/output-SPE11_b

echo "make report for JOB_ID=$JOB_ID"

# create the $SCRATCH directory
srun -m arbitrary -w $SLURM_NODELIST mkdir -p $SCRATCH
# copy from home to local
srun -m arbitrary -w $SLURM_NODELIST cp -r $RESULT $SCRATCH
srun -m arbitrary -w $SLURM_NODELIST cp pv_extract.py $SCRATCH
srun -m arbitrary -w $SLURM_NODELIST cp make_report.py $SCRATCH
echo "copies to $SCRATCH are done"
# go into the node local directory
cd $SCRATCH

conda run -n ComPASS-$conda_evt python -m ComPASS.postprocess -s -C -t "second" $SCRATCH_RESULT
# copy the paraview directory to the results directory
srun -m arbitrary -w $SLURM_NODELIST rsync -av $SCRATCH_RESULT/paraview $RESULT
echo "states postprocess is done"
# Run the pv_extract.py script using pvbatch
conda run -n ComPASS-$conda_evt mpirun -n $SLURM_NTASKS pvbatch pv_extract.py "$SCRATCH_RESULT/paraview/states.pvd" $SCRATCH_RESULT/cells_rkt.vtu -o "$SCRATCH_RESULT/extracts"
srun -m arbitrary -w $SLURM_NODELIST rsync -av $SCRATCH_RESULT/extracts $RESULT
echo "paraview extract is done"

# Make report scripts
conda run -n ComPASS-$conda_evt python make_report.py 'b' --src $SCRATCH_RESULT/extracts --out $SCRATCH_RESULT/report
srun -m arbitrary -w $SLURM_NODELIST rsync -av $SCRATCH_RESULT/report $RESULT
echo "make report is done"

# copy back the output directory (if other files have been modified than the paraview,
# extracts and report directories)
srun -m arbitrary -w $SLURM_NODELIST rsync -av $SCRATCH_RESULT/* $RESULT
echo "copies back (if any) are done"
# remove the scratch local directory on all nodes
srun -m arbitrary -w $SLURM_NODELIST rm -rf $SCRATCH
echo "cleaning scratch is done"
