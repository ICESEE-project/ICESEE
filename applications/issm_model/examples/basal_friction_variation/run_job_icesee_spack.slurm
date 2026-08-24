#!/bin/sh
#SBATCH -t 50:00:00
#SBATCH -J ISMIP
#SBATCH -N 2
#SBATCH -n 48
#SBATCH --ntasks-per-node=24
#SBATCH --partition=cpu-large
#SBATCH --mem=64G
#SBATCH -p inferno
#SBATCH -A gts-arobel3-atlas
#SBATCH -o steadystate-%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=bankyanjo@gmail.com

cd $SLURM_SUBMIT_DIR

# --- User parameters ---
NENS=80
MODEL_NPROCS=1
NP=16   # match -n above

echo "SLURM JOB INFO"
echo "Nodes: $SLURM_JOB_NUM_NODES"
echo "Tasks: $SLURM_NTASKS"

# --- Call run script ---
./run_icesee_spack.sh $NENS $MODEL_NPROCS $NP
