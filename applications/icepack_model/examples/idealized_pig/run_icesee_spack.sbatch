#!/bin/sh
#SBATCH -t 50:00:00
#SBATCH -J idealized_pig
#SBATCH -N 2
#SBATCH -n 48
#SBATCH --ntasks-per-node=24
#SBATCH --partition=cpu-large
#SBATCH --mem=128G
#SBATCH -p inferno
#SBATCH -A gts-arobel3-atlas
#SBATCH -o idealized_pig-%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jcruz89@gatech.edu

cd $SLURM_SUBMIT_DIR

# --- User parameters ---
NENS=32
MODEL_NPROCS=1
NP=8   # match -n above

echo "SLURM JOB INFO"
echo "Nodes: $SLURM_JOB_NUM_NODES"
echo "Tasks: $SLURM_NTASKS"

# --- Call run script ---
./run_icesee_spack.sh $NENS $MODEL_NPROCS $NP