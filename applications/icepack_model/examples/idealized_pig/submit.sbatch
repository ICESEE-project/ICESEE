#!/bin/bash
#SBATCH -J Icepack
#SBATCH --account=gts-arobel3-atlas
#SBATCH -N 2
#SBATCH -n 40
#SBATCH --ntasks-per-node=20
#SBATCH --mem-per-cpu=12G
#SBATCH -t 10:00:00
#SBATCH -q inferno
#SBATCH -o EnKf-%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=bkyanjo3@gatech.edu

# Run code
srun -n 40 --mpi=pmi2 apptainer exec icepack_working.sif python3 run_da_icepack.py --Nens=40
