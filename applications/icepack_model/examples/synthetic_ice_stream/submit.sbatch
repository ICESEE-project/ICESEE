#!/bin/bash
#SBATCH -JIcepack                   # Job name
#SBATCH --account=gts-arobel3-atlas     # charge account
#SBATCH -N2 -n20  # --gres=gpu:V100:1 --ntasks-per-node=2 
#SBATCH --mem-per-cpu=64G 
#SBATCH -t 20:00:00          # Duration of the job 
#SBATCH -qinferno            # QOS Name inferno, embers
#SBATCH -oEnKf-%j.out        # Combined output and error messages file
#SBATCH --mail-type=BEGIN,END,FAIL        # Mail preferences
#SBATCH --mail-user=bkyanjo3@gatech.edu  # E-mail address for notifications

# Run code
srun --mpi=pmi2 apptainer exec icepack_working.sif python3 run_da_icepack.py --Nens=40
