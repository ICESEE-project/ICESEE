#!/bin/bash
#SBATCH -JIcepack                   # Job name
#SBATCH --account=gts-arobel3-atlas     # charge account
#SBATCH -N1 --tasks-per-node=2  # --gres=gpu:V100:1 --ntasks-per-node=2 
# # SBATCH --mem-per-gpu=12G 
#SBATCH -t 72:00:00          # Duration of the job 
#SBATCH -qinferno            # QOS Name inferno, embers
#SBATCH -oEnKf100-%j.out        # Combined output and error messages file
#SBATCH --mail-type=BEGIN,END,FAIL        # Mail preferences
#SBATCH --mail-user=bkyanjo3@gatech.edu  # E-mail address for notifications

# Run code
apptainer exec  icepack_working.sif python3 run_da_icepack.py
