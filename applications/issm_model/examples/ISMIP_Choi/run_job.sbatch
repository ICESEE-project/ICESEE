#!/bin/sh
#SBATCH -t230:00:00
#SBATCH -JISMIP
#SBATCH -N 4
#SBATCH -n 64
#SBATCH --ntasks-per-node=22
#SBATCH --mem-per-cpu=8G
#SBATCH -pinferno
#SBATCH -A gts-arobel3-atlas
#SBATCH -osteadystate-%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=bkyanjo3@gatech.edu

cd $SLURM_SUBMIT_DIR
module purge
module load matlab/r2024b
module load gcc/12.3.0
module load openmpi

# Activate the virtual environment
source ~/icesee-env/bin/activate

# ISSM environment
source $ISSM_DIR/etc/environment.sh
export LD_PRELOAD=/usr/local/pace-apps/spack/packages/linux-rhel9-x86_64_v3/gcc-11.3.1/gcc-12.3.0-ukkkutsxfl5kpnnaxflpkq2jtliwthfz/lib64/libstdc++.so.6

srun -n 8 python run_da_issm.py --Nens=40 --model_nprocs=4 --verbose

