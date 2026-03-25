#!/bin/sh
#SBATCH -t230:00:00
#SBATCH -JISMIP
#SBATCH -N 5
#SBATCH -n 100
#SBATCH --ntasks-per-node=22
#SBATCH --partition=cpu-large         # Partition name
#SBATCH --mem=40G                     # Total memory per node
#SBATCH -pinferno
#SBATCH -A gts-arobel3-atlas
#SBATCH -osteadystate-%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=bkyanjo3@gatech.edu

cd $SLURM_SUBMIT_DIR
module purge
module load matlab/r2024b
module load gcc/13
#module load openmpi
#module load hdf5/1.12.3-ompi

# Activate the virtual environment
source ~/icesee-env/bin/activate

export ISSM_DIR=/storage/home/hcoda1/8/bkyanjo3/r-arobel3-0/ISSM
export LD_LIBRARY_PATH=/storage/home/hcoda1/8/bkyanjo3/r-arobel3-0/openmpi/lib:$LD_LIBRARY_PATH
export PATH=/storage/home/hcoda1/8/bkyanjo3/r-arobel3-0/openmpi/bin:$PATH

alias matlab='matlab -r "addpath ${ISSM_DIR}/bin ${ISSM_DIR}/lib"'

# ISSM environment
#source $ISSM_DIR/etc/environment.sh
#export LD_PRELOAD=/usr/local/pace-apps/spack/packages/linux-rhel9-x86_64_v3/gcc-11.3.1/gcc-12.3.0-ukkkutsxfl5kpnnaxflpkq2jtliwthfz/lib64/libstdc++.so.6

mpirun -np 24 python run_da_issm.py --Nens=40 --model_nprocs=4 

