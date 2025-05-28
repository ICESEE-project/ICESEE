#!/bin/sh
#SBATCH -t30:00:00
#SBATCH -JISMIP
#SBATCH -N 3
#SBATCH -n 64
#SBATCH --ntasks-per-node=22
#SBATCH --mem-per-cpu=12G
#SBATCH -pinferno
#SBATCH -A gts-arobel3-atlas
#SBATCH -oTranisentSteadystate-%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=bkyanjo3@gatech.edu

cd $SLURM_SUBMIT_DIR
module purge
module load matlab/r2023b
module load gcc/12.3.0
module load openmpi

# ISSM environment
source $ISSM_DIR/etc/environment.sh
export LD_PRELOAD=/usr/local/pace-apps/spack/packages/linux-rhel9-x86_64_v3/gcc-11.3.1/gcc-12.3.0-ukkkutsxfl5kpnnaxflpkq2jtliwthfz/lib64/libstdc++.so.6

# activate ICESEE environment
source ~/icesee-env/bin/activate

srun -n 4 --mpi=pmi2 python run_da_issm.py --Nens=4 --model_nprocs=4 --verbose

# matlab -nodesktop -nosplash -r "run('issm_env'); run('runme')"
