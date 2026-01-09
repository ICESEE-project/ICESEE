#!/bin/bash
#!/bin/bash
#SBATCH -JIcepack                   # Job name
#SBATCH --account=gts-gburdell3     # charge account
#SBATCH -N1 --ntasks-per-node=4     # Number of nodes and cores per node required
#SBATCH --mem-per-cpu=1G                  # Memory per core
#SBATCH -t15                              # Duration of the job (Ex: 15 mins)
#SBATCH -qinferno                         # QOS Name
#SBATCH -oReport-%j.out                   # Combined output and error messages file
#SBATCH --mail-type=BEGIN,END,FAIL        # Mail preferences
#SBATCH --mail-user=bkyanjo33@gatech.edu  # E-mail address for notifications
cd $SLURM_SUBMIT_DIR                      # Change to working directory

# Load module dependencies 
module load gcc
module load mvapich2

mpiexec -n 4 apptainer exec  icepack.sif python3 test.py