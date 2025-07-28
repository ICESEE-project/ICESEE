#!/bin/bash
# Script for installing ISSM on a Linux-based HPC cluster with Slurm (e.g., PACE)
# This script assumes a module system and dependencies like GCC, MPI, and MATLAB
# Adjust module names, versions, and paths as needed for your cluster environment

# Usage: Run this script on the cluster's compute node 
# author: brian KYNAJO
# date: 2025-07-09

# Exit on error
set -e

# Define GCC library directory for PACE
export GCC_LIB_DIR="/usr/local/pace-apps/spack/packages/linux-rhel9-x86_64_v3/gcc-11.3.1/gcc-12.3.0-ukkkutsxfl5kpnnaxflpkq2jtliwthfz/lib64"

# Define installation directory, for pace -- we install in $HOME/p-arobel3-0/ISSM
export ISSM_DIR="$HOME/p-arobel3-0/ISSM"

# check if ISSM_DIR already exists and if so remove it
if [ -d "$ISSM_DIR" ]; then
    echo "Directory $ISSM_DIR already exists. Removing it..."
    rm -rf "$ISSM_DIR"
fi

# Create ISSM installation directory
mkdir -p "$ISSM_DIR"

# Clone ISSM repository
git clone https://github.com/ISSMteam/ISSM.git "$ISSM_DIR"
cd "$ISSM_DIR"

# Set up ISSM environment
source "$ISSM_DIR/etc/environment.sh"

# Download Slurm cluster configuration for ISSM
wget -O src/m/classes/clusters/generic.m \
    https://raw.githubusercontent.com/ICESEE-project/ICESEE/refs/heads/main/applications/issm_model/issm_utils/slurm_cluster/generic.m

# Load required modules (adjust module names/versions to match your cluster)
module load gcc/12
module load mvapich2
module load matlab/r2024b

# Install PETSc (v3.22.3)
cd "$ISSM_DIR/externalpackages/petsc"
mkdir -p src
git clone -b v3.22.3 https://gitlab.com/petsc/petsc.git src
cd src
./configure \
    --prefix="$ISSM_DIR/externalpackages/petsc/install" \
    --PETSC_DIR="$ISSM_DIR/externalpackages/petsc/src" \
    --with-mpi-dir="${MPI_ROOT}" \
    --with-debugging=0 \
    --with-valgrind=0 \
    --with-x=0 \
    --with-ssl=0 \
    --with-pic=1 \
    --download-fblaslapack=1 \
    --download-metis=1 \
    --download-mumps=1 \
    --download-parmetis=1 \
    --download-scalapack=1 \
    --download-zlib=1
make -j$(nproc) && make install


# Install external package: Triangle
cd "$ISSM_DIR/externalpackages/triangle"
./install-linux.sh

# Install external package: m1qn3
cd "$ISSM_DIR/externalpackages/m1qn3"
./install-linux.sh

# Refresh ISSM environment
source "$ISSM_DIR/etc/environment.sh"

# Compile ISSM
cd "$ISSM_DIR"
autoreconf -ivf

# Set compiler environment variables
export CC=mpicc
export CXX=mpicxx
export FC=mpifort

# Configure ISSM (adjust paths and flags as needed for your cluster)
./configure \
    --prefix="$ISSM_DIR" \
    --with-matlab-dir="${MATLABROOT}" \
    --with-fortran-lib="-L${GCC_LIB_DIR} -lgfortran" \
    --with-mpi-include="${MPI_ROOT}/include" \
    --with-mpi-libflags="-L${MPI_ROOT}/lib -lmpi -lmpichcxx -lmpifort" \
    --with-triangle-dir="$ISSM_DIR/externalpackages/triangle/install" \
    --with-petsc-dir="$ISSM_DIR/externalpackages/petsc/install" \
    --with-metis-dir="$ISSM_DIR/externalpackages/petsc/install" \
    --with-parmetis-dir="$ISSM_DIR/externalpackages/petsc/install" \
    --with-blas-lapack-dir="$ISSM_DIR/externalpackages/petsc/install" \
    --with-scalapack-dir="$ISSM_DIR/externalpackages/petsc/install" \
    --with-mumps-dir="$ISSM_DIR/externalpackages/petsc/install" \
    --with-m1qn3-dir="$ISSM_DIR/externalpackages/m1qn3/install" \
    --with-numthreads=20

# Compile and install ISSM
make -j$(nproc)
make install

echo "ISSM installation completed successfully in $ISSM_DIR"