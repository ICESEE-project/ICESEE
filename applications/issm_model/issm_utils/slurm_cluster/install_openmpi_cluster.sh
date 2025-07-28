#!/bin/bash

# Universal build script for OpenMPI
# Usage: ./build_openmpi.sh [prefix_dir]

# Default values
PREFIX_DIR=${1:-/storage/home/hcoda1/8/bkyanjo3/r-arobel3-0/openmpi}
SLURM_DIR=/opt/slurm/current
PMIX_DIR=/opt/pmix/5.0.1
MAKE_JOBS=16

# Check for required commands
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# for cmd in gcc g++ gfortran make; do
#     if ! command_exists "$cmd"; then
#         echo "Error: $cmd is not installed or not in PATH"
#         exit 1
#     fi
# done
module load gcc/13

# Check for required directories
for dir in "$SLURM_DIR" "$PMIX_DIR"; do
    if [ ! -d "$dir" ]; then
        echo "Error: Directory $dir does not exist"
        exit 1
    fi
done

# Ensure prefix directory exists
mkdir -p "$PREFIX_DIR" || {
    echo "Error: Cannot create prefix directory $PREFIX_DIR"
    exit 1
}

# Run configure with specified options
echo "Running configure..."
./configure \
    --prefix="$PREFIX_DIR" \
    --with-libevent \
    --with-hwloc \
    --with-ucx \
    --with-slurm="$SLURM_DIR" \
    --with-pmix="$PMIX_DIR" \
    --enable-mpi1-compatibility \
    CC=gcc \
    CXX=g++ \
    FC=gfortran || {
    echo "Error: Configure failed"
    exit 1
}

# Build with parallel jobs
echo "Building with $MAKE_JOBS parallel jobs..."
make -j"$MAKE_JOBS" || {
    echo "Error: Make failed"
    exit 1
}

# Install
echo "Installing..."
make install || {
    echo "Error: Make install failed"
    exit 1
}

echo "OpenMPI successfully built and installed to $PREFIX_DIR"