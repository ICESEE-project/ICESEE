#!/bin/sh

# Define directories
PERSISTENT_DIR="/scratch/icepack_cache"
APPTAINER_CACHEDIR="/scratch/apptainer_cache"
APPTAINER_TMPDIR="/scratch/apptainer_tmp"

# Ensure directories exist or clean them
for dir in "$APPTAINER_TMPDIR" "$APPTAINER_CACHEDIR" "$PERSISTENT_DIR"; do
    if [ ! -d "$dir" ]; then
        mkdir -p "$dir"
    else
        rm -rf "$dir"/*
    fi

    # Create persistent subdirectories only for PERSISTENT_DIR
    if [ "$dir" = "$PERSISTENT_DIR" ]; then
        mkdir -p "$PERSISTENT_DIR/pyop2" "$PERSISTENT_DIR/tsfc" "$PERSISTENT_DIR/xdg"
    fi
done

# Export environment variables
export APPTAINER_TMPDIR
export APPTAINER_CACHEDIR
export PYOP2_CACHE_DIR="$PERSISTENT_DIR/pyop2"
export FIREDRAKE_TSFC_KERNEL_CACHE_DIR="$PERSISTENT_DIR/tsfc"
export XDG_CACHE_HOME="$PERSISTENT_DIR/xdg"
# export PETSC_DIR=/opt/firedrake/petsc
# export PETSC_ARCH=complex
# export SLEPC_DIR=/opt/firedrake/slepc
# export MPICH_DIR=/opt/firedrake/petsc/packages/bin
# export HDF5_DIR=/opt/firedrake/petsc/packages
export HDF5_MPI=ON
