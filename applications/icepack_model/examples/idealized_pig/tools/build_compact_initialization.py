#!/usr/bin/env python3
# ==============================================================================
# @des: Portable, read-only extraction utility that builds a compact Icepack
# production-initialization checkpoint from the full 1000-year spin-up
# history checkpoint (~39 GB).
#
# Root cause this addresses (Mode-3 scalability investigation): the
# production ``initFile`` stores a ``thickness``/``velocity``/``surface``
# TIME SERIES (idx 0..20000, one entry per spin-up step), but
# ``_icepack_model.py::initializeRun`` only ever reads exactly ONE index
# (``idx=20000``, the final/steady-state snapshot) for those three fields,
# plus five non-timestepped fields (``bed``, ``grounded``, ``floating``,
# ``fluidity``, ``extended_beta``) that were never time series to begin
# with. This script copies the mesh plus exactly those 8 fields into a new,
# standalone ``firedrake.CheckpointFile``.
#
# IMPORTANT, measured finding: ``CheckpointFile.save_function(f, idx=N)``
# ("timestepping mode") allocates on-disk storage proportional to N, not to
# how many indices are actually written -- writing ONLY idx=20000 into a
# fresh file still produced a 34.67 GB file (barely smaller than the 39.46
# GB source) in direct testing. The three previously-timestepped fields are
# therefore written here in NORMAL (non-timestepping) mode -- no ``idx=``
# -- which measured 22.2 MB total for mesh + one field in the same test.
# This means the compact file's on-disk layout for these three fields
# differs from the source (no idx dimension), so a small, explicitly
# opt-in read-side change is required to consume it -- see
# ``initializeRun``'s ``compact_initialization`` kwarg (default False,
# preserving the existing idx=20000 read path against the full file
# byte-for-byte unchanged).
#
# NEVER deletes or modifies the source file -- read-only open throughout.
# Portable: no hardcoded paths, no platform-specific shell calls; every
# location is a CLI argument, and the MPI communicator is whatever this
# process's own COMM_WORLD is (works identically under a bare
# ``python``/serial HDF5 build or a real ``mpirun`` launch with an
# MPI-enabled HDF5/PETSc build) -- exactly the existing
# ``initializeRun``/``getMeshFromCheckPoint`` convention of forwarding
# ``comm=`` rather than assuming COMM_WORLD or any specific filesystem.
#
# Usage:
#   python build_compact_initialization.py --source <full.h5> --dest <compact.h5> --idx 20000
#   mpirun -n P python build_compact_initialization.py --source ... --dest ... --idx 20000
# ==============================================================================
from __future__ import annotations

import argparse
import os
import time

# Timestepped fields initializeRun reads at exactly one idx, and the
# non-timestepped fields it reads unconditionally. Field names must match
# _icepack_model.py::initializeRun exactly -- these are not free choices.
TIMESTEPPED_FIELDS = ("velocity", "thickness", "surface")
STATIC_FIELDS = ("bed", "grounded", "floating", "fluidity", "extended_beta")


def build_compact_initialization(source, dest, idx, comm=None):
    """Copy the mesh and the 8 fields ``initializeRun`` needs from
    ``source`` (the full spin-up-history checkpoint) into a new, compact
    checkpoint at ``dest``, preserving field names and values. The three
    timestepped fields are read at history index ``idx`` and written
    without a history index (see the module header), so the destination
    holds only that one state. Returns a dict of measured timings/sizes.
    """
    import firedrake

    if comm is None:
        comm = firedrake.COMM_WORLD

    source = os.path.abspath(source)
    dest = os.path.abspath(dest)
    if os.path.abspath(source) == os.path.abspath(dest):
        raise ValueError("source and dest must be different files")

    timings = {}

    t0 = time.time()
    with firedrake.CheckpointFile(source, "r", comm=comm) as src:
        t_open_source = time.time() - t0
        timings["open_source_s"] = t_open_source

        t1 = time.time()
        mesh = src.load_mesh()
        timings["load_mesh_s"] = time.time() - t1

        fields = {}
        t2 = time.time()
        for name in TIMESTEPPED_FIELDS:
            fields[name] = src.load_function(mesh, name, idx=idx)
        for name in STATIC_FIELDS:
            fields[name] = src.load_function(mesh, name)
        timings["load_fields_s"] = time.time() - t2

    t3 = time.time()
    # Overwrite mode ("w"): dest is a NEW compact file this script owns,
    # never the source. If a stale dest from a prior attempt exists, start
    # clean rather than silently appending onto a possibly-inconsistent
    # file.
    #
    # All 8 fields are written in NORMAL (non-timestepping) mode -- no
    # idx= -- including the three that were timestepped in the source.
    # Measured directly: save_function(f, idx=20000) allocates on-disk
    # storage proportional to the index value even when only that one
    # index is ever written (34.67 GB for one field), while normal mode
    # for the same field measured 22.2 MB. The read side must match (see
    # initializeRun's compact_initialization kwarg).
    with firedrake.CheckpointFile(dest, "w", comm=comm) as out:
        out.save_mesh(mesh)
        for name in TIMESTEPPED_FIELDS:
            out.save_function(fields[name], name=name)
        for name in STATIC_FIELDS:
            out.save_function(fields[name], name=name)
    timings["write_dest_s"] = time.time() - t3

    timings["total_s"] = time.time() - t0

    if comm.rank == 0 and os.path.exists(dest):
        timings["dest_bytes"] = os.path.getsize(dest)
        timings["source_bytes"] = (
            os.path.getsize(source) if os.path.exists(source) else None
        )

    return timings


def _main():
    parser = argparse.ArgumentParser(
        description=(
            "Build a compact Icepack production-initialization checkpoint "
            "from a full spin-up-history checkpoint, without altering the "
            "source file."
        )
    )
    parser.add_argument("--source", required=True, help="path to the full checkpoint (.h5)")
    parser.add_argument("--dest", required=True, help="path to write the compact checkpoint (.h5)")
    parser.add_argument(
        "--idx", type=int, default=20000,
        help="timestep index to extract for velocity/thickness/surface (default: 20000, "
             "the production initializeRun value)",
    )
    args = parser.parse_args()

    timings = build_compact_initialization(args.source, args.dest, args.idx)

    from mpi4py import MPI
    if MPI.COMM_WORLD.rank == 0:
        print("[build_compact_initialization] timings/sizes:")
        for key, value in timings.items():
            print(f"  {key}: {value}")


if __name__ == "__main__":
    _main()
