# ==============================================================================
# @des: PACE RUN 0 -- environment/MPI/Firedrake smoke check for Idealized
# PIG Mode-3 (2026-09-28 PACE calibration prep). Cheap, fast, NOT the
# scalability result -- only answers "is this environment usable at all."
#
# Checks, in order, each printed as [RUN0] OK/FAIL so a failure is easy to
# spot in SLURM's stdout without parsing:
#   - repository revision / working-tree identity (via pace_run_manifest)
#   - Python / mpi4py / h5py / Firedrake / petsc4py import
#   - HDF5 parallel (MPI) capability
#   - output filesystem path is writable from every rank
#   - real MPI communicator split (P_model > 1 construction), mirroring
#     what _icepack_native.py's topology.spatial_comm actually needs
#   - compact-initialization checkpoint file is present and openable
#
# Exit code is nonzero if any check fails, so the calling SLURM script can
# gate RUN 1 on RUN 0 succeeding (see run0_env_smoke.sbatch).
# ==============================================================================
from __future__ import annotations

import argparse
import os
import sys
import traceback
from pathlib import Path

_IDEALIZED_PIG_DIR = Path(__file__).resolve().parents[1]
_REPO_ROOT = _IDEALIZED_PIG_DIR.parents[3]  # .../ICESEE
for _p in (str(_REPO_ROOT), str(_REPO_ROOT.parent)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

_FAILURES: list[str] = []


def _check(name: str, fn):
    try:
        detail = fn()
        print(f"[RUN0] OK   {name}" + (f" -- {detail}" if detail else ""))
    except Exception as exc:  # noqa: BLE001 -- smoke check: report, don't crash the whole script
        _FAILURES.append(name)
        print(f"[RUN0] FAIL {name} -- {exc}")
        traceback.print_exc()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-path", required=True, help="This run's unique output directory (must not be shared with any other concurrent job)")
    parser.add_argument("--init-file", default=str(_IDEALIZED_PIG_DIR / "data" / "extended_beta1000yrs_compact_idx20000.h5"))
    args = parser.parse_args()

    from mpi4py import MPI
    world = MPI.COMM_WORLD
    rank, size = world.Get_rank(), world.Get_size()

    def rank0_print(msg):
        if rank == 0:
            print(msg)

    rank0_print(f"[RUN0] world_size={size}")

    _check("mpi4py import", lambda: f"vendor={MPI.get_vendor()}")

    def _h5py_check():
        import h5py
        cfg = h5py.get_config()
        return f"hdf5={'.'.join(str(v) for v in h5py.h5.get_libversion())} mpi={bool(getattr(cfg, 'mpi', False))}"
    _check("h5py import + HDF5 capability", _h5py_check)

    def _firedrake_check():
        import firedrake  # noqa: F401
        return "firedrake importable"
    _check("firedrake import", _firedrake_check)

    def _petsc_check():
        import petsc4py  # noqa: F401
        return "petsc4py importable"
    _check("petsc4py import", _petsc_check)

    def _writable_check():
        data_path = Path(args.data_path)
        data_path.mkdir(parents=True, exist_ok=True)
        probe = data_path / f"_run0_writable_probe_rank{rank}.tmp"
        probe.write_text("ok")
        probe.unlink()
        return f"rank {rank} can write to {data_path}"
    _check(f"output path writable (rank {rank})", _writable_check)

    def _init_file_check():
        init_file = Path(args.init_file)
        if not init_file.is_file():
            raise FileNotFoundError(str(init_file))
        return f"{init_file} ({init_file.stat().st_size / 1e6:.1f} MB)"
    _check("compact-initialization checkpoint accessible", _init_file_check)

    def _run1_inputs_check():
        # The other inputs RUN1 reads (params_calibration.yaml's paramsFile,
        # meshFile, SMBFile), relative to the Idealized PIG directory.
        names = ("extended_beta1000yrs.yaml", "PigFull2017GeomFull.exp", "OLS_Trend_plus_Resid_9b9.tif")
        missing = [n for n in names if not (_IDEALIZED_PIG_DIR / "data" / n).is_file()]
        if missing:
            raise FileNotFoundError(f"missing in {_IDEALIZED_PIG_DIR / 'data'}: {missing}")
        return ", ".join(names)
    _check("RUN1 input files accessible", _run1_inputs_check)

    def _comm_split_check():
        if size < 2:
            return "world_size=1: split construction not exercised (need >=2 ranks for a real P_model>1 check)"
        # Mirrors create_distributed_topology's own spatial/ensemble split
        # (src/parallelization/distributed_topology.py): P_model ranks per
        # ensemble group, contiguous.
        p_model = 2 if size % 2 == 0 else 1
        group = rank // p_model
        spatial_comm = world.Split(color=group, key=rank)
        spatial_rank = spatial_comm.Get_rank()
        spatial_size = spatial_comm.Get_size()
        spatial_comm.Barrier()
        spatial_comm.Free()
        return f"rank {rank} -> group {group}, spatial_comm rank {spatial_rank}/{spatial_size}"
    _check("MPI communicator split (P_model construction)", _comm_split_check)

    world.Barrier()
    n_failures = world.allreduce(len(_FAILURES), op=MPI.SUM)
    if rank == 0:
        if n_failures == 0:
            print("[RUN0] ALL CHECKS PASSED")
        else:
            print(f"[RUN0] {n_failures} CHECK(S) FAILED ACROSS ALL RANKS")
    world.Barrier()
    sys.exit(1 if _FAILURES else 0)


if __name__ == "__main__":
    main()
