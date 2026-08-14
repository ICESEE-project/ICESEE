# =============================================================================
# @author: Brian Kyanjo
# @date: 2025-09-07
# @description: - Class to handle parallel I/O operations for Ensemble Kalman Filter (EnKF) data.
#                 This class is designed to work with MPI for parallel processing and supports both
#                 serial and parallel file batch creation modes.
#               - It extends the EnKFIO_zarr class to provide additional functionality specific to
#                 parallel I/O operations with zarr format.
#               - EnkF analysis step utils have been added to this class including generation of synthetic
#                 observations, and observation operator.
#               - The analysis step for the EnKF and its mean have also been parallelized and added here.
# =============================================================================

import h5py
import numpy as np
from mpi4py import MPI
import copy
import os
import re
import glob
import gc
import zarr
import traceback
import sys
import shutil
# from numcodecs import blosc
import time
import functools
from concurrent.futures import ThreadPoolExecutor, as_completed

# blosc.use_threads = False

from typing import Callable, Optional, TypeVar, Any
from ICESEE.src.utils.tools import icesee_get_index
from ICESEE.src.utils.localization import iter_generated_observation_columns

T = TypeVar("T")


def normalize_hdf5_compression(value):
    """Return a valid h5py compression selector from configuration input.

    YAML 1.1 loaders do not consistently interpret the spelling ``None`` as
    a null value.  Several ICESEE application files historically use that
    spelling, so mode 2 may receive the literal string ``"None"``.  Passing
    it to h5py asks HDF5 for a filter named ``None`` and fails while creating
    the first streamed state shard.
    """
    if value is None or value is False:
        return None
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in {"", "none", "null", "false", "off"}:
            return None
        if normalized in {"gzip", "lzf", "szip"}:
            return normalized
    raise ValueError(
        "h5_file_compression must be None, 'gzip', 'lzf', or 'szip'; "
        f"received {value!r}."
    )


def hdf5_compression_options(compression, level):
    """Return filter options accepted by h5py for ICESEE's configuration."""
    if compression == "gzip":
        level = int(level)
        if not 0 <= level <= 9:
            raise ValueError("h5_file_compression_level must be in [0, 9].")
        return level
    # lzf accepts no options.  For szip, let h5py use its portable default;
    # an integer gzip-style level is not a valid szip option tuple.
    return None


def _partition_1d(n, size, rank):
    """Return a balanced half-open interval for ``rank``."""
    q, r = divmod(int(n), int(size))
    start = rank * q + min(rank, r)
    return start, start + q + (1 if rank < r else 0)


def resolve_row_chunk_size(nrows, nens, budget_mb=256.0, explicit_rows=0,
                           working_buffers=3):
    """Choose a row chunk that bounds the main dense working set.

    The estimate is intentionally conservative: it accounts for several
    ``rows x ensemble`` float64 buffers that BLAS and NumPy may hold at once.
    """
    nrows = int(nrows)
    if nrows <= 0:
        return 0
    explicit_rows = int(explicit_rows or 0)
    if explicit_rows > 0:
        return max(1, min(nrows, explicit_rows))
    bytes_per_row = max(1, int(nens)) * 8 * max(1, int(working_buffers)) + 8
    budget = max(1.0, float(budget_mb)) * 1024.0 * 1024.0
    return max(1, min(nrows, int(budget // bytes_per_row)))


def legacy_transform_from_gram(gram, rhs, energy_fraction=0.999):
    """Compute the legacy ICESEE ensemble transform from streamed products.

    ``gram`` is ``A.T @ A`` and ``rhs`` is ``A.T @ D`` for
    ``A = 2 * Yprime``.  This is algebraically equivalent to the previous
    observation-space SVD, but only stores ``Nens x Nens`` arrays.
    """
    gram = 0.5 * (np.asarray(gram) + np.asarray(gram).T)
    rhs = np.asarray(rhs)
    values, vectors = np.linalg.eigh(gram)
    order = np.argsort(values)[::-1]
    values = np.maximum(values[order], 0.0)
    vectors = vectors[:, order]
    total = float(values.sum())
    inverse = np.zeros_like(values)
    if total > 0.0:
        cumulative = np.cumsum(values)
        nkeep = int(np.searchsorted(cumulative, energy_fraction * total,
                                    side="left") + 1)
        tolerance = np.finfo(values.dtype).eps * max(gram.shape) * values[0]
        keep = (np.arange(values.size) < nkeep) & (values > tolerance)
        inverse[keep] = 1.0 / values[keep]
    gram_pinv = (vectors * inverse) @ vectors.T
    # HAprime.T = A.T / 2 in the legacy_prior_anomalies mode.
    return np.eye(gram.shape[0], dtype=gram.dtype) + 0.5 * gram_pinv @ rhs


def ensemble_transform_from_products(cross, gram, rhs,
                                     energy_fraction=0.999):
    """Compute the Evensen transform from streamed ensemble-space products.

    For ``A = Y' + Eta`` the three inputs are ``Y'.T @ A``, ``A.T @ A``
    and ``A.T @ D'``.  The identity

    ``(A A.T)^dagger = A (A.T A)^dagger (A.T A)^dagger A.T``

    makes this algebraically equivalent to the observation-space SVD while
    retaining only ``Nens x Nens`` matrices.  This is valid for both the
    standard stochastic and legacy analysis-factor modes.
    """
    cross = np.asarray(cross, dtype=np.float64)
    gram = np.asarray(gram, dtype=np.float64)
    rhs = np.asarray(rhs, dtype=np.float64)
    gram = 0.5 * (gram + gram.T)
    values, vectors = np.linalg.eigh(gram)
    order = np.argsort(values)[::-1]
    values = np.maximum(values[order], 0.0)
    vectors = vectors[:, order]
    total = float(values.sum())
    inverse = np.zeros_like(values)
    if total > 0.0:
        cumulative = np.cumsum(values)
        nkeep = int(np.searchsorted(
            cumulative, float(energy_fraction) * total, side="left"
        ) + 1)
        tolerance = np.finfo(values.dtype).eps * max(gram.shape) * values[0]
        keep = (np.arange(values.size) < nkeep) & (values > tolerance)
        inverse[keep] = 1.0 / values[keep]
    gram_pinv = (vectors * inverse) @ vectors.T
    return np.eye(gram.shape[0]) + cross @ gram_pinv @ gram_pinv @ rhs


def _read_rows_preserve_order(dataset, indices, column_slice=slice(None)):
    """Read arbitrary HDF5 rows, preserving duplicates and caller order."""
    indices = np.asarray(indices, dtype=np.int64)
    if indices.size == 0:
        shape = (0,) + tuple(dataset.shape[1:])
        return np.empty(shape, dtype=dataset.dtype)
    unique, inverse = np.unique(indices, return_inverse=True)
    return np.asarray(dataset[unique, column_slice])[inverse]


def _global_variable_index_map(nd, vec_inputs, var_nd=None):
    """Build state-vector variable indices without rank-local MPI slicing."""
    vec_inputs = list(vec_inputs or [])
    if not vec_inputs:
        return {}
    if isinstance(var_nd, dict):
        sizes = [int(var_nd[name]) for name in vec_inputs]
    elif isinstance(var_nd, (list, tuple, np.ndarray)) and len(var_nd) == len(vec_inputs):
        sizes = [int(value) for value in var_nd]
    else:
        q, r = divmod(int(nd), len(vec_inputs))
        sizes = [q + (1 if i < r else 0) for i in range(len(vec_inputs))]
    if sum(sizes) != int(nd):
        raise ValueError(
            f"State block sizes sum to {sum(sizes)}, but the state dimension is {nd}."
        )
    result = {}
    offset = 0
    for name, size in zip(vec_inputs, sizes):
        result[name] = np.arange(offset, offset + size, dtype=np.int64)
        offset += size
    return result


def _inversion_variable_info(nd, icesee_kwargs):
    """Resolve the full-state friction block used by hybrid inversion.

    Mode 1 temporarily removes this block from the EnKF state.  Mode 2 keeps
    the file-backed state shape fixed and obtains the same result by excluding
    its observations and restoring its forecast rows before inversion.
    """
    if not bool(icesee_kwargs.get("inversion_flag", False)):
        return None
    vec_inputs = list(
        icesee_kwargs.get("vec_inputs_old")
        or icesee_kwargs.get("vec_inputs")
        or []
    )
    if not vec_inputs:
        raise ValueError("Hybrid inversion requires vec_inputs_old or vec_inputs.")
    friction_idx = icesee_kwargs.get("friction_idx")
    if friction_idx is None:
        aliases = {
            "coefficient", "friction", "friction_coefficient", "fcoef",
            "frictioncoefficient",
        }
        matches = [
            index for index, name in enumerate(vec_inputs)
            if str(name).lower() in aliases
        ]
        if len(matches) != 1:
            raise ValueError(
                "Hybrid inversion could not uniquely identify the friction "
                "block; set friction_idx explicitly."
            )
        friction_idx = matches[0]
    friction_idx = int(friction_idx)
    if friction_idx < 0 or friction_idx >= len(vec_inputs):
        raise IndexError(
            f"friction_idx={friction_idx} is outside {len(vec_inputs)} blocks."
        )
    variable_map = _global_variable_index_map(
        nd, vec_inputs, icesee_kwargs.get("var_nd_old", icesee_kwargs.get("var_nd"))
    )
    name = vec_inputs[friction_idx]
    indices = variable_map[name]
    if not indices.size:
        raise ValueError(f"The inversion block {name!r} is empty.")
    return {
        "name": name,
        "index": friction_idx,
        "indices": indices,
        "start": int(indices[0]),
        "stop": int(indices[-1]) + 1,
    }


def _exclude_inversion_observations(obs_indices, icesee_kwargs, nd):
    """Return the observation-row mask matching mode-1 state reduction."""
    info = _inversion_variable_info(nd, icesee_kwargs)
    obs_indices = np.asarray(obs_indices, dtype=np.int64)
    if info is None:
        return np.ones(obs_indices.size, dtype=bool)
    return (obs_indices < info["start"]) | (obs_indices >= info["stop"])

def retry_on_failure(
    max_attempts: int = 5,
    delay: float = 1.0,
    mpi_comm: Optional[Any] = None
) -> Callable:
    """
    A decorator to retry a function or method up to max_attempts times with a delay between attempts.

    Args:
        max_attempts (int): Maximum number of retry attempts (default: 5).
        delay (float): Seconds to wait between retries (default: 1.0).
        mpi_comm: Optional MPI communicator object for distributed environments (default: None).

    Returns:
        Callable: The wrapped function with retry logic.
    """
    def decorator(func: Callable[..., T]) -> Callable[..., T]:
        @functools.wraps(func)
        def wrapper(*args, **icesee_kwargs) -> T:
            rank = mpi_comm.Get_rank() if mpi_comm is not None else "N/A"
            for attempt in range(max_attempts):
                try:
                    return func(*args, **icesee_kwargs)
                except IndexError as e:
                    if attempt < max_attempts - 1:
                        print(f"[Rank {rank}] Attempt {attempt + 1} failed with IndexError: {e}. Retrying in {delay}s...")
                        time.sleep(delay)
                    else:
                        print(f"[Rank {rank}] Final attempt failed with IndexError: {e}. Aborting.")
                        raise
                except Exception as e:
                    if attempt < max_attempts - 1:
                        print(f"[Rank {rank}] Attempt {attempt + 1} failed with {type(e).__name__}: {e}. Retrying in {delay}s...")
                        time.sleep(delay)
                    else:
                        print(f"[Rank {rank}] Final attempt failed with {type(e).__name__}: {e}. Aborting.")
                        raise
        return wrapper
    return decorator

class EnKF_fully_parallel_IO:
    def __init__(self, file_prefix, nd, nens, nt, subcomm, mpi_comm, icesee_kwargs, \
                 serial_file_creation=False, base_path="enkf_data", batch_size=50,\
                 h5_file_compression=None, h5_file_compression_level=4, h5_file_chunk_size=1000):
        try:
            self.nd = nd
            self.nens = nens
            # ``nt`` is the number of model advances.  Ensemble history also
            # contains the initial condition, so valid state shards are
            # 0, ..., nt (nt + 1 snapshots), matching execution mode 1.
            self.nt = int(nt) + 1
            self.icesee_kwargs = icesee_kwargs
            self.base_path = base_path
            self.file_prefix = file_prefix
            self.batch_size = max(1, int(batch_size))
            requested_history_mode = str(
                icesee_kwargs.get("ensemble_history_mode", "auto")
            ).strip().lower()
            if requested_history_mode not in {"auto", "full", "rolling"}:
                raise ValueError(
                    "ensemble_history_mode must be 'auto', 'full', or 'rolling'."
                )
            self.mpi_comm = mpi_comm
            # self.comm = subcomm if nens >= mpi_comm.Get_size() else mpi_comm
            # self.rank = self.comm.Get_rank() if nens >= mpi_comm.Get_size() else mpi_comm.Get_rank()
            # self.size = self.comm.Get_size() if nens >= mpi_comm.Get_size() else mpi_comm.Get_size()
            self.comm = mpi_comm
            self.rank = mpi_comm.Get_rank()
            self.size = mpi_comm.Get_size()

            self.subcomm = subcomm if subcomm is not None else MPI.COMM_SELF
            self.sub_rank = self.subcomm.Get_rank()
            self.sub_size = self.subcomm.Get_size()
            self.subcomm = subcomm
            self.serial_file_creation = serial_file_creation
            self.h5_file_compression = normalize_hdf5_compression(
                h5_file_compression
            )
            self.h5_file_compression_level = h5_file_compression_level
            self.h5_file_compression_options = hdf5_compression_options(
                self.h5_file_compression,
                h5_file_compression_level,
            )
            self.h5_file_chunk_size = h5_file_chunk_size
            self.storage_dtype = np.dtype(
                icesee_kwargs.get("ensemble_storage_dtype", "float64")
            )
            if self.storage_dtype not in (np.dtype("float32"), np.dtype("float64")):
                raise ValueError(
                    "ensemble_storage_dtype must be 'float32' or 'float64'."
                )
            self.reuse_existing = bool(
                icesee_kwargs.get("restart_enabled", True)
                and not icesee_kwargs.get("force_fresh_start", False)
            )

            # collective I/O threshold (32 is a reasonable default for many systems)
            self.collective_threshold = int(icesee_kwargs.get("collective_threshold", 16))
            self.use_collective_io = (self.mpi_comm.Get_size() >= self.collective_threshold)

            def partition_1d(n, size, rank):
                q, r = divmod(n, size)
                start = rank*q + min(rank, r)
                stop  = start + q + (1 if rank < r else 0)
                return start, stop   # [start, stop) half-open
            self.nd_start, self.nd_end = partition_1d(self.nd, self.size, self.rank)
            self.nd_local = self.nd_end - self.nd_start
            size_world = mpi_comm.Get_size()
            rank_world = mpi_comm.Get_rank()
            self.nd_start_world, self.nd_end_world = partition_1d(self.nd, size_world, rank_world)
            self.nd_local_world = self.nd_end_world - self.nd_start_world

            # Create directory.  Existing shards are retained for restart;
            # force_fresh_start is the explicit opt-in for deleting them.
            if mpi_comm.Get_rank() == 0:
                os.makedirs(base_path, exist_ok=True)
                if not self.reuse_existing:
                    patterns = [f"{self.base_path}/{self.file_prefix}_*.h5"]
                    for pattern in patterns:
                        for file_path in glob.glob(pattern):
                            try:
                                os.remove(file_path)
                            except OSError as e:
                                print(f"Error deleting {file_path}: {e}")
            self.mpi_comm.Barrier()

            bytes_per_step = self.nd * self.nens * self.storage_dtype.itemsize
            safety = float(icesee_kwargs.get("storage_safety_factor", 1.10))
            if self.rank == 0:
                free_bytes = shutil.disk_usage(self.base_path).free
                full_history_bytes = bytes_per_step * self.nt
                if requested_history_mode == "auto":
                    selected_history_mode = (
                        "full"
                        if full_history_bytes * safety <= free_bytes
                        else "rolling"
                    )
                    if selected_history_mode == "rolling":
                        print(
                            "[ICESEE][full_parallel] Full ensemble history requires "
                            f"about {full_history_bytes / 2**30:.2f} GiB; selecting "
                            "rolling history automatically."
                        )
                else:
                    selected_history_mode = requested_history_mode
            else:
                free_bytes = None
                selected_history_mode = None

            self.history_mode = self.mpi_comm.bcast(
                selected_history_mode, root=0
            )
            icesee_kwargs["ensemble_history_mode"] = self.history_mode

            if self.rank == 0:
                retained_steps = self.nt if self.history_mode == "full" else min(2, self.nt)
                history_bytes = bytes_per_step * retained_steps
                print(
                    "[ICESEE][full_parallel] File-backed ensemble: "
                    f"{bytes_per_step / 2**30:.2f} GiB per time step; "
                    f"about {history_bytes / 2**30:.2f} GiB retained in "
                    f"{self.history_mode!r} history mode ({self.storage_dtype.name})."
                )
                if history_bytes * safety > free_bytes:
                    message = (
                        "[ICESEE][full_parallel] Estimated ensemble history exceeds "
                        f"available storage ({free_bytes / 2**30:.2f} GiB free). "
                        "Use a larger data_path, ensemble_history_mode: rolling, "
                        "ensemble_storage_dtype: float32, or reduce the run size."
                    )
                    if icesee_kwargs.get("fail_on_insufficient_storage", False):
                        raise OSError(message)
                    print("WARNING: " + message)

            # Initialize file and dataset lists
            self.files = []
            self.datasets = []
            self.current_batch_start = -1

            # Timestep shards are opened lazily.  This is essential on restart:
            # a rolling-history run may legitimately no longer contain step 0.
        except Exception as e:
            print(f"Error occurred in __init__: {e}")
            tb_str = "".join(traceback.format_exception(*sys.exc_info()))
            print(f"Traceback details:\n{tb_str}")
            self.mpi_comm.Abort(1)

    def _state_file_name(self, timestep):
        """Return the canonical path for one ensemble-state shard."""
        return os.path.join(
            self.base_path,
            f"{self.file_prefix}_{int(timestep):04d}.h5",
        )

    def _analysis_file_name(self, timestep):
        """Return the unpublished, transactional analysis-shard path."""
        return self._state_file_name(timestep) + ".analysis.tmp"

    def _create_analysis_file(self, timestep):
        """Create and collectively open a temporary analysis shard."""
        path = self._analysis_file_name(timestep)
        self.mpi_comm.Barrier()
        if self.rank == 0:
            try:
                os.remove(path)
            except FileNotFoundError:
                pass
            with h5py.File(path, "w") as handle:
                handle.create_dataset(
                    "states",
                    shape=(self.nd, self.nens),
                    chunks=(max(1, min(self.h5_file_chunk_size, self.nd)), 1),
                    compression=self.h5_file_compression,
                    compression_opts=self.h5_file_compression_options,
                    dtype=self.storage_dtype,
                )
        self.mpi_comm.Barrier()
        handle = h5py.File(path, "r+", driver="mpio", comm=self.mpi_comm)
        return path, handle, handle["states"]

    @staticmethod
    def _requires_shared_finalization(icesee_kwargs):
        """Whether publishing requires the shared ensemble-level contract."""
        if str(icesee_kwargs.get("model_name", "")).lower() == "issm":
            return True
        if icesee_kwargs.get("inference_plugin_enabled", False):
            return True
        if icesee_kwargs.get("physics_bed_inference", False):
            return True
        if str(icesee_kwargs.get("bed_update_domain", "all")).lower() != "all":
            return True
        return bool(icesee_kwargs.get("bed_observation_anchor_enabled", False))

    @staticmethod
    def _requires_ensemble_coupled_finalization(icesee_kwargs):
        """Return whether post-analysis physics couples ensemble columns.

        The ordinary ISSM geometry projection, bed support gate, and direct
        observation anchors are member-local.  The optional bed/SMB inference
        routines retain ensemble-valued temporal/reference state and must see
        the complete ensemble to reproduce execution mode 1 exactly.
        """
        if not bool(icesee_kwargs.get("inference_plugin_enabled", False)):
            return False
        return bool(
            icesee_kwargs.get("physics_bed_inference", False)
            or icesee_kwargs.get("physics_smb_inference", False)
        )

    def _finalize_analysis_member_slabs(
        self, canonical, temp_path, timestep, icesee_kwargs
    ):
        """Finalize member-local physics with bounded rank-zero memory."""
        from ICESEE.src.parallelization._parallel_i_o import (
            finalize_analysis_ensemble,
        )

        budget = float(
            icesee_kwargs.get("analysis_finalize_memory_budget_mb", 2048.0)
        ) * 1024.0 * 1024.0
        # Forecast, analysis, and finalized float64 work arrays may coexist.
        bytes_per_member = max(1, 3 * self.nd * np.dtype(np.float64).itemsize)
        derived = max(1, int(budget // bytes_per_member))
        requested = int(
            icesee_kwargs.get("analysis_finalize_member_chunk_size", 0) or 0
        )
        member_chunk = min(
            self.nens,
            requested if requested > 0 else derived,
        )
        if member_chunk < 1 or bytes_per_member > budget:
            raise MemoryError(
                "One member cannot fit within the exact analysis-finalization "
                f"budget: approximately {bytes_per_member / 2**30:.2f} GiB "
                "is required. Increase analysis_finalize_memory_budget_mb or "
                "provide a model-native spatially partitioned finalizer."
            )

        with h5py.File(canonical, "r") as forecast_handle, h5py.File(
            temp_path, "r+"
        ) as analysis_handle:
            forecast_states = forecast_handle["states"]
            analysis_states = analysis_handle["states"]
            for member0 in range(0, self.nens, member_chunk):
                member1 = min(self.nens, member0 + member_chunk)
                forecast = np.asarray(
                    forecast_states[:, member0:member1], dtype=np.float64
                )
                analysis = np.asarray(
                    analysis_states[:, member0:member1], dtype=np.float64
                )
                # Member-local finalization must not leak mutable diagnostic
                # state from one slab into the next.
                slab_kwargs = dict(icesee_kwargs)
                finalized = finalize_analysis_ensemble(
                    analysis_vec=analysis,
                    forecast_vec=forecast,
                    timestep=timestep,
                    icesee_kwargs=slab_kwargs,
                )
                analysis_states[:, member0:member1] = np.asarray(
                    finalized, dtype=self.storage_dtype
                )
            analysis_handle.flush()

    def _publish_analysis_file(self, timestep, temp_path, icesee_kwargs):
        """Finalize and atomically publish one completed analysis shard.

        The forecast shard is never modified until every enabled analysis
        control and physical projection has succeeded.  Ensemble-coupled
        finalization currently uses a bounded rank-zero gather; exceeding the
        configured budget is an explicit unsupported capability rather than a
        silent change of algorithm.
        """
        canonical = self._state_file_name(timestep)
        self._close_batch()
        self.current_batch_start = -1
        self.mpi_comm.Barrier()

        error = None
        if self.rank == 0:
            try:
                if self._requires_shared_finalization(icesee_kwargs):
                    if not self._requires_ensemble_coupled_finalization(
                        icesee_kwargs
                    ):
                        self._finalize_analysis_member_slabs(
                            canonical, temp_path, timestep, icesee_kwargs
                        )
                        return_error = None
                    else:
                        return_error = "ensemble-coupled"
                else:
                    return_error = None

                if return_error == "ensemble-coupled":
                    itemsize = np.dtype(self.storage_dtype).itemsize
                    required = 2 * self.nd * self.nens * itemsize
                    budget = float(icesee_kwargs.get(
                        "analysis_finalize_memory_budget_mb", 2048.0
                    )) * 1024.0 * 1024.0
                    if required > budget:
                        raise MemoryError(
                            "Exact mode-1 ensemble-coupled inference requires "
                            "forecast and analysis ensembles together, but "
                            "their estimated "
                            f"working set is {required / 2**30:.2f} GiB and "
                            f"analysis_finalize_memory_budget_mb permits "
                            f"{budget / 2**30:.2f} GiB. Increase the budget or "
                            "disable the ensemble-coupled hook; ICESEE will not "
                            "silently publish a non-equivalent mode-2 analysis."
                        )
                    from ICESEE.src.parallelization._parallel_i_o import (
                        finalize_analysis_ensemble,
                    )
                    with h5py.File(canonical, "r") as forecast_handle:
                        forecast = np.asarray(forecast_handle["states"])
                    with h5py.File(temp_path, "r+") as analysis_handle:
                        analysis = np.asarray(analysis_handle["states"])
                        analysis = finalize_analysis_ensemble(
                            analysis_vec=analysis,
                            forecast_vec=forecast,
                            timestep=timestep,
                            icesee_kwargs=icesee_kwargs,
                        )
                        analysis_handle["states"][:, :] = analysis
                        analysis_handle.flush()
            except Exception as exc:
                error = f"{type(exc).__name__}: {exc}"

        error = self.mpi_comm.bcast(error, root=0)
        self.mpi_comm.Barrier()
        if error is not None:
            raise RuntimeError(
                "Transactional full-parallel analysis was not published: " + error
            )

        # The temporary shard is now physically valid but remains unpublished.
        # Hybrid inversion is applied member by member to this shard.  A failed
        # inversion therefore leaves the canonical forecast and restart point
        # untouched.
        if bool(icesee_kwargs.get("inversion_flag", False)):
            self._apply_memberwise_inversion(temp_path, timestep, icesee_kwargs)

        error = None
        if self.rank == 0:
            try:
                os.replace(temp_path, canonical)
            except Exception as exc:
                error = f"{type(exc).__name__}: {exc}"
        error = self.mpi_comm.bcast(error, root=0)
        self.mpi_comm.Barrier()
        if error is not None:
            raise RuntimeError(
                "Transactional full-parallel analysis was not published: " + error
            )

    def _apply_memberwise_inversion(self, temp_path, timestep, icesee_kwargs):
        """Apply mode-1-equivalent inversion without broadcasting the ensemble."""
        info = _inversion_variable_info(self.nd, icesee_kwargs)
        model_module = icesee_kwargs.get("model_module")
        inverse_step = getattr(model_module, "inverse_step_single", None)
        if inverse_step is None:
            raise RuntimeError(
                "inversion_flag is active, but model_module.inverse_step_single "
                "is unavailable."
            )

        local_updates = []
        local_error = None
        try:
            # Independent read-only member access avoids the mode-1 full-state
            # broadcast. Each rank retains only one member plus its recovered
            # friction values at a time.
            with h5py.File(temp_path, "r") as handle:
                states = handle["states"]
                for ens_id in range(self.rank, self.nens, self.size):
                    member = np.asarray(states[:, ens_id], dtype=np.float64)
                    inv_kwargs = dict(icesee_kwargs)
                    inv_kwargs.update({
                        "comm": MPI.COMM_SELF,
                        "ens_id": ens_id,
                        "nd": self.nd,
                        "nd_old": self.nd,
                        "vec_inputs": list(
                            icesee_kwargs.get("vec_inputs_old")
                            or icesee_kwargs.get("vec_inputs")
                        ),
                    })
                    result = inverse_step(ensemble=member, **inv_kwargs)
                    value = None
                    for key, candidate in result.items():
                        if str(key).lower() in {
                            "coefficient", "friction", "friction_coefficient",
                            "fcoef", "frictioncoefficient",
                            str(info["name"]).lower(),
                        }:
                            value = np.asarray(candidate, dtype=np.float64).reshape(-1)
                            break
                    if value is None:
                        raise KeyError(
                            "inverse_step_single did not return a recognized "
                            "friction/coefficient field."
                        )
                    if value.size != info["stop"] - info["start"]:
                        raise ValueError(
                            f"Inversion returned {value.size} friction values; "
                            f"expected {info['stop'] - info['start']}."
                        )
                    local_updates.append((ens_id, value))
        except Exception as exc:
            local_error = (
                f"rank {self.rank}: {type(exc).__name__}: {exc}\n"
                + traceback.format_exc()
            )

        errors = self.mpi_comm.allgather(local_error)
        failures = [message for message in errors if message]
        if failures:
            raise RuntimeError(
                "Member-wise hybrid inversion failed before publication:\n"
                + "\n".join(failures)
            )

        self.mpi_comm.Barrier()
        if self.size == 1:
            handle = h5py.File(temp_path, "r+")
        else:
            handle = h5py.File(
                temp_path, "r+", driver="mpio", comm=self.mpi_comm
            )
        try:
            states = handle["states"]
            for ens_id, value in local_updates:
                states[info["start"]:info["stop"], ens_id] = value
            handle.flush()
        finally:
            handle.close()
        self.mpi_comm.Barrier()

    def _create_file_collective(self, fname):
        import h5py, h5py.h5p as h5p, h5py.h5f as h5f, h5py.h5fd as h5fd
        fapl = h5p.create(h5p.FILE_ACCESS)
        fapl.set_fapl_mpio(self.comm, MPI.Info.Create())
        # collective metadata flags
        fapl.set_all_coll_metadata_ops(1)
        fapl.set_coll_metadata_write(1)
        fid = h5f.create(bytes(fname, 'utf-8'), flags=h5f.ACC_TRUNC, fapl=fapl)
        return h5py.File(fid)

    def _create_batch_serial(self, t_start):
        try:
            self._close_batch()
            self.files = []
            self.datasets = []
            self.current_batch_start = t_start

            remaining = max(0, self.nt - t_start)
            if remaining == 0:
                return

            nfiles = min(self.batch_size, remaining)

            if self.mpi_comm.Get_rank() == 0:
                for t in range(t_start, t_start + nfiles):
                    fname = self._state_file_name(t)

                    # Fresh-start cleanup is performed once in __init__.  From
                    # this point onward an existing shard belongs to the
                    # current run and must never be truncated when the sliding
                    # window reopens it.
                    mode = "a" if os.path.exists(fname) else "w"
                    with h5py.File(fname, mode) as f:
                        row_chunk = max(1, min(self.h5_file_chunk_size, self.nd))
                        col_chunk = 1
                        if "states" not in f:
                            f.create_dataset(
                                "states",
                                shape=(self.nd, self.nens),
                                chunks=(row_chunk, col_chunk),
                                compression=self.h5_file_compression,
                                compression_opts=self.h5_file_compression_options,
                                dtype=self.storage_dtype,
                            )
                        elif tuple(f["states"].shape) != (self.nd, self.nens):
                            raise ValueError(
                                f"Restart shard {fname} has shape {f['states'].shape}; "
                                f"expected {(self.nd, self.nens)}. Use force_fresh_start."
                            )

            self.mpi_comm.Barrier()

            for t in range(t_start, t_start + nfiles):
                fname = self._state_file_name(t)
                f = h5py.File(fname, "a", driver="mpio", comm=self.mpi_comm)
                f.atomic = False
                self.files.append(f)
                self.datasets.append(f["states"])

            self.mpi_comm.Barrier()

        except Exception as e:
            print(f"Error occurred in _create_batch_serial: {e}")
            tb_str = "".join(traceback.format_exception(*sys.exc_info()))
            print(f"Traceback details:\n{tb_str}")
            self.mpi_comm.Abort(1)

    def _create_batch_parallel(self, t_start):
        try:
            self._close_batch()
            self.files = []
            self.datasets = []
            self.current_batch_start = t_start

            remaining = max(0, self.nt - t_start)
            if remaining == 0:
                return

            nfiles = min(self.batch_size, remaining)

            for t in range(t_start, t_start + nfiles):
                fname = self._state_file_name(t)

                exists = os.path.exists(fname) if self.rank == 0 else None
                exists = self.mpi_comm.bcast(exists, root=0)
                # As in the serial-creation path, reopening a shard during the
                # same run must preserve data regardless of restart settings.
                mode = "a" if exists else "w"
                f = h5py.File(fname, mode, driver="mpio", comm=self.mpi_comm)
                f.atomic = False

                row_chunk = max(1, min(self.h5_file_chunk_size, self.nd))
                col_chunk = 1

                if "states" not in f:
                    dset = f.create_dataset(
                        "states",
                        shape=(self.nd, self.nens),
                        chunks=(row_chunk, col_chunk),
                        compression=self.h5_file_compression,
                        compression_opts=self.h5_file_compression_options,
                        dtype=self.storage_dtype,
                    )
                else:
                    dset = f["states"]
                    if tuple(dset.shape) != (self.nd, self.nens):
                        raise ValueError(
                            f"Restart shard {fname} has shape {dset.shape}; "
                            f"expected {(self.nd, self.nens)}. Use force_fresh_start."
                        )

                self.files.append(f)
                self.datasets.append(dset)

            self.mpi_comm.Barrier()

        except Exception as e:
            print(f"Error occurred in _create_batch_parallel: {e}")
            tb_str = "".join(traceback.format_exception(*sys.exc_info()))
            print(f"Traceback details:\n{tb_str}")
            self.mpi_comm.Abort(1)

    def _make_dxpl(self):
        import h5py.h5p as h5p, h5py.h5fd as h5fd
        dxpl = h5p.create(h5p.DATASET_XFER)
        mode = h5fd.MPIO_COLLECTIVE if self.use_collective_io else h5fd.MPIO_INDEPENDENT
        dxpl.set_dxpl_mpio(mode)
        return dxpl

    def _close_batch(self):
        try:
            for f in self.files:
                # try:
                #     f.flush()
                # except Exception:
                #     pass
                f.close()
            self.files = []
            self.datasets = []
        except Exception as e:
            print(f"Error occurred in _close_batch: {e}")
            tb_str = "".join(traceback.format_exception(*sys.exc_info()))
            print(f"Traceback details:\n{tb_str}")
            self.mpi_comm.Abort(1)

    def _ensure_batch(self, t):
        try:
            t = int(t)
            batch_end = self.current_batch_start + len(self.datasets)
            if not (self.current_batch_start <= t < batch_end):
                # A sliding window keeps the input and output timestep open at
                # once.  Unlike aligned batches, it does not thrash at every
                # batch boundary during read(k) -> write(k+1).
                batch_start = max(0, t - self.batch_size + 1)
                # self._create_batch(batch_start)
                if self.serial_file_creation:
                    self._create_batch_serial(batch_start)
                else:
                    self._create_batch_parallel(batch_start)
        except Exception as e:
            print(f"Error occurred in _ensure_batch: {e}")
            tb_str = "".join(traceback.format_exception(*sys.exc_info()))
            print(f"Traceback details:\n{tb_str}")
            self.mpi_comm.Abort(1)

    def prune_history(self, keep_t):
        """Retain only ``keep_t`` when rolling ensemble history is enabled.

        The ensemble mean time series is stored separately and is not removed.
        All ranks close the active HDF5 window before rank zero unlinks older
        shards, making the retained shard a complete restart point.
        """
        if self.history_mode != "rolling":
            return
        keep_t = int(keep_t)
        self._close_batch()
        self.current_batch_start = -1
        self.mpi_comm.Barrier()
        if self.rank == 0:
            keep_name = os.path.abspath(self._state_file_name(keep_t))
            pattern = f"{self.base_path}/{self.file_prefix}_*.h5"
            for path in glob.glob(pattern):
                # The mean file shares the prefix but is not a state shard.
                if os.path.abspath(path) == keep_name or path.endswith("_mean.h5"):
                    continue
                if re.search(r"_\d+\.h5$", os.path.basename(path)):
                    try:
                        os.remove(path)
                    except FileNotFoundError:
                        pass
        self.mpi_comm.Barrier()

    def _rw_select(self, ds, start_row, nrows, col0, ncols, buf=None, write=False):
        import numpy as np, h5py.h5s as h5s
        file_space = ds.id.get_space()
        file_space.select_hyperslab((start_row, col0), (nrows, ncols))
        mem_shape = (nrows,) if ncols == 1 else (nrows, ncols)
        mem_space = h5s.create_simple(mem_shape)
        dxpl = self._make_dxpl()
        if write:
            ds.id.write(mem_space, file_space, np.ascontiguousarray(buf), dxpl=dxpl)
        else:
            out = np.empty(mem_shape, dtype=self.storage_dtype, order='C')
            ds.id.read(mem_space, file_space, out, dxpl=dxpl)
            return out

    @retry_on_failure(max_attempts=3, delay=0.5, mpi_comm=MPI.COMM_WORLD)
    def read_forecast(self, t, ens_idx):
        """
        Forecast path: return the full ensemble vector.

        This matches _mpi_forecast_functions.py where the MPI model receives
        a complete state vector for one ensemble member.
        """
        self._ensure_batch(t)
        batch_idx = t - self.current_batch_start

        data = self.datasets[batch_idx][:, ens_idx]
        return np.asarray(data, dtype=np.float64).reshape(-1)


    @retry_on_failure(max_attempts=3, delay=0.5, mpi_comm=MPI.COMM_WORLD)
    def write_forecast(self, t, data, ens_idx):
        """
        Forecast path: write one full ensemble vector.

        This should be called only by sub_rank == 0 from
        parallel_forecast_step_default_full_parallel_run.
        """
        self._ensure_batch(t)
        batch_idx = t - self.current_batch_start

        data = np.asarray(data, dtype=self.storage_dtype).reshape(-1)

        if data.size != self.nd:
            raise ValueError(
                f"write_forecast expected full vector of size {self.nd}, "
                f"got {data.size}"
            )

        self.datasets[batch_idx][:, ens_idx] = data

    def write_forecast_from_root(self, t, data, ens_idx, root=0):
        """Collectively write one member initialized on ``root``.

        Sequential model initialization is useful when a model instance is
        expensive or not safe to construct concurrently.  The old mode-2 path
        let rank zero call :meth:`write_forecast` alone, which deadlocked in
        the collective batch-opening barriers.  This routine keeps only one
        full member on the root, scatters disjoint row slabs, and has every
        world rank participate in the parallel-HDF5 write.
        """
        self._ensure_batch(t)
        batch_idx = int(t) - self.current_batch_start
        counts = np.asarray(
            [
                (self.nd // self.size) + (1 if rank < self.nd % self.size else 0)
                for rank in range(self.size)
            ],
            dtype=np.int64,
        )
        displacements = np.zeros(self.size, dtype=np.int64)
        if self.size > 1:
            displacements[1:] = np.cumsum(counts[:-1])

        if self.rank == root:
            source = np.asarray(data, dtype=np.float64).reshape(-1)
            if source.size != self.nd:
                raise ValueError(
                    f"write_forecast_from_root expected {self.nd} values, "
                    f"got {source.size}"
                )
        else:
            source = None

        local = np.empty(int(counts[self.rank]), dtype=np.float64)
        self.mpi_comm.Scatterv(
            [source, counts, displacements, MPI.DOUBLE] if self.rank == root else None,
            local,
            root=root,
        )
        if local.size:
            self._rw_select(
                self.datasets[batch_idx],
                int(displacements[self.rank]),
                int(local.size),
                int(ens_idx),
                1,
                buf=local.astype(self.storage_dtype, copy=False),
                write=True,
            )
        self.mpi_comm.Barrier()

    @retry_on_failure(max_attempts=3, delay=0.5, mpi_comm=MPI.COMM_WORLD)
    def read_analysis(self, t, ens_idx):
        """
        Analysis path: return this world-rank's row block.
        """
        self._ensure_batch(t)
        batch_idx = t - self.current_batch_start

        data = self.datasets[batch_idx][
            self.nd_start_world:self.nd_end_world,
            ens_idx,
        ]

        return np.asarray(data, dtype=np.float64).reshape(-1)


    @retry_on_failure(max_attempts=3, delay=0.5, mpi_comm=MPI.COMM_WORLD)
    def write_analysis(self, t, data, ens_idx):
        """
        Analysis path: write this world-rank's row block.
        """
        self._ensure_batch(t)
        batch_idx = t - self.current_batch_start

        data = np.asarray(data, dtype=self.storage_dtype).reshape(-1)

        expected = self.nd_local_world
        if data.size != expected:
            raise ValueError(
                f"write_analysis rank {self.rank} expected local vector size "
                f"{expected}, got {data.size}"
            )

        self.datasets[batch_idx][
            self.nd_start_world:self.nd_end_world,
            ens_idx,
        ] = data

    def compute_forecast_mean_chunked_v2(self, k, flag=None):
        """
        Simple & hang-free:
        - running sum in RAM (length = local_rows)
        - collective dataset creation
        - collective write with empty selection for zero-row ranks
        """
        from mpi4py import MPI
        import numpy as np, h5py
        import h5py.h5p as h5p, h5py.h5s as h5s, h5py.h5fd as h5fd
        import sys, traceback

        comm = self.mpi_comm
        rank = comm.Get_rank()
        size = comm.Get_size()

        nt = self.nt

        try:
            nd0, nd1 = self.nd_start_world, self.nd_end_world
            local_rows = nd1 - nd0
            if local_rows < 0:
                raise ValueError("Invalid local row bounds")

            # ---- Optional per-rank Zarr cache (safe; not required) ---------------
            # If you keep this, it's fine—just per-rank path to avoid contention.
            # import zarr
            # zarr.create_array(f"{self.base_path}/{self.file_prefix}_forecast_updates_{rank}.zarr",
            #                   shape=(local_rows, self.nens),
            #                   chunks=(min(local_rows, 1000), 1), dtype='f8', overwrite=True)

            # ---- Running sum while reading ensembles ------------------------------
            local_sum = np.zeros(max(local_rows, 0), dtype='f8')
            batch_idx = k - self.current_batch_start
            for ens_idx in range(self.nens):
                if local_rows > 0:
                    # v = self.read_analysis(k, ens_idx)
                    v = self.datasets[batch_idx][self.nd_start_world:self.nd_end_world, ens_idx]
                    v = np.asarray(v, dtype='f8')
                    if v.ndim != 1 or v.size != local_rows:
                        v = v.reshape(-1)
                        assert v.size == local_rows, "read_analysis must return (local_rows,)"
                    local_sum += v

            local_mean = (local_sum / float(self.nens)) if local_rows > 0 else np.empty((0,), dtype='f8')

            # ---- Parallel HDF5: collective create + collective write -------------
            file_path = f"{self.base_path}/{self.file_prefix}_mean.h5"

            if flag == 'initial':
                if rank == 0 and k==0: # remove old file if any
                    try:
                        os.remove(file_path)
                    except OSError:
                        pass
                comm.Barrier()

            # Dataset transfer property list: COLLECTIVE for the write
            dxpl = h5p.create(h5p.DATASET_XFER)
            dxpl.set_dxpl_mpio(h5fd.MPIO_COLLECTIVE)

            with h5py.File(file_path, 'a', driver='mpio', comm=comm) as f:

                # --- Collective dataset creation: ALL ranks must take the same branch
                exists_local = ('mean' in f)
                # If any rank sees it, treat as exists for all to avoid split branches
                exists_any = comm.allreduce(1 if exists_local else 0, op=MPI.SUM) > 0

                if not exists_any:
                    # ALL ranks call create_dataset with identical args (collective)
                    chunk_rows = min(self.nd, 4096)
                    f.create_dataset('mean', (self.nd, self.nt),
                                    chunks=(chunk_rows, 1),
                                    dtype=self.storage_dtype)
                    # Ensure all ranks see the new metadata
                    comm.Barrier()
                else:
                    # Ensure all ranks take the same path
                    comm.Barrier()

                dset = f['mean']

                # --- Collective write: ALL ranks must participate
                file_space = dset.id.get_space()
                if local_rows > 0:
                    # Select this rank's row slab for column k
                    file_space.select_hyperslab((nd0, k), (local_rows, 1))
                    mem_space = h5s.create_simple((local_rows,))
                    buf = np.ascontiguousarray(local_mean)
                else:
                    # Empty (NULL) selection for ranks with no rows
                    mem_space = h5s.create_simple((0,))
                    file_space.select_none()
                    buf = np.empty((0,), dtype='f8')

                dset.id.write(mem_space, file_space, buf, dxpl=dxpl)

            comm.Barrier()

        except Exception as e:
            print(f"Error in compute_forecast_mean_chunked_v2: {e}")
            tb_str = "".join(traceback.format_exception(*sys.exc_info()))
            print(f"Traceback details:\n{tb_str}")
            self.mpi_comm.Abort(1)

    def compare_forecast_means(self, t, ens_chunk_size=1, rtol=1e-5, atol=1e-8):
        """Debug check without materializing the state-by-ensemble matrix."""
        try:
            comm = self.mpi_comm
            rank = comm.Get_rank()

            mean_original = None
            if rank == 0:
                file_path = f"{self.base_path}/{self.file_prefix}_mean.h5"
                self._ensure_batch(t)
                batch_idx = t - self.current_batch_start
                mean_original = np.empty(self.nd, dtype=np.float64)
                rows = max(1, resolve_row_chunk_size(
                    self.nd,
                    self.nens,
                    explicit_rows=int(self.icesee_kwargs.get(
                        "analysis_row_chunk_size", 0
                    )),
                    memory_budget_mb=float(self.icesee_kwargs.get(
                        "analysis_memory_budget_mb", 256.0
                    )),
                ))
                for row0 in range(0, self.nd, rows):
                    row1 = min(self.nd, row0 + rows)
                    block = np.asarray(
                        self.datasets[batch_idx][row0:row1, :],
                        dtype=np.float64,
                    )
                    mean_original[row0:row1] = np.mean(block, axis=1)

            self.compute_forecast_mean_chunked_v2(t)
            mean_chunked = None
            if rank == 0:
                with h5py.File(file_path, 'r') as f:
                    mean_chunked = f['mean'][:, t].copy()

            if rank == 0:
                if mean_original is None or mean_chunked is None:
                    print("Error: One or both means could not be read from file.")
                    return False

                are_equal = np.allclose(mean_original, mean_chunked, rtol=rtol, atol=atol)
                if are_equal:
                    print(f"Means for time step {t} are equivalent within tolerance (rtol={rtol}, atol={atol}).")
                else:
                    print(f"Means for time step {t} differ beyond tolerance (rtol={rtol}, atol={atol}).")
                    max_diff = np.max(np.abs(mean_original - mean_chunked))
                    print(f"Maximum absolute difference: {max_diff}")
                return are_equal
            return None
        except Exception as e:
            print(f"Error occurred in compare_forecast_means: {e}")
            tb_str = "".join(traceback.format_exception(*sys.exc_info()))
            print(f"Traceback details:\n{tb_str}")
            self.mpi_comm.Abort(1)

    @retry_on_failure(max_attempts=5, delay=1.0, mpi_comm=MPI.COMM_WORLD)
    def generate_observation_schedule(self, **icesee_kwargs):
        import numpy as np

        t = np.asarray(icesee_kwargs["t"], dtype=float)
        if t.ndim != 1 or t.size == 0:
            raise ValueError("`t` must be a 1D non-empty array of times.")
        t_min, t_max = float(t[0]), float(t[-1])

        freq_obs = float(self.icesee_kwargs["freq_obs"])
        obs_start = float(self.icesee_kwargs["obs_start_time"])
        obs_max_cfg = float(self.icesee_kwargs["obs_max_time"])

        obs_start = max(obs_start, t_min)
        obs_max = min(obs_max_cfg, t_max)

        if freq_obs <= 0.0 or obs_start > obs_max:
            return np.array([]), np.array([], dtype=int), 0

        # --- Build ideal observation times (requested) ---
        n_obs = int(np.floor((obs_max - obs_start) / freq_obs)) + 1
        obs_t_req = obs_start + np.arange(n_obs, dtype=float) * freq_obs

        # --- Match model time points to requested obs times (same as partial) ---
        dt_grid = float(np.min(np.diff(t))) if t.size > 1 else 1.0

        obs_idx = []
        for tobs in obs_t_req:
            i = int(np.argmin(np.abs(t - tobs)))
            # accept the nearest model time if it's close enough
            if abs(t[i] - tobs) <= 0.5 * dt_grid:
                obs_idx.append(i)

        obs_idx = np.array(sorted(set(obs_idx)), dtype=int)
        num_observations = int(obs_idx.size)

        return obs_t_req, obs_idx, num_observations

    @retry_on_failure(max_attempts=5, delay=1.0, mpi_comm=MPI.COMM_WORLD)
    def _create_synthetic_observations(self, **icesee_kwargs):
        import os
        import re
        import h5py
        import numpy as np
        import traceback
        import sys

        nd = self.nd
        nt = self.nt

        obs_t, ind_m, m_obs = self.generate_observation_schedule(**icesee_kwargs)
        ind_m = np.asarray(ind_m, dtype=int)  # your code treats these as step indices
        obs_t = np.asarray(obs_t, dtype=float)

        rank = self.mpi_comm.Get_rank()

        total_state_param_vars = self.icesee_kwargs["total_state_param_vars"]
        hdim = nd // total_state_param_vars

        # ----------------------------
        # Bed snapshots: interpret as YEARS like partial-parallel
        # ----------------------------
        bed_snaps = np.asarray(icesee_kwargs.get("bed_obs_snapshot", []), dtype=float)
        bed_snap_cols = []
        bed_time_to_col = {}

        if bed_snaps.size > 0 and obs_t.size > 0:
            for bed_time in bed_snaps:
                diffs = np.abs(obs_t - bed_time)
                j = int(np.argmin(diffs))  # 0-based obs column
                bed_snap_cols.append(j)
                bed_time_to_col[bed_time] = j
            bed_snap_cols = sorted(set(bed_snap_cols))

        # The snapshot columns are analysis metadata, not a writer-local
        # implementation detail.  Persist them in both dictionaries before
        # rank 0 writes the observations so every later analysis/restart path
        # sees the same gate schedule as execution mode 1.
        icesee_kwargs["bed_snap_cols"] = list(bed_snap_cols)
        self.icesee_kwargs["bed_snap_cols"] = list(bed_snap_cols)

        if rank == 0:
            try:
                print("[ICESEE] Generating synthetic observations ...")
                print(f"[ICESEE] observation times: {obs_t}, indices: {ind_m}, total: {m_obs}")
                print(f"[ICESEE] bed_snaps (years): {bed_snaps}")
                print(f"[ICESEE] bed_snap_cols (obs columns): {bed_snap_cols}")
                if len(bed_snap_cols) > 0:
                    print(f"[ICESEE] obs_t at bed_snap_cols: {obs_t[bed_snap_cols]}")

                obs_file = f"{self.base_path}/synthetic_obs.h5"

                # Build index map once
                vec_inputs = list(icesee_kwargs["vec_inputs"])
                indx_map = _global_variable_index_map(
                    nd, vec_inputs, icesee_kwargs.get("var_nd")
                )

                # Helpers / options
                num_state_vars = icesee_kwargs.get("num_state_vars", self.icesee_kwargs.get("num_state_vars"))
                observed_params = set(icesee_kwargs.get("observed_params", []))

                bed_aliases = {"bed", "bedrock", "bed_topography", "bedtopo", "bedtopography"}
                key_is_bed = {k: (k in bed_aliases) for k in vec_inputs}
                key_idx_map = {k: np.asarray(indx_map[k], dtype=int) for k in vec_inputs}

                # ---- Bed sparsity controls (accept BOTH naming conventions safely) ----
                Ly = icesee_kwargs.get("Ly", self.icesee_kwargs.get("Ly", None))
                Lx = icesee_kwargs.get("Lx", self.icesee_kwargs.get("Lx", None))
                model_name = icesee_kwargs.get("model_name", None)

                # allow either key name:
                bed_stride_km = icesee_kwargs.get("bed_obs_stride", None)
                if bed_stride_km is None:
                    bed_stride_km = icesee_kwargs.get("bed_obs_stride_km", None)

                bed_spacing_pts = icesee_kwargs.get("bed_obs_spacing", None)
                bed_indices_user = icesee_kwargs.get("bed_obs_indices", None)
                bed_mask_user = icesee_kwargs.get("bed_obs_mask", None)

                # ---- Build bed_mask_map using the SAME priority/shape handling as partial ----
                bed_mask_map = {}
                for k in vec_inputs:
                    if not key_is_bed.get(k, False):
                        continue

                    idx = key_idx_map[k]
                    local_len = idx.size

                    mask = np.ones(local_len, dtype=bool)  # default observe all

                    # Priority 1: explicit mask
                    if isinstance(bed_mask_user, (list, np.ndarray)):
                        msk = np.asarray(bed_mask_user, dtype=bool)
                        if msk.ndim > 1 and msk.shape[0] == 1:
                            msk = msk[0]
                        msk = msk.ravel()

                        if msk.size != local_len:
                            if msk.ndim == 1 and msk.size > 0:
                                rep = int(np.ceil(local_len / msk.size))
                                msk = np.tile(msk, rep)[:local_len]
                            else:
                                msk = np.ones(local_len, dtype=bool)

                        mask = msk

                    # Priority 2: explicit indices
                    elif isinstance(bed_indices_user, (list, np.ndarray)):
                        mask = np.zeros(local_len, dtype=bool)
                        idxs = np.asarray(bed_indices_user, dtype=int)
                        idxs = idxs[(idxs >= 0) & (idxs < local_len)]
                        mask[idxs] = True

                    # Priority 3: spacing in points
                    elif isinstance(bed_spacing_pts, (int, np.integer)) and bed_spacing_pts > 1:
                        n = int(bed_spacing_pts)
                        mask = np.zeros(local_len, dtype=bool)
                        mask[::n] = True

                    # Priority 4: LiDAR-like / stride in km (2D grid or ISSM mesh)
                    elif (bed_stride_km is not None) and (Lx is not None) and (Ly is not None):
                        if re.match(r"(?i)^issm$", str(model_name)):
                            import h5py  # already imported, but harmless

                            icesee_path = icesee_kwargs.get("icesee_path")
                            data_path = icesee_kwargs.get("data_path")

                            file_path = f"{icesee_path}/{data_path}/mesh_idxy_{0}.h5"
                            try:
                                with h5py.File(file_path, "r") as f:
                                    x_param = f["/fric_x"][:]
                                    y_param = f["/fric_y"][:]
                            except FileNotFoundError:
                                raise FileNotFoundError(
                                    f"ISSM mesh file '{file_path}' not found. "
                                    "Please generate the mesh indicies before running ICESEE."
                                )

                            y_param = np.asarray(y_param / 1000.0, dtype=float).reshape(-1)
                            x_param = np.asarray(x_param / 1000.0, dtype=float).reshape(-1)

                            y_min, y_max = np.min(y_param), np.max(y_param)
                            x_min, x_max = np.min(x_param), np.max(x_param)

                            local_len = x_param.size
                            bed_stride_km_local = float(bed_stride_km) / 1000.0  # keep your partial conversion

                            x_lines = np.arange(x_min, x_max + 1e-6, bed_stride_km_local)

                            if x_lines.size > 1:
                                dx_nom = (x_max - x_min) / (x_lines.size - 1)
                            else:
                                dx_nom = bed_stride_km_local
                            band = 0.5 * dx_nom

                            mask = np.zeros(local_len, dtype=bool)
                            for x_line in x_lines:
                                mask |= np.abs(x_param - x_line) <= band

                            print(f"[ICESEE<-ISSM] bed LiDAR mask for '{k}': {mask.sum()} of {local_len} points observed")

                        else:
                            # 2D grid assumption like partial
                            Nx = int(hdim)
                            Ny = int(local_len // Nx) if Nx > 0 else 1

                            if Nx * Ny != local_len:
                                intervals = max(hdim - 1, 1)
                                dx = float(Lx) / intervals
                                n = max(int(round(float(bed_stride_km) / max(dx, 1e-12))), 1)
                                mask = np.zeros(local_len, dtype=bool)
                                mask[::n] = True
                            else:
                                intervals_y = max(Ny - 1, 1)
                                dy = float(Ly) / intervals_y
                                stride_y_pts = max(int(round(float(bed_stride_km) / max(dy, 1e-12))), 1)

                                mask2d = np.zeros((Ny, Nx), dtype=bool)
                                for j in range(0, Ny, stride_y_pts):
                                    mask2d[j, :] = True  # keep all along-track points by default

                                mask = mask2d.ravel(order="C")
                                print(f"[ICESEE] bed 2D mask for key '{k}': Ny={Ny}, Nx={Nx}, stride_y_pts={stride_y_pts}")

                    bed_mask_map[k] = mask

                # H_matrix and the streamed analysis must see the same masks.
                self.icesee_kwargs["bed_mask_map"] = bed_mask_map
                icesee_kwargs["bed_mask_map"] = bed_mask_map

                # ---- Create / overwrite output datasets ----
                with h5py.File(obs_file, "a") as f:
                    for name in (
                        "hu_obs", "error_R", "hu_obs_compact", "obs_indices",
                        "obs_std", "obs_active", "obs_index", "obs_t", "obs_max_time",
                        "bed_snap_cols",
                    ):
                        if name in f:
                            del f[name]

                    sig_obs = np.asarray(self.icesee_kwargs["sig_obs"]).reshape(-1)
                    if sig_obs.size < len(vec_inputs):
                        sig_obs = np.pad(
                            sig_obs, (0, len(vec_inputs) - sig_obs.size), mode="edge"
                        )

                    observed = list(dict.fromkeys(
                        list(icesee_kwargs.get("observed_vars", []))
                        + list(icesee_kwargs.get("observed_params", []))
                    ))
                    obs_indices_parts = []
                    obs_std_parts = []
                    obs_is_bed_parts = []
                    for key in observed:
                        if key not in key_idx_map:
                            raise KeyError(f"Observed key '{key}' is not in vec_inputs")
                        idx = key_idx_map[key]
                        bed_flag = key_is_bed.get(key, False)
                        if bed_flag:
                            idx = idx[bed_mask_map.get(key, np.ones(idx.size, dtype=bool))]
                        sigma = float(sig_obs[vec_inputs.index(key)])
                        obs_indices_parts.append(idx)
                        obs_std_parts.append(np.full(idx.size, sigma, dtype=np.float64))
                        obs_is_bed_parts.append(np.full(idx.size, bed_flag, dtype=bool))

                    obs_indices = np.concatenate(obs_indices_parts).astype(np.int64)
                    obs_std = np.concatenate(obs_std_parts)
                    obs_is_bed = np.concatenate(obs_is_bed_parts)
                    nobs = int(obs_indices.size)
                    col_chunk = max(1, min(50, m_obs))
                    row_chunk = max(1, min(
                        nobs,
                        resolve_row_chunk_size(
                            nobs, max(1, m_obs),
                            budget_mb=icesee_kwargs.get(
                                "analysis_memory_budget_mb", 256.0
                            ),
                            explicit_rows=icesee_kwargs.get(
                                "observation_row_chunk_size", 0
                            ),
                            working_buffers=2,
                        ),
                    ))
                    if m_obs:
                        hu_obs = f.create_dataset(
                            "hu_obs_compact", (nobs, m_obs),
                            chunks=(row_chunk, col_chunk), dtype="f8"
                        )
                        active = f.create_dataset(
                            "obs_active", (nobs, m_obs),
                            chunks=(row_chunk, col_chunk), dtype="bool"
                        )
                    else:
                        # HDF5 chunks may not exceed a zero-length axis.  A
                        # contiguous empty dataset preserves mode-1 semantics
                        # for shortened/no-observation benchmark runs.
                        hu_obs = f.create_dataset(
                            "hu_obs_compact",
                            data=np.empty((nobs, 0), dtype=np.float64),
                        )
                        active = f.create_dataset(
                            "obs_active",
                            data=np.empty((nobs, 0), dtype=bool),
                        )
                    f.create_dataset("obs_indices", data=obs_indices, dtype="i8")
                    f.create_dataset("obs_std", data=obs_std, dtype="f8")
                    f.create_dataset("obs_index", data=ind_m, dtype="i8")
                    f.create_dataset("obs_t", data=obs_t, dtype="f8")
                    f.create_dataset(
                        "bed_snap_cols",
                        data=np.asarray(bed_snap_cols, dtype=np.int64),
                        dtype="i8",
                    )
                    f.create_dataset("obs_max_time", data=np.asarray([obs_t[-1] if obs_t.size else 0.0]))
                    f.attrs["observation_storage"] = "compact_rows"
                    f.attrs["state_dimension"] = int(nd)

                    # Synthetic observations must be an execution-mode
                    # invariant.  Do not inherit the process-global RNG (or a
                    # rank-local forecast RNG), because modes 1 and 2 consume
                    # those streams in a different order before reaching this
                    # routine.
                    observation_seed = int(
                        icesee_kwargs.get(
                            "synthetic_observation_seed",
                            int(
                                icesee_kwargs.get(
                                    "base_seed", icesee_kwargs.get("seed", 42)
                                )
                            )
                            + 2000003,
                        )
                    )
                    observation_rng = np.random.default_rng(observation_seed)
                    self.icesee_kwargs["synthetic_observation_seed"] = observation_seed
                    with h5py.File(
                        f"{self.base_path}/true_nurged_states.h5", "r"
                    ) as true_file:
                        truth = true_file["true_state"]
                        for km, step in enumerate(ind_m):
                            # Match the partial-parallel observation contract:
                            # ``ind_m`` contains zero-based columns of the
                            # stored reference trajectory.
                            tcol = int(step)
                            active_col = (~obs_is_bed) | np.isin(km, bed_snap_cols)
                            active[:, km] = active_col
                            for start in range(0, nobs, row_chunk):
                                stop = min(nobs, start + row_chunk)
                                use = active_col[start:stop]
                                values = np.zeros(stop - start, dtype=np.float64)
                                if np.any(use):
                                    rows = obs_indices[start:stop][use]
                                    values[use] = _read_rows_preserve_order(
                                        truth, rows, tcol
                                    ) + observation_rng.normal(
                                        0.0, obs_std[start:stop][use]
                                    )
                                hu_obs[start:stop, km] = values

                    if str(icesee_kwargs.get(
                        "synthetic_observation_storage", "compact"
                    )).lower() in {"dense", "both"}:
                        dense = f.create_dataset(
                            "hu_obs", (nd, m_obs),
                            chunks=(min(4096, nd), col_chunk), dtype="f8",
                            fillvalue=0.0,
                        )
                        for start in range(0, nobs, row_chunk):
                            stop = min(nobs, start + row_chunk)
                            dense[obs_indices[start:stop], :] = hu_obs[start:stop, :]
                        print(
                            "[ICESEE] Warning: dense synthetic-observation export "
                            "is enabled and scales with the full state dimension."
                        )

            except Exception as e:
                print(f"[ICESEE] Error in full-parallel _create_synthetic_observations: {e}")
                tb_str = "".join(traceback.format_exception(*sys.exc_info()))
                print(f"Traceback details:\n{tb_str}")
                raise

        self.mpi_comm.Barrier()
        return ind_m, m_obs

    def H_matrix(self, **icesee_kwargs):
        """
        Fully-parallel version: write H directly to Zarr, but use the SAME
        observation-index logic as the partial-parallel run (including bed masks).

        H contains ONLY real observation points (and bed subsampling/masking),
        with rows ordered exactly as `icesee_kwargs["all_observed"]` concatenation.
        """
        try:
            icesee_kwargs = self.icesee_kwargs
            observed = icesee_kwargs["all_observed"]  # e.g. ['h','u','v','smb','bed']
            vec_inputs = icesee_kwargs.get("vec_inputs", [])

            nd = int(self.nd)  # state size
            indx_map = _global_variable_index_map(
                nd, vec_inputs, icesee_kwargs.get("var_nd")
            )

            # --- Retrieve bed masks (must already exist) ---
            bed_mask_map = icesee_kwargs.get("bed_mask_map", {})

            bed_aliases = {"bed", "bedrock", "bed_topography", "bedtopo", "bedtopography"}
            key_is_bed = {k: (k in bed_aliases) for k in vec_inputs}

            compact_file = f"{self.base_path}/synthetic_obs.h5"
            if os.path.exists(compact_file):
                with h5py.File(compact_file, "r") as obs_handle:
                    if "obs_indices" in obs_handle:
                        obs_indices = np.asarray(
                            obs_handle["obs_indices"], dtype=np.int64
                        )
                    else:
                        obs_indices = None
            else:
                obs_indices = None

            if obs_indices is None:
                # Legacy observation files: reconstruct the selection operator.
                all_obs_indices = []
                for key in observed:
                    if key not in indx_map:
                        raise KeyError(f"Observed key '{key}' not found in indx_map")
                    idx = np.asarray(indx_map[key], dtype=int)
                    if key_is_bed.get(key, False) and key in bed_mask_map:
                        mask = np.asarray(bed_mask_map[key], dtype=bool).ravel()
                        if mask.size != idx.size:
                            raise ValueError(
                                f"bed_mask_map['{key}'] length {mask.size} does "
                                f"not match bed vector length {idx.size}"
                            )
                        idx = idx[mask]
                    all_obs_indices.append(idx)
                obs_indices = np.concatenate(all_obs_indices).astype(np.int64)

            if obs_indices.size == 0:
                raise ValueError("No observation indices found (obs_indices is empty).")

            # --- Safety ---
            if obs_indices.max() >= nd:
                raise ValueError(
                    f"H_matrix error: obs index {obs_indices.max()} >= state size {nd}"
                )

            # m = number of real observations (rows)
            m = int(obs_indices.size)

            # ---- Output path (HDF5) ----
            H_matrix_file = icesee_kwargs.get("H_matrix_h5_path", f"{self.base_path}/H_matrix.h5")

            rank = self.mpi_comm.Get_rank()
            if rank == 0:
                out_dir = os.path.dirname(H_matrix_file)
                if out_dir:
                    os.makedirs(out_dir, exist_ok=True)
                try:
                    os.remove(H_matrix_file)
                except FileNotFoundError:
                    pass

                with h5py.File(H_matrix_file, "w") as f:
                    f.create_dataset("obs_indices", data=obs_indices, dtype="i8")
                    f.attrs["state_dimension"] = nd
                    f.attrs["operator_format"] = "row_indices"

                    # Dense H is retained only as an opt-in compatibility export.
                    storage = str(icesee_kwargs.get(
                        "observation_operator_storage", "indices"
                    )).lower()
                    if storage in {"dense", "both"}:
                        row_chunk = max(1, min(
                            int(icesee_kwargs.get("H_row_chunk", min(512, m))), m
                        ))
                        col_chunk = max(1, min(
                            int(icesee_kwargs.get("H_col_chunk", min(2048, nd))), nd
                        ))
                        dset = f.create_dataset(
                            "H", shape=(m, nd), dtype="f8",
                            chunks=(row_chunk, col_chunk), fillvalue=0.0
                        )
                        for row, column in enumerate(obs_indices):
                            dset[row, column] = 1.0

            self.mpi_comm.Barrier()
        except Exception as e:
            print(f"Error in H_matrix: {e}")
            tb_str = "".join(traceback.format_exception(*sys.exc_info()))
            print(f"Traceback details:\n{tb_str}")
            self.mpi_comm.Abort(1)

    def compute_X5_utils_(self, **icesee_kwargs):
        """Accumulate exact mode-1 analysis products in bounded memory.

        Because ICESEE's observation operator selects state-vector rows, the
        dense ``m x nd`` matrix is unnecessary.  Observation rows are divided
        among MPI ranks and streamed from the timestep and observation files.
        Only three ``Nens x Nens`` products are globally reduced.  Standard
        stochastic observation perturbations are stored as an observation-
        sized temporary factor so every MPI partition consumes the exact same
        seeded sequence as mode 1 without replicating the state ensemble.
        """
        k = icesee_kwargs.get('k')
        k = k + 1 if k < self.nt - 1 else k
        km = icesee_kwargs.get('km')
        try:
            Nens = self.nens
            with h5py.File(f"{self.base_path}/H_matrix.h5", "r") as handle:
                obs_indices = np.asarray(handle["obs_indices"], dtype=np.int64)
            m = obs_indices.size
            inversion_keep_all = _exclude_inversion_observations(
                obs_indices, icesee_kwargs, self.nd
            )
            rank = self.mpi_comm.Get_rank()
            size = self.mpi_comm.Get_size()
            obs_start, obs_end = _partition_1d(m, size, rank)
            local_indices = obs_indices[obs_start:obs_end]
            chunk_rows = resolve_row_chunk_size(
                local_indices.size, Nens,
                icesee_kwargs.get("analysis_memory_budget_mb", 256),
                icesee_kwargs.get("observation_row_chunk_size", 0),
                working_buffers=4,
            )

            error_mode = str(icesee_kwargs.get(
                "enkf_observation_error_mode", "stochastic_R"
            )).lower()
            if error_mode not in {
                "stochastic_r", "legacy_prior_anomalies", "generated_r"
            }:
                raise ValueError(
                    "enkf_observation_error_mode must be stochastic_R, "
                    "legacy_prior_anomalies, or generated_R"
                )

            cross_local = np.zeros((Nens, Nens), dtype=np.float64)
            gram_local = np.zeros((Nens, Nens), dtype=np.float64)
            rhs_local = np.zeros((Nens, Nens), dtype=np.float64)
            batch_idx = k - self.current_batch_start
            state_dset = self.datasets[batch_idx]
            obs_file = f"{self.base_path}/synthetic_obs.h5"
            io_start = MPI.Wtime()
            with h5py.File(obs_file, "r") as obs_handle:
                compact = "hu_obs_compact" in obs_handle
                obs_dset = obs_handle[
                    "hu_obs_compact" if compact else "hu_obs"
                ]
                local_positions = np.arange(obs_start, obs_end, dtype=np.int64)
                if compact and "obs_active" in obs_handle:
                    keep = np.asarray(
                        obs_handle["obs_active"][obs_start:obs_end, int(km)],
                        dtype=bool,
                    )
                    local_positions = local_positions[keep]
                    local_indices = local_indices[keep]
                keep_inversion = inversion_keep_all[local_positions]
                local_positions = local_positions[keep_inversion]
                local_indices = local_indices[keep_inversion]

                eta_handle = None
                eta_dset = None
                eta_path = os.path.join(
                    self.base_path, f"analysis_eta_{int(km):04d}.h5"
                )
                if error_mode in {"stochastic_r", "generated_r"}:
                    # Rank zero writes the compact perturbation factor once.
                    # Drawing active rows in canonical order reproduces the
                    # mode-1 global RNG call exactly and is chunk invariant.
                    if rank == 0:
                        with h5py.File(obs_file, "r") as source:
                            active_all = (
                                np.asarray(source["obs_active"][:, int(km)], dtype=bool)
                                if compact and "obs_active" in source
                                else np.ones(m, dtype=bool)
                            )
                            active_all &= inversion_keep_all
                            sigma_all = (
                                np.asarray(source["obs_std"][:], dtype=float)
                                if "obs_std" in source
                                else None
                            )
                        if error_mode == "stochastic_r" and sigma_all is None:
                            raise ValueError(
                                "stochastic_R in execution_mode=2 requires "
                                "compact obs_std metadata"
                            )
                        try:
                            os.remove(eta_path)
                        except FileNotFoundError:
                            pass
                        active_positions = np.flatnonzero(active_all)
                        with h5py.File(eta_path, "w") as eta_file:
                            eta_out = eta_file.create_dataset(
                                "eta", (m, Nens), dtype="f8",
                                chunks=(max(1, min(1024, m)), Nens),
                                fillvalue=0.0,
                            )
                            if error_mode == "stochastic_r":
                                seed = int(icesee_kwargs.get("base_seed", 42)) + 1000003 * (
                                    int(km) + 1
                                )
                                rng = np.random.default_rng(seed)
                                draw_rows = max(1, min(
                                    int(icesee_kwargs.get(
                                        "observation_perturbation_row_chunk_size", 1024
                                    )), active_positions.size or 1
                                ))
                                for draw0 in range(0, active_positions.size, draw_rows):
                                    positions = active_positions[draw0:draw0 + draw_rows]
                                    eta = rng.standard_normal((positions.size, Nens))
                                    eta *= sigma_all[positions, None]
                                    eta -= eta.mean(axis=1, keepdims=True)
                                    eta_out[positions, :] = eta
                                eta_file.attrs["seed"] = seed
                            else:
                                selected_rows = obs_indices[active_positions]
                                for member, column in iter_generated_observation_columns(
                                    icesee_kwargs, selected_rows, Nens, self.nd, km
                                ):
                                    eta_out[active_positions, member] = column
                                center_rows = max(1, min(
                                    int(icesee_kwargs.get(
                                        "observation_perturbation_row_chunk_size", 1024
                                    )), active_positions.size or 1
                                ))
                                for row0 in range(0, active_positions.size, center_rows):
                                    positions = active_positions[row0:row0 + center_rows]
                                    values = np.asarray(eta_out[positions, :])
                                    values -= values.mean(axis=1, keepdims=True)
                                    eta_out[positions, :] = values
                    self.mpi_comm.Barrier()
                    eta_handle = h5py.File(eta_path, "r")
                    eta_dset = eta_handle["eta"]

                local_terms = [] if icesee_kwargs.get("local_analysis", False) else None
                try:
                    for row0 in range(0, local_indices.size, max(1, chunk_rows)):
                        rows = local_indices[row0:row0 + chunk_rows]
                        positions = local_positions[row0:row0 + chunk_rows]
                        forecast = _read_rows_preserve_order(state_dset, rows)
                        if compact:
                            observations = _read_rows_preserve_order(
                                obs_dset, positions, int(km)
                            ).reshape(-1)
                        else:
                            observations = _read_rows_preserve_order(
                                obs_dset, rows, int(km)
                            ).reshape(-1)
                        yprime = forecast - forecast.mean(axis=1, keepdims=True)
                        if error_mode == "legacy_prior_anomalies":
                            eta = yprime.copy()
                            innovations = observations[:, None] - forecast
                        else:
                            eta = _read_rows_preserve_order(eta_dset, positions)
                            innovations = observations[:, None] + eta - forecast
                        factor = yprime + eta
                        cross_local += yprime.T @ factor
                        gram_local += factor.T @ factor
                        rhs_local += factor.T @ innovations
                        if local_terms is not None:
                            local_terms.append((
                                rows.copy(), yprime.copy(), eta.copy(),
                                innovations.copy()
                            ))
                finally:
                    if eta_handle is not None:
                        eta_handle.close()
            icesee_kwargs["time_analysis_file_writing"] += MPI.Wtime() - io_start

            cross = np.empty_like(cross_local)
            gram = np.empty_like(gram_local)
            rhs = np.empty_like(rhs_local)
            self.mpi_comm.Allreduce(cross_local, cross, op=MPI.SUM)
            self.mpi_comm.Allreduce(gram_local, gram, op=MPI.SUM)
            self.mpi_comm.Allreduce(rhs_local, rhs, op=MPI.SUM)
            return cross, gram, rhs, local_terms, icesee_kwargs
        except Exception as e:
            print(f"Error in compute_X5_utils: {e}")
            tb_str = "".join(traceback.format_exception(*sys.exc_info()))
            print(f"Traceback details:\n{tb_str}")
            self.mpi_comm.Abort(1)

    def compute_X5_modified(self, **icesee_kwargs):
        # Eta = HA-Hmean where HA = H*state and Hmean = H*mean(state)
        # Dprime[:ens_idx] = d - Hmean
        try:
            cross, gram, rhs, local_terms, icesee_kwargs = self.compute_X5_utils_(
                **icesee_kwargs
            )
            X5 = ensemble_transform_from_products(
                cross, gram, rhs,
                float(icesee_kwargs.get("analysis_svd_energy", 0.999)),
            )
            if not np.allclose(np.sum(X5, axis=0), 1.0, rtol=1e-8, atol=1e-8):
                print(f"\n[ICESEE] X5 column sums differ from 1: {np.sum(X5, axis=0)}\n")

            local_patches = None
            if icesee_kwargs.get("local_analysis", False):
                gathered = self.mpi_comm.gather(local_terms, root=0)
                if self.rank == 0:
                    flattened = [item for rank_terms in gathered for item in rank_terms]
                    local_rows = sum(item[0].size for item in flattened)
                    required = local_rows * self.nens * 8 * 3
                    budget = float(icesee_kwargs.get(
                        "local_analysis_memory_budget_mb", 2048.0
                    )) * 1024.0 * 1024.0
                    if required > budget:
                        raise MemoryError(
                            "Exact grouped local analysis needs its active "
                            "observation factors on rank zero; estimated "
                            f"working set {required / 2**30:.2f} GiB exceeds "
                            f"local_analysis_memory_budget_mb="
                            f"{budget / 2**20:.0f}."
                        )
                    obs_rows = np.concatenate([item[0] for item in flattened])
                    yprime = np.vstack([item[1] for item in flattened])
                    eta = np.vstack([item[2] for item in flattened])
                    innovations = np.vstack([item[3] for item in flattened])
                    from ICESEE.src.utils.localization import compute_local_patches_X5
                    local_patches = compute_local_patches_X5(
                        vec_inputs=icesee_kwargs.get("vec_inputs", []),
                        hdim=self.nd // icesee_kwargs["total_state_param_vars"],
                        HAprime=yprime, Eta=eta,
                        Dprime=innovations, Nens=self.nens,
                        obs_indices=obs_rows,
                        icesee_kwargs=icesee_kwargs,
                    )
                local_patches = self.mpi_comm.bcast(local_patches, root=0)

            return X5, local_patches, icesee_kwargs

        except Exception as e:
            print(f"Error in compute_X5_modified: {e}")
            tb_str = "".join(traceback.format_exception(*sys.exc_info()))
            print(f"Traceback details:\n{tb_str}")
            self.mpi_comm.Abort(1)


    # compute analysis mean
    def compute_analysis_update(self, **icesee_kwargs):
        # Compute the analysis update for each rank
        k = icesee_kwargs.get('k')
        k = k + 1 if k < self.nt - 1 else k
        nt = self.nt
        frozen_was_configured = "frozen_analysis_vars" in icesee_kwargs
        configured_frozen = copy.deepcopy(
            icesee_kwargs.get("frozen_analysis_vars", None)
        )
        try:
            # Mode 1 removes friction from the EnKF state at hybrid cycles.
            # Keep mode 2's fixed file shape, but freeze the same rows after
            # deriving the transform from non-friction observations.
            if bool(icesee_kwargs.get("inversion_flag", False)):
                inversion_info = _inversion_variable_info(self.nd, icesee_kwargs)
                frozen = list(icesee_kwargs.get("frozen_analysis_vars", []) or [])
                if inversion_info["name"] not in frozen:
                    frozen.append(inversion_info["name"])
                icesee_kwargs["frozen_analysis_vars"] = frozen
            self._ensure_batch(k)
            comm = self.mpi_comm
            rank = comm.Get_rank()
            size = comm.Get_size()

            batch_idx = k - self.current_batch_start

            start = MPI.Wtime()

            # call the compute X5 function
            # X5 = self.compute_X5_(k, **icesee_kwargs)
            X5, local_patches, icesee_kwargs = self.compute_X5_modified(
                **icesee_kwargs
            )
            # compute column sums for X5

            local_dim = self.nd_end_world - self.nd_start_world
            Nens = self.nens
            from ICESEE.src.parallelization._mpi_analysis_functions import (
                apply_analysis_controls_local,
            )
            from ICESEE.src.utils.localization import apply_local_patches
            THICKNESS_VARS = {
                "h", "ice_thickness", "thickness", "ice_thick",
                "hi", "h_ice", "h_ice_thickness", "H"
            }
            min_thickness = 1
            variable_map = _global_variable_index_map(
                self.nd, icesee_kwargs.get("vec_inputs", []),
                icesee_kwargs.get("var_nd"),
            )
            thickness_ranges = []
            for name, indices in variable_map.items():
                if name in THICKNESS_VARS and indices.size:
                    thickness_ranges.append((int(indices[0]), int(indices[-1]) + 1))

            row_chunk = resolve_row_chunk_size(
                local_dim, Nens,
                icesee_kwargs.get("analysis_memory_budget_mb", 256),
                icesee_kwargs.get("analysis_row_chunk_size", 0),
                working_buffers=3,
            )
            io_start = MPI.Wtime()
            state_dset = self.datasets[batch_idx]
            temp_path, temp_handle, temp_dset = self._create_analysis_file(k)
            try:
                for local0 in range(0, local_dim, max(1, row_chunk)):
                    global0 = self.nd_start_world + local0
                    global1 = min(self.nd_end_world, global0 + row_chunk)
                    states = np.asarray(state_dset[global0:global1, :])
                    raw_analysis = states @ X5
                    raw_analysis = apply_local_patches(
                        raw_analysis, states,
                        np.arange(global0, global1, dtype=np.int64),
                        local_patches,
                    )
                    updated = apply_analysis_controls_local(
                        analysis_vec=raw_analysis,
                        forecast_vec=states,
                        global_rows=np.arange(global0, global1, dtype=np.int64),
                        icesee_kwargs=icesee_kwargs,
                    )
                    # Generic mode-1 analysis does not impose a thickness
                    # floor here; model-specific geometry projection happens
                    # in the shared finalization stage. Retain the former
                    # mode-2 safeguard only as an explicit opt-in.
                    if icesee_kwargs.get("analysis_project_min_thickness", False):
                        for thick0, thick1 in thickness_ranges:
                            clip0, clip1 = max(global0, thick0), min(global1, thick1)
                            if clip0 < clip1:
                                updated[clip0-global0:clip1-global0, :] = np.maximum(
                                    updated[clip0-global0:clip1-global0, :], min_thickness
                                )
                    temp_dset[global0:global1, :] = updated
                temp_handle.flush()
            finally:
                temp_handle.close()

            self._publish_analysis_file(k, temp_path, icesee_kwargs)

            # Means are derived only after the finalized analysis is published;
            # otherwise diagnostics and restart state can describe different
            # ensembles when geometry/inference modifies the raw transform.
            if bool(icesee_kwargs.get("compute_analysis_mean", False)):
                self._ensure_batch(k)
                published_dset = self.datasets[k - self.current_batch_start]
                file_path = f"{self.base_path}/{self.file_prefix}_mean.h5"
                with h5py.File(file_path, "a", driver="mpio", comm=comm) as mean_handle:
                    if "mean" not in mean_handle:
                        mean_handle.create_dataset(
                            "mean", (self.nd, self.nt),
                            chunks=(min(self.nd, 1000), 1), dtype="f8"
                        )
                    mean_handle["mean"][self.nd_start_world:self.nd_end_world, k] = (
                        np.asarray(
                            published_dset[
                                self.nd_start_world:self.nd_end_world, :
                            ]
                        ).mean(axis=1)
                    )
                    mean_handle.flush()
            elapsed_io = MPI.Wtime() - io_start
            icesee_kwargs["time_analysis_file_writing"] += elapsed_io
            icesee_kwargs["time_analysis_ensemble_mean_generation"] += 0.0
            # print(f"\n[ICESEE] Rank {rank} completed analysis update for time step {k+1}/{nt} analysis_mean norm {np.linalg.norm(analysis_mean)}\n")

            # clean up zarr file
            # if os.path.exists(zarr_path):
            #     shutil.rmtree(zarr_path)
            # for path in [allstates_sate_zarr_path, mean_params_zarr_path, pertubations_zarr_path, analysis_updates_zarr_path]:
            #     if os.path.exists(path):
            #         shutil.rmtree(path)

            self.mpi_comm.Barrier()
            # del all_states_zarr
            gc.collect()
            # ``coefficient`` is frozen only for the hybrid cycle whose EnKF
            # update is followed by inversion.  Do not leak that temporary
            # control into later ordinary EnKF cycles.
            if frozen_was_configured:
                icesee_kwargs["frozen_analysis_vars"] = configured_frozen
            else:
                icesee_kwargs.pop("frozen_analysis_vars", None)
            return icesee_kwargs
            # ----**** --------------------------------------------------------

        except Exception as e:
            print(f"Error in compute_analysis_mean: {e}")
            tb_str = "".join(traceback.format_exception(*sys.exc_info()))
            print(f"Traceback details:\n{tb_str}")
            self.mpi_comm.Abort(1)

    def close(self):
        try:
            self._close_batch()
        except Exception as e:
            print(f"Error occurred in close: {e}")
            tb_str = "".join(traceback.format_exception(*sys.exc_info()))
            print(f"Traceback details:\n{tb_str}")
            self.mpi_comm.Abort(1)
