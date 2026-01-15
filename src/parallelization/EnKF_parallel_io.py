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
import os
import re
import glob
import gc
import zarr
import traceback
import sys
import shutil
from numcodecs import blosc
import time
import functools
from concurrent.futures import ThreadPoolExecutor, as_completed

blosc.use_threads = False

from typing import Callable, Optional, TypeVar, Any
from ICESEE.src.utils.tools import icesee_get_index

T = TypeVar("T")

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
        def wrapper(*args, **kwargs) -> T:
            rank = mpi_comm.Get_rank() if mpi_comm is not None else "N/A"
            for attempt in range(max_attempts):
                try:
                    return func(*args, **kwargs)
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
    def __init__(self, file_prefix, nd, nens, nt, subcomm, mpi_comm, params, \
                 serial_file_creation=False, base_path="enkf_data", batch_size=50,\
                 h5_file_compression=None, h5_file_compression_level=4, h5_file_chunk_size=1000):
        try:
            self.nd = nd
            self.nens = nens
            self.nt = nt
            self.params = params
            self.base_path = base_path
            self.file_prefix = file_prefix
            self.batch_size = batch_size
            self.mpi_comm = mpi_comm
            self.comm = subcomm if nens >= mpi_comm.Get_size() else mpi_comm
            self.rank = self.comm.Get_rank() if nens >= mpi_comm.Get_size() else mpi_comm.Get_rank()
            self.size = self.comm.Get_size() if nens >= mpi_comm.Get_size() else mpi_comm.Get_size()
            self.subcomm = subcomm
            self.serial_file_creation = serial_file_creation
            self.h5_file_compression = h5_file_compression
            self.h5_file_compression_level = h5_file_compression_level
            self.h5_file_chunk_size = h5_file_chunk_size

            # collective I/O threshold (32 is a reasonable default for many systems)
            self.collective_threshold = int(params.get("collective_threshold", 16))
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

            # Create directory and clean up old files
            if mpi_comm.Get_rank() == 0:
                os.makedirs(base_path, exist_ok=True)
                patterns = [f"{self.base_path}/{self.file_prefix}_*.h5"]
                for pattern in patterns:
                    for file_path in glob.glob(pattern):
                        try:
                            os.remove(file_path)
                        except OSError as e:
                            print(f"Error deleting {file_path}: {e}")
            self.mpi_comm.Barrier()

            # Initialize file and dataset lists
            self.files = []
            self.datasets = []
            self.current_batch_start = -1

            # Create initial batch
            # self._create_batch(0)
            if self.serial_file_creation:
                self._create_batch_serial(0)
            else:
                self._create_batch_parallel(0)
        except Exception as e:
            print(f"Error occurred in __init__: {e}")
            tb_str = "".join(traceback.format_exception(*sys.exc_info()))
            print(f"Traceback details:\n{tb_str}")
            self.mpi_comm.Abort(1)

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
            start = MPI.Wtime()
            self._close_batch()
            self.files = []
            self.datasets = []
            self.current_batch_start = t_start
            # nfiles = min(self.batch_size, self.nt - t_start)
            # --- PATCH: create at least one file and clamp to nt ---
            remaining = max(0, self.nt - t_start)
            if remaining == 0:
                print(f"[Rank {self.rank}] WARNING: no files to create (t_start={t_start}, nt={self.nt})")
                return
            
            nfiles = min(self.batch_size, remaining)
            t_end = t_start + nfiles
            print(f"[Rank {self.rank}] Creating batch from {t_start} to {t_end - 1} ({nfiles} files, nt={self.nt})")

            if MPI.COMM_WORLD.Get_rank() == 0:
                for t in range(t_start, t_start + nfiles):
                    fname = f"{self.base_path}/{self.file_prefix}_ens_{t:04d}.h5"
                    with h5py.File(fname, 'w') as f:
                        # row_chunk = min(1024, self.nd)
                        row_chunk = self.nd_local_world
                        # col_chunk = min(32, self.nens)
                        col_chunk = 1
                        f.create_dataset(
                            'states', (self.nd, self.nens),
                            chunks=(row_chunk, col_chunk),
                            # compression="gzip", compression_opts=4,
                            # compression="lzf",
                            compression=None,
                            dtype='f8'
                        )
                        
            self.mpi_comm.Barrier()

            for t in range(t_start, t_start + nfiles):
                fname = f"{self.base_path}/{self.file_prefix}_ens_{t:04d}.h5"
                f = h5py.File(fname, 'a', driver='mpio', comm=self.comm)
                # f = h5py.File(fname, 'a', driver='mpio', comm=self.mpi_comm)
                f.atomic = False
                dset = f['states']
                self.files.append(f)
                self.datasets.append(dset)
        except Exception as e:
            print(f"Error occurred in _create_batch_serial: {e}")
            tb_str = "".join(traceback.format_exception(*sys.exc_info()))
            print(f"Traceback details:\n{tb_str}")
            self.mpi_comm.Abort(1)

    def _create_batch_parallel(self, t_start):
        try:
            start = MPI.Wtime()
            self._close_batch()
            self.files = []
            self.datasets = []
            self.current_batch_start = t_start
            nfiles = min(self.batch_size, self.nt - t_start)
            for t in range(t_start, t_start + nfiles):
                fname = f"{self.base_path}/{self.file_prefix}_{t:04d}.h5"
                f = h5py.File(fname, 'w', driver='mpio', comm=self.comm)
                # f = h5py.File(fname, 'w', driver='mpio', comm=self.mpi_comm)
                f.atomic = False
                # row_chunk = min(1024, self.nd)
                row_chunk = self.nd_local_world
                # col_chunk = min(32, self.nens)
                col_chunk = 1
                dset = f.create_dataset(
                    'states', (self.nd, self.nens),
                    chunks=(row_chunk, col_chunk),
                    # compression="gzip", compression_opts=4,
                    # compression="lzf",
                    compression=None,
                    dtype='f8'
                )
                self.files.append(f)
                self.datasets.append(dset)
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

    @retry_on_failure(max_attempts=3, delay=0.5, mpi_comm=MPI.COMM_WORLD)  # Reduce retries/delays for efficiency
    def _create_batch(self, t_start):
        self._close_batch()
        self.files = []
        self.datasets = []
        self.current_batch_start = t_start
        nfiles = min(self.batch_size, self.nt - t_start)

        # All ranks collectively create files and datasets
        for t in range(t_start, t_start + nfiles):
            fname = f"{self.base_path}/{self.file_prefix}_{t:04d}.h5"
            f = h5py.File(fname, 'w', driver='mpio', comm=self.mpi_comm)
            f.atomic = True  # Enable atomic writes for consistency
            row_chunk = min(1024, self.nd_local)  # Align with local partition
            col_chunk = min(32, self.nens)  # Chunk ensembles for better access
            dset = f.create_dataset(
                'states', (self.nd, self.nens),
                chunks=(row_chunk, col_chunk),
                compression=None,  # Disable for now; test blosc if space is needed
                dtype='f8'
            )
            self.files.append(f)
            self.datasets.append(dset)
        self.mpi_comm.Barrier()  # Single barrier after all creations

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
            batch_start = (t // self.batch_size) * self.batch_size
            if batch_start != self.current_batch_start:
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
            out = np.empty(mem_shape, dtype='f8', order='C')
            ds.id.read(mem_space, file_space, out, dxpl=dxpl)
            return out

    @retry_on_failure(max_attempts=3, delay=0.5, mpi_comm=MPI.COMM_WORLD) 
    def read_forecast(self, t, ens_idx):
        self._ensure_batch(t)
        batch_idx = t - self.current_batch_start
        start = MPI.Wtime()
        data = self.datasets[batch_idx][self.nd_start:self.nd_end, ens_idx]
        # data = self._rw_select(self.datasets[batch_idx],
        #                self.nd_start, self.nd_local, ens_idx, 1, write=False).reshape(-1)

        return data
       

    @retry_on_failure(max_attempts=3, delay=0.5, mpi_comm=MPI.COMM_WORLD) 
    def write_forecast(self, t, data, ens_idx):
        self._ensure_batch(t)
        batch_idx = t - self.current_batch_start
        ds = self.datasets[batch_idx]
        ds[self.nd_start:self.nd_end, ens_idx] = data
        # self._rw_select(self.datasets[batch_idx],
        #         self.nd_start, self.nd_local, ens_idx, 1, buf=data, write=True)

    @retry_on_failure(max_attempts=3, delay=0.5, mpi_comm=MPI.COMM_WORLD) 
    def read_analysis(self, t, ens_idx):
        # try:
        self._ensure_batch(t)
        batch_idx = t - self.current_batch_start
        start = MPI.Wtime()
        data = self.datasets[batch_idx][self.nd_start_world:self.nd_end_world, ens_idx]
        # data = self._rw_select(self.datasets[batch_idx],
        #                self.nd_start_world, self.nd_local_world, ens_idx, 1, write=False).reshape(-1)
        # print(f"[ICESEE] Finished reading analysisensemble {ens_idx} ensemble shape: {data.shape} norm {np.linalg.norm(data)}")
        read_time = MPI.Wtime() - start
        return data
       

    @retry_on_failure(max_attempts=3, delay=0.5, mpi_comm=MPI.COMM_WORLD) 
    def write_analysis(self, t, data, ens_idx):
        # try:
        self._ensure_batch(t)
        batch_idx = t - self.current_batch_start
        start = MPI.Wtime()
        # self.datasets[batch_idx][self.nd_start:self.nd_end, ens_idx] = data
        self.datasets[batch_idx][self.nd_start_world:self.nd_end_world, ens_idx] = data
        # self._rw_select(self.datasets[batch_idx],
        #         self.nd_start_world, self.nd_local_world, ens_idx, 1, buf=data, write=True)
    
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
                        shutil.rmtree(file_path)
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
                                    chunks=(chunk_rows, 1), dtype='f8')
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
        try:
            comm = self.mpi_comm
            rank = comm.Get_rank()

            mean_original = None
            if rank == 0:
                file_path = f"{self.base_path}/{self.file_prefix}_mean.h5"
                self._ensure_batch(t)
                batch_idx = t - self.current_batch_start
                ens_mean = self.datasets[batch_idx][:,:]
                mean_original = np.mean(ens_mean, axis=1)

            self.compute_forecast_mean_chunked(t, ens_chunk_size=ens_chunk_size)
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
    def generate_observation_schedule(self, **kwargs):
        import numpy as np

        t = np.asarray(kwargs["t"], dtype=float)
        if t.ndim != 1 or t.size == 0:
            raise ValueError("`t` must be a 1D non-empty array of times.")
        t_min, t_max = float(t[0]), float(t[-1])

        freq_obs = float(self.params["freq_obs"])
        obs_start = float(self.params["obs_start_time"])
        obs_max_cfg = float(self.params["obs_max_time"])

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
    def _create_synthetic_observations(self, **kwargs):
        import os
        import re
        import h5py
        import numpy as np
        import traceback
        import sys

        nd = self.nd
        nt = self.nt

        obs_t, ind_m, m_obs = self.generate_observation_schedule(**kwargs)
        ind_m = np.asarray(ind_m, dtype=int)  # your code treats these as step indices
        obs_t = np.asarray(obs_t, dtype=float)

        rank = self.mpi_comm.Get_rank()

        total_state_param_vars = self.params["total_state_param_vars"]
        hdim = nd // total_state_param_vars

        # ----------------------------
        # Bed snapshots: interpret as YEARS like partial-parallel
        # ----------------------------
        bed_snaps = np.asarray(kwargs.get("bed_obs_snapshot", []), dtype=float)
        bed_snap_cols = []
        bed_time_to_col = {}

        if bed_snaps.size > 0 and obs_t.size > 0:
            for bed_time in bed_snaps:
                diffs = np.abs(obs_t - bed_time)
                j = int(np.argmin(diffs))  # 0-based obs column
                bed_snap_cols.append(j)
                bed_time_to_col[bed_time] = j
            bed_snap_cols = sorted(set(bed_snap_cols))

        if rank == 0:
            try:
                print("[ICESEE] Generating synthetic observations ...")
                print(f"[ICESEE] observation times: {obs_t}, indices: {ind_m}, total: {m_obs}")
                print(f"[ICESEE] bed_snaps (years): {bed_snaps}")
                print(f"[ICESEE] bed_snap_cols (obs columns): {bed_snap_cols}")
                if len(bed_snap_cols) > 0:
                    print(f"[ICESEE] obs_t at bed_snap_cols: {obs_t[bed_snap_cols]}")

                obs_file = f"{self.base_path}/synthetic_obs.h5"

                # Load true state
                with h5py.File(f"{self.base_path}/true_nurged_states.h5", "r") as f:
                    statevec_true = f["true_state"][:]  # (nd, nt or nt+1)

                # Build index map once
                _, indx_map, _ = icesee_get_index(**kwargs)
                vec_inputs = list(kwargs["vec_inputs"])

                # Helpers / options
                num_state_vars = kwargs.get("num_state_vars", self.params.get("num_state_vars"))
                observed_params = set(kwargs.get("observed_params", []))

                bed_aliases = {"bed", "bedrock", "bed_topography", "bedtopo", "bedtopography"}
                key_is_bed = {k: (k in bed_aliases) for k in vec_inputs}
                key_idx_map = {k: np.asarray(indx_map[k], dtype=int) for k in vec_inputs}

                # ---- Bed sparsity controls (accept BOTH naming conventions safely) ----
                Ly = kwargs.get("Ly", self.params.get("Ly", None))
                Lx = kwargs.get("Lx", self.params.get("Lx", None))
                model_name = kwargs.get("model_name", None)

                # allow either key name:
                bed_stride_km = kwargs.get("bed_obs_stride", None)
                if bed_stride_km is None:
                    bed_stride_km = kwargs.get("bed_obs_stride_km", None)

                bed_spacing_pts = kwargs.get("bed_obs_spacing", None)
                bed_indices_user = kwargs.get("bed_obs_indices", None)
                bed_mask_user = kwargs.get("bed_obs_mask", None)

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

                            icesee_path = kwargs.get("icesee_path")
                            data_path = kwargs.get("data_path")

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

                # ---- Create / overwrite output datasets ----
                with h5py.File(obs_file, "a") as f:
                    if "hu_obs" in f or "error_R" in f:
                        print(f"[ICESEE] Warning: {obs_file} already contains 'hu_obs' or 'error_R'. Overwriting datasets.")
                        if "hu_obs" in f:
                            del f["hu_obs"]
                        if "error_R" in f:
                            del f["error_R"]

                    hu_obs = f.create_dataset(
                        "hu_obs", (nd, m_obs),
                        chunks=(min(1000, nd), min(50, m_obs)),
                        dtype="f8"
                    )
                    error_R = f.create_dataset(
                        "error_R", (nd, m_obs * 2 + 1),
                        chunks=(min(1000, nd), min(50, m_obs * 2 + 1)),
                        dtype="f8"
                    )

                    # Fill error_R by blocks (same idea)
                    sig_obs = np.asarray(self.params["sig_obs"]).reshape(-1)
                    if sig_obs.size < total_state_param_vars:
                        sig_obs = np.pad(sig_obs, (0, total_state_param_vars - sig_obs.size), mode="edge")
                    elif sig_obs.size > total_state_param_vars:
                        sig_obs = sig_obs[:total_state_param_vars]

                    for i, sig in enumerate(sig_obs):
                        s = i * hdim
                        e = s + hdim
                        error_R[s:e, :] = sig

                    # ==== Main observation loop (columns km), aligned to partial ====
                    obs_set = set(kwargs.get("observed_vars", [])) | set(kwargs.get("observed_params", []))

                    km = 0
                    for step in range(nt):
                        if km >= m_obs:
                            break

                        if step == ind_m[km]:
                            # guard for nt vs nt+1 storage
                            tcol = step + 1 if (step + 1) < statevec_true.shape[1] else step

                            for ii, key in enumerate(vec_inputs):
                                idx = key_idx_map[key]
                                bed_flag = key_is_bed.get(key, False)

                                # ---- non-bed vars: normal obs at every obs time ----
                                if (key in obs_set) and (not bed_flag):
                                    sigma = error_R[idx, km]
                                    hu_obs[idx, km] = statevec_true[idx, tcol] + \
                                                    np.random.normal(0.0, sigma, size=idx.size)
                                else:
                                    hu_obs[idx, km] = 0.0

                                # ---- bed vars: only at bed snapshot columns, on masked subset ----
                                if bed_flag and (km in bed_snap_cols):
                                    mask = bed_mask_map.get(key, np.ones(idx.size, dtype=bool))
                                    idx_obs = idx[mask]
                                    if idx_obs.size > 0:
                                        sigma_obs = error_R[idx_obs, km]
                                        hu_obs[idx_obs, km] = statevec_true[idx_obs, tcol] + \
                                                            np.random.normal(0.0, sigma_obs, size=idx_obs.size)

                            km += 1

            except Exception as e:
                print(f"[ICESEE] Error in full-parallel _create_synthetic_observations: {e}")
                tb_str = "".join(traceback.format_exception(*sys.exc_info()))
                print(f"Traceback details:\n{tb_str}")
                raise

        self.mpi_comm.Barrier()
        return ind_m, m_obs

    def H_matrix(self, **kwargs):
        """
        Fully-parallel version: write H directly to Zarr, but use the SAME
        observation-index logic as the partial-parallel run (including bed masks).

        H contains ONLY real observation points (and bed subsampling/masking),
        with rows ordered exactly as `params["all_observed"]` concatenation.
        """
        try:
            zarr_path = kwargs.get("H_matrix_zarr_path")
            if zarr_path is None:
                raise ValueError("H_matrix_zarr_path is required in kwargs")

            params = self.params
            observed = params["all_observed"]  # e.g. ['h','u','v','smb','bed']
            vec_inputs = kwargs.get("vec_inputs", [])

            # --- Recompute index map (same as partial parallel) ---
            vecs, indx_map, _ = icesee_get_index(**kwargs)

            nd = int(self.nd)  # state size

            # --- Retrieve bed masks (must already exist) ---
            bed_mask_map = kwargs.get("bed_mask_map", {})

            bed_aliases = {"bed", "bedrock", "bed_topography", "bedtopo", "bedtopography"}
            key_is_bed = {k: (k in bed_aliases) for k in vec_inputs}

            # --- Collect observation indices, applying masks to bed ---
            all_obs_indices = []
            for key in observed:
                if key not in indx_map:
                    raise KeyError(f"Observed key '{key}' not found in indx_map")

                idx = np.asarray(indx_map[key], dtype=int)

                # Apply bed mask if present (keep behavior consistent with your partial-parallel code)
                if key == "bed" and key in bed_mask_map:
                    mask = np.asarray(bed_mask_map[key], dtype=bool)

                    # tolerate mask stored as shape (1, n) like your partial version
                    if mask.ndim > 1 and mask.shape[0] == 1:
                        mask = mask[0]
                    else:
                        mask = mask.ravel()

                    if mask.size != idx.size:
                        raise ValueError(
                            f"bed_mask_map['{key}'] length {mask.size} does not "
                            f"match bed vector length {idx.size}"
                        )

                    idx = idx[mask]

                all_obs_indices.append(idx)

            # Flatten
            obs_indices = np.concatenate(all_obs_indices).astype(int)

            if obs_indices.size == 0:
                raise ValueError("No observation indices found (obs_indices is empty).")

            # --- Safety ---
            if obs_indices.max() >= nd:
                raise ValueError(
                    f"H_matrix error: obs index {obs_indices.max()} >= state size {nd}"
                )

            # m = number of real observations (rows)
            m = int(obs_indices.size)

            # --- Allocate Zarr H ---
            # Use chunking  with correct row count.
            H_matrix_file = zarr.create_array(
                store=zarr_path,
                shape=(m, nd),
                chunks=(min(50, m), min(1000, nd)),
                dtype="f8",
                overwrite=True,
            )

            # --- Write identity rows into Zarr in blocks (memory-safe) ---
            block = int(kwargs.get("H_write_block_rows", 5000))
            for r0 in range(0, m, block):
                r1 = min(m, r0 + block)
                rows = np.arange(r0, r1, dtype=int)
                cols = obs_indices[r0:r1]
            #     # Set the 1's for this block; everything else stays 0
            #     H_matrix_file[rows, cols] = 1.0
                for rr, cc in zip(rows, cols):
                    H_matrix_file[rr, cc] = 1.0
                

            # --- Joint estimation: zero out parameter columns (preserve your original behavior) ---
            if params.get("joint_estimation", False):
                ndim = nd // params["total_state_param_vars"]
                state_variables_size = ndim * params["num_state_vars"]
                # if any obs accidentally landed in param part, this forces them off
                H_matrix_file[:, state_variables_size:] = 0.0

        except Exception as e:
            print(f"Error in H_matrix: {e}")
            tb_str = "".join(traceback.format_exception(*sys.exc_info()))
            print(f"Traceback details:\n{tb_str}")
            self.mpi_comm.Abort(1)

    def compute_X5_utils_(self, **kwargs):
        # Eta = HA-Hmean where HA = H*state and Hmean = H*mean(state)
        # Dprime[:ens_idx] = d - Hmean
        k = kwargs.get('k')
        k = k + 1 if k < self.nt - 1 else k
        km = kwargs.get('km')
        try:
            H_matrix_zarr_path = kwargs.get('H_matrix_zarr_path', f"{self.base_path}/H_matrix.zarr")
            synthetic_obs_zarr_path = kwargs.get('synthetic_obs_zarr_path', f"{self.base_path}/synthetic_obs.zarr")
            m = kwargs.get('m_obs')
            Nens = self.nens

            # --- Open H (read-only) once and slice local columns
            H_matrix = zarr.open_array(H_matrix_zarr_path, mode='r')
            # H has shape (m, nd_total); we take our local column block
            H_local = H_matrix[:, self.nd_start_world:self.nd_end_world]  # (m, local_nd)

            local_nd = self.nd_end_world - self.nd_start_world

            # --- Read ensemble mean ONCE (parallel)
            mean_file_path = f"{self.base_path}/{self.file_prefix}_mean.h5"
            with h5py.File(mean_file_path, 'r', driver='mpio', comm=self.mpi_comm) as f:
                ens_mean_local = f['mean'][self.nd_start_world:self.nd_end_world, k ]  # (local_nd,)

            # --- Synthetic obs (once)
            # synthetic_obs = zarr.open_array(synthetic_obs_zarr_path, mode='r')
            # synthetic_obs_local = synthetic_obs[self.nd_start_world:self.nd_end_world, km]  # (local_nd,)

            # *--with open synthetic obs h5file *---rememdy for now---*
            obs_file = f"{self.base_path}/synthetic_obs.h5"
            with h5py.File(obs_file, 'r', driver='mpio', comm=self.mpi_comm) as f:
                synthetic_obs_local = f['hu_obs'][self.nd_start_world:self.nd_end_world, km]  # (local_nd,)
            # *---rememdy for now---*

            # print(f"\n[Rank {self.mpi_comm.Get_rank()}] H_local norm : {np.linalg.norm(H_local)}, ens_mean_local norm: {np.linalg.norm(ens_mean_local)}, synthetic_obs_local norm: {np.linalg.norm(synthetic_obs_local)}\n")

            # --- d = H * y_obs (one GEMV + one Allreduce)
            d_local = H_local @ synthetic_obs_local              # (m,)
            d_global = np.empty_like(d_local)
            self.mpi_comm.Allreduce(d_local, d_global, op=MPI.SUM)  # 1st collective

            # --- Hmean = H * ens_mean (one GEMV + (no extra open calls))
            Hmean_local = H_local @ ens_mean_local               # (m,)
            Hmean_global = np.empty_like(Hmean_local)
            self.mpi_comm.Allreduce(Hmean_local, Hmean_global, op=MPI.SUM)  # 2nd collective

            # --- Read all ensemble states locally and batch into a matrix
            _local_analysis_time = MPI.Wtime()
            # Shape: (local_nd, Nens)
            # States_local = np.empty((local_nd, Nens), dtype=H_local.dtype, order='C')
            States_local =zarr.zeros((local_nd, Nens), dtype=H_local.dtype, store=f"{self.base_path}/States_local_{self.mpi_comm.Get_rank()}.zarr", chunks=(local_nd, 1), overwrite=True)
            for j in range(Nens):
                # States_local[:, j] = self.read_analysis(k, j)  # each returns local slice (local_nd,) 
                States_local[:, j] = self.read_analysis(k, j)  # each returns local slice (local_nd,)
            _local_analysis_time = MPI.Wtime() - _local_analysis_time
            kwargs["time_analysis_file_writing"] += _local_analysis_time


            # --- HA for all ensemble members in one GEMM + one Allreduce on the whole matrix
            # local (m, local_nd) @ (local_nd, Nens) -> (m, Nens)
            HA_local = H_local @ States_local                   # matrix-matrix multiply
            # Single Allreduce over full (m, Nens) buffer
            HA = np.empty_like(HA_local, order='C')
            # Allreduce on a contiguous buffer; flatten for safety
            self.mpi_comm.Allreduce(
                [HA_local, MPI.DOUBLE],
                [HA, MPI.DOUBLE],
                op=MPI.SUM
            )


            # --- Eta = HA - Hmean[:, None]
            Eta = HA - Hmean_global[:, None]             # (m, Nens)

            # --- D' = (d - Hmean)[:, None], same for all ensemble members
            Dprime = (d_global - Hmean_global)[:, None] * np.ones((1, Nens), dtype=HA.dtype)

            # print(f"\n[Rank {self.mpi_comm.Get_rank()}] norms H: {np.linalg.norm(H_matrix)},  ens_mean:{np.linalg.norm(np.mean(States_local, axis=1))}, d: {np.linalg.norm(d_local)} D: {np.linalg.norm(d_global.reshape(-1,1) + Eta)}, HA: {np.linalg.norm(HA)}, Eta: {np.linalg.norm(Eta)}, ensemble_vec: {np.linalg.norm(States_local)} \n")

            return Dprime, Eta, Eta, kwargs
        except Exception as e:
            print(f"Error in compute_X5_utils: {e}")
            tb_str = "".join(traceback.format_exception(*sys.exc_info()))
            print(f"Traceback details:\n{tb_str}")
            self.mpi_comm.Abort(1)

    def compute_X5_modified(self, **kwargs):
        # Eta = HA-Hmean where HA = H*state and Hmean = H*mean(state)
        # Dprime[:ens_idx] = d - Hmean
        try:
            Nens = self.nens

            # Dprime, Eta, HA, kwargs = self.compute_X5_utils_(km, **kwargs)
            Dprime, Eta, HAprime, kwargs  = self.compute_X5_utils_(**kwargs)
            # print(f"\n [Rank {self.mpi_comm.Get_rank()}] Dprime norm: {np.linalg.norm(Dprime)}, Eta norm: {np.linalg.norm(Eta)}, HA norm: {np.linalg.norm(HA)} \n")
            # Dprime, Eta, HA = self.compute_X5_utils_batch(**kwargs)

            # compute the HAbar
            # HAbar = np.mean(HA, axis=1)
            # HAprime = HA - HAbar[:, np.newaxis]
            # one_N = np.ones((Nens,Nens))/Nens
            # HAprime= HA@(np.eye(Nens) - one_N) # mxNens

            # get m
            m = HAprime.shape[0]

            # compute HAprime + Eta
            HAprime_Eta = HAprime + Eta
            # HAprime_Eta = Eta  # since HAprime = Eta
           
            # print(f"\n[Rank {self.mpi_comm.Get_rank()}] HAprime_Eta norm: {np.linalg.norm(HAprime_Eta)}, shape: {HAprime_Eta.shape}\n")
            # print(f"\n [Rank {self.mpi_comm.Get_rank()}] Dprime_local shape: {Dprime_local.shape} HAprime_local shape: {HAprime_local.shape} HAprime_Eta_local shape: {HAprime_Eta_local.shape}\n ")

            # compute SVD of HAprime_Eta
            U, sig, Vt = np.linalg.svd(HAprime_Eta, full_matrices=False)

            # get the min (m Nens)
            nrmin = min(Nens, m)
            
            # convert S to eigenvalues
            sig = sig**2

            sigsum = np.sum(sig[:nrmin])
            # print(f"[Rank {self.mpi_comm.Get_rank()}] sigsum: {sigsum}, sig: {sig[:nrmin]}")
            sigsum1 = 0.0
            nrsigma = 0
            if sigsum == 0:
                print(f"[Rank {self.mpi_comm.Get_rank()}] Warning: sigsum is zero, setting nrsigma to 0")
                nrsigma = 0
                sig[:] = 0.0
            else:
                for i in range(nrmin):
                    if sigsum1 / sigsum < 0.999:
                        nrsigma += 1
                        sigsum1 += sig[i]
                        sig[i] = 1.0 / sig[i]
                    else:
                        sig[i:nrmin] = 0.0
                        break

            # compute X1 = sig*UT (Nens x m)
            X1 = np.empty((nrmin, m))
            for j in range(m):
                for i in range(nrmin):
                    X1[i,j] =sig[i]*U[j,i]

            # compute X2 = X1*Dprime # Nens x Nens
            X2 = np.dot(X1, Dprime)
            # Dprime_mat = np.broadcast_to(Dprime, (m, Nens))
            # X2 = X1 @ Dprime_mat

            #  compute X3 = U*X2 # m_obs x Nens
            X3 = np.dot(U, X2)

            # print(f"[ICESEE] Rank: {rank_world} X3 shape: {X3.shape}")
            # compute X4 = (HAprime.T)*X3 # Nens x Nens
            X4 = np.dot(HAprime.T, X3)
            # del X2, X3, U, HAprime, HA, Eta, Dprime
            del X2, X3, U, HAprime, HAprime_Eta, Eta, Dprime
            gc.collect()

            # compute X5 = X4 + I
            X5 = X4 + np.eye(Nens)
            # sum of each column of X5 should be 1
            if np.sum(X5, axis=0).all() != 1.0:
                print(f"\n[ICESEE] Sum of each X5 column is not 1.0: {np.sum(X5, axis=0)}\n")
            # print(f"[ICESEE] Rank: {self.mpi_comm.Get_rank()} X5 sum: {np.sum(X5, axis=0)}")
            del X4; gc.collect()

            return X5, kwargs

        except Exception as e:
            print(f"Error in compute_X5_modified: {e}")
            tb_str = "".join(traceback.format_exception(*sys.exc_info()))
            print(f"Traceback details:\n{tb_str}")
            self.mpi_comm.Abort(1)


    # compute analysis mean
    def compute_analysis_update(self, **kwargs):
        # Compute the analysis update for each rank
        k = kwargs.get('k')  
        k = k + 1 if k < self.nt - 1 else k
        nt = self.nt
        try:
            self._ensure_batch(k)
            comm = self.mpi_comm
            rank = comm.Get_rank()
            size = comm.Get_size()

            batch_idx = k - self.current_batch_start

            start = MPI.Wtime()

            # call the compute X5 function
            # X5 = self.compute_X5_(k, **kwargs)
            X5, kwargs = self.compute_X5_modified(**kwargs)
            # compute column sums for X5

            # ---
            local_dim = self.nd_end_world - self.nd_start_world
            Nens = self.nens

            # ----works----
            # Read all ensemble data at once
            # all_states = np.zeros((local_dim, Nens))
            # write all_states to zarr file *--------------------------------
            allstates_sate_zarr_path = f"{self.base_path}/all_states_rank_{rank}.zarr"
            mean_params_zarr_path = f"{self.base_path}/mean_params_rank_{rank}.zarr"
            pertubations_zarr_path = f"{self.base_path}/pertubations_rank_{rank}.zarr"
            analysis_updates_zarr_path = f"{self.base_path}/analysis_updates_rank_{rank}.zarr"

            for path in [mean_params_zarr_path, pertubations_zarr_path, analysis_updates_zarr_path]:
                if os.path.exists(path):
                    shutil.rmtree(path)

            # if os.path.exists(allstates_sate_zarr_path):
            #     shutil.rmtree(allstates_sate_zarr_path)
            all_states_zarr = zarr.create_array(store=allstates_sate_zarr_path, shape=(local_dim, Nens), chunks=(min(1000, local_dim), min(50, Nens)), dtype='f8', overwrite=True)
            mean_params = zarr.create_array(store=mean_params_zarr_path, shape=(local_dim, 1), chunks=(min(1000, local_dim), 1), dtype='f8', overwrite=True)
            pertubations = zarr.create_array(store=pertubations_zarr_path, shape=(local_dim, Nens), chunks=(min(1000, local_dim), min(50, Nens)), dtype='f8', overwrite=True)
            analysis_updates = zarr.create_array(store=analysis_updates_zarr_path, shape=(local_dim, Nens), chunks=(min(1000, local_dim), min(50, Nens)), dtype='f8', overwrite=True)

            _local_analysis_time0 = MPI.Wtime()
            for i in range(Nens):
                # all_states_zarr[:, i] = self.read_analysis(k, i)
                all_states_zarr[:, i] = self.read_analysis(k, i)
            _local_analysis_time0 = MPI.Wtime() - _local_analysis_time0

            # Compute analysis updates for all paucall ensembles using matrix multiplication
            analysis_updates = all_states_zarr @ X5  # Matrix multiplication

            # performm inflation
            inflation_factor = self.params.get('inflation_factor', 1.0)
            ndim = analysis_updates.shape[0]//self.params["total_state_param_vars"]
            state_block_size = ndim * self.params["num_state_vars"]
            mean_params = np.mean(analysis_updates[state_block_size:, :], axis=1).reshape(-1, 1)
            pertubations = analysis_updates[state_block_size:, :] - mean_params
            # inflated_pertubations = pertubations * inflation_factor
            analysis_updates[state_block_size:, :] = mean_params + (pertubations * inflation_factor)

            # check for negative thicknes and set to 1e-3 if vec_input contains h
            # Define valid thickness variable names
            THICKNESS_VARS = {
                "h", "ice_thickness", "thickness", "ice_thick", 
                "hi", "h_ice", "h_ice_thickness", "H"
            }
            min_thickness = 1
            vec_inputs = kwargs.get("vec_inputs", None)
         
            for i, input_var in enumerate(vec_inputs or []):
                if input_var in THICKNESS_VARS:
                    start = i * ndim
                    end = start + ndim
                    analysis_updates[start:end, :] = np.maximum(analysis_updates[start:end, :], min_thickness)

            # Write back all analysis updates
            _local_analysis_time1 = MPI.Wtime()
            # _k = k + 1 if k < self.nt - 1 else k
            for j in range(Nens):
                self.write_analysis(k, analysis_updates[:, j], j)
            _local_analysis_time1 = MPI.Wtime() - _local_analysis_time1
            kwargs["time_analysis_file_writing"] += (_local_analysis_time0 + _local_analysis_time1)

            # compute the anlysis mean and write to h5 file
            _time_analysis_mean = MPI.Wtime()
            if kwargs.get('compute_analysis_mean', False):
                yi = np.sum(X5, axis=1)
                analysis_mean = np.dot(all_states_zarr, yi) / Nens
                # print(f"[Rank {self.mpi_comm.Get_rank()}]  shape: {analysis_mean.shape} shape_ {self.nd_end_world - self.nd_start_world}")
                file_path = f"{self.base_path}/{self.file_prefix}_mean.h5"
                with h5py.File(file_path, 'a', driver='mpio', comm=self.mpi_comm) as f:
                    if 'mean' not in f:
                        if rank == 0:
                            f.create_dataset(
                                'mean', (self.nd, self.nt),
                                chunks=(min(self.nd, 1000), 1),
                                dtype='f8'
                            )
                    comm.Barrier()
                    # _k = k + 1 if k < self.nt - 1 else k
                    f['mean'][self.nd_start_world:self.nd_end_world, k] = analysis_mean
            kwargs["time_analysis_ensemble_mean_generation"] += (MPI.Wtime() - _time_analysis_mean)
            # print(f"\n[ICESEE] Rank {rank} completed analysis update for time step {k+1}/{nt} analysis_mean norm {np.linalg.norm(analysis_mean)}\n")

            # clean up zarr file
            # if os.path.exists(zarr_path):
            #     shutil.rmtree(zarr_path)
            for path in [allstates_sate_zarr_path, mean_params_zarr_path, pertubations_zarr_path, analysis_updates_zarr_path]:
                if os.path.exists(path):
                    shutil.rmtree(path)
            
            self.mpi_comm.Barrier()
            # del all_states_zarr
            gc.collect()
            return kwargs
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

    
    # def compute_X5_utils_batch(self, **kwargs):
    #     """
    #     Returns:
    #     Eta   : (m, Nens)
    #     Dprime: (m, Nens)   # constant across Nens (each column equal)
    #     HA    : (m, Nens)
    #     """
    #     k = kwargs.get('k')
    #     k = k + 1 if k < self.nt - 1 else k
    #     km = kwargs.get('km', k)
    #     nt = self.nt
    #     try:
    #         comm = self.mpi_comm
    #         rank = comm.Get_rank()

    #         # ---- Config / inputs
    #         H_matrix_zarr_path = kwargs.get('H_matrix_zarr_path', f"{self.base_path}/H_matrix.zarr")
    #         synthetic_obs_zarr_path = kwargs.get('synthetic_obs_zarr_path', f"{self.base_path}/synthetic_obs.zarr")
    #         m = kwargs.get('m_obs')
    #         Nens = int(self.nens)
    #         block_size = int(kwargs.get('block_size', max(16, min(64, Nens))))  # tuneable batch size

    #         # ---- Open H once and slice local columns
    #         H_matrix = zarr.open_array(H_matrix_zarr_path, mode='r')  # shape (m_total, nd_total) or (m, nd_total) as per your layout
    #         H_local = H_matrix[:, self.nd_start_world:self.nd_end_world]  # shape (m, local_nd)
    #         H_local = np.ascontiguousarray(H_local, dtype=np.float64)
    #         m_local, local_nd = H_local.shape  # expect m_local == m

    #         # ---- Read local ensemble mean once
    #         mean_file_path = f"{self.base_path}/{self.file_prefix}_mean.h5"
    #         with h5py.File(mean_file_path, 'r', driver='mpio', comm=comm) as f:
    #             ens_mean_local = f['mean'][self.nd_start_world:self.nd_end_world, k ]
    #         ens_mean_local = np.ascontiguousarray(ens_mean_local, dtype=np.float64)  # (local_nd,)

    #         # ---- Synthetic observations (local slice)
    #         # synthetic_obs = zarr.open_array(synthetic_obs_zarr_path, mode='r')
    #         # synthetic_obs_local = synthetic_obs[self.nd_start_world:self.nd_end_world, km]
    #         # synthetic_obs_local = np.ascontiguousarray(synthetic_obs_local, dtype=np.float64)  # (local_nd,)
    #          # *--with open synthetic obs h5file *---rememdy for now---*
    #         obs_file = f"{self.base_path}/synthetic_obs.h5"
    #         with h5py.File(obs_file, 'r', driver='mpio', comm=self.mpi_comm) as f:
    #             synthetic_obs_local = f['hu_obs'][self.nd_start_world:self.nd_end_world, km]  # (local_nd,)
    #         # *---rememdy for now---*

    #         # ---- Fuse d and Hmean into a single GEMM + single Allreduce
    #         # Build a 2-column local matrix [y_obs, ens_mean]
    #         V_local = np.empty((local_nd, 2), dtype=np.float64, order='C')
    #         V_local[:, 0] = synthetic_obs_local
    #         V_local[:, 1] = ens_mean_local

    #         Y_local = H_local @ V_local                   # (m, 2)
    #         Y_global = np.empty_like(Y_local, order='C')  # (m, 2)

    #         # Single collective for both vectors
    #         comm.Allreduce([Y_local, MPI.DOUBLE], [Y_global, MPI.DOUBLE], op=MPI.SUM)
    #         d_global     = Y_global[:, 0]                 # (m,)
    #         Hmean_global = Y_global[:, 1]                 # (m,)

    #         # ---- Compute HA for all ensemble members in batches
    #         HA_global = np.empty((m_local, Nens), dtype=np.float64, order='C')

    #         for j0 in range(0, Nens, block_size):
    #             j1 = min(j0 + block_size, Nens)
    #             B = j1 - j0

    #             # Load a contiguous local block of states: shape (local_nd, B)
    #             States_local_blk = np.empty((local_nd, B), dtype=np.float64, order='C')
    #             for jj, ens_idx in enumerate(range(j0, j1)):
    #                 States_local_blk[:, jj] = self.read_analysis(k, ens_idx)  # must return local slice (local_nd,)

    #             # Local GEMM then one Allreduce for this batch
    #             HA_local_blk = H_local @ States_local_blk             # (m, B)
    #             HA_global_blk = np.empty_like(HA_local_blk, order='C')
    #             comm.Allreduce([HA_local_blk, MPI.DOUBLE], [HA_global_blk, MPI.DOUBLE], op=MPI.SUM)

    #             HA_global[:, j0:j1] = HA_global_blk

    #         # ---- Eta and D'
    #         # Eta = HA - Hmean[:, None]
    #         Eta = HA_global - Hmean_global[:, None]                   # (m, Nens)
    #         # print(f"\n[Rank {self.mpi_comm.Get_rank()}] H_local norm : {np.linalg.norm(H_local)}, ens_mean_local norm: {np.linalg.norm(ens_mean_local)}, synthetic_obs_local norm: {np.linalg.norm(synthetic_obs_local)}\n")

    #         # D' = (d - Hmean) broadcast across columns
    #         d_minus_Hmean = (d_global - Hmean_global)                 # (m,)
    #         # Make every column identical, no extra collectives
    #         Dprime = np.broadcast_to(d_minus_Hmean[:, None], (m_local, Nens)).copy()

    #         return Eta, Dprime, HA_global

    #     except Exception as e:
    #         print(f"Error in compute_X5_utils: {e}")
    #         tb_str = "".join(traceback.format_exception(*sys.exc_info()))
    #         print(f"Traceback details:\n{tb_str}")
    #         self.mpi_comm.Abort(1)

    # def compute_X5_utils(self, **kwargs):
    #     """
    #     Parallel EnKF X5 minimal utilities:
    #     d_global   = H @ y_obs
    #     Hbar       = H @ ensemble_mean
    #     Dprime     = d_global - Hbar
    #     HAprime    = Eta (from H @ ensemble_perturbations)
    #     """

    #     import numpy as np, h5py, zarr, sys, traceback
    #     from mpi4py import MPI

    #     comm  = self.mpi_comm
    #     rank  = comm.Get_rank()
    #     size  = comm.Get_size()
    #     k     = kwargs.get('k')
    #     k = k + 1 if k < self.nt - 1 else k
    #     km = kwargs.get('km', k)

    #     try:
    #         # ----------------------------------------------------------------------
    #         # Parameters and file paths
    #         # ----------------------------------------------------------------------
    #         H_matrix_zarr_path = kwargs.get('H_matrix_zarr_path', f"{self.base_path}/H_matrix.zarr")
    #         # synthetic_obs_zarr_path = kwargs.get('synthetic_obs_zarr_path', f"{self.base_path}/synthetic_obs.zarr")
    #         Nens = self.nens

    #         if kwargs.get("inversion_flag", False):
    #             # exlude friction from assimilation
    #             friction_idx = kwargs.get("friction_idx", None)
    #             excluded_indices = [friction_idx]

    #             # entire state vector on each rank
    #             i0, i1 = self.nd_start_world, self.nd_end_world
    #             local_nd = i1 - i0

    #             print(f"[rank {rank}] Excluding friction index {local_nd} from assimilation"); exit(0)

    #         else:
    #             i0, i1 = self.nd_start_world, self.nd_end_world
    #             local_nd = i1 - i0

    #             # Choose MPI datatype dynamically (float32 or float64)
    #             mpi_type = MPI._typedict[np.dtype(np.float64).char]

    #             # ----------------------------------------------------------------------
    #             # Load local column block of H and local slices of ensemble_mean, obs
    #             # ----------------------------------------------------------------------
    #             H_matrix = zarr.open_array(H_matrix_zarr_path, mode='r')
    #             H_local  = np.asarray(H_matrix[:, i0:i1], dtype=np.float64, order='C')  # (m, local_nd)

    #             # Synthetic observation local slice (y)
    #             # synthetic_obs = zarr.open_array(synthetic_obs_zarr_path, mode='r')
    #             # y_local = np.asarray(synthetic_obs[i0:i1, km], dtype=np.float64)  # (local_nd,)
    #             obs_file = f"{self.base_path}/synthetic_obs.h5"
    #             with h5py.File(obs_file, 'r', driver='mpio', comm=self.mpi_comm) as f:
    #                 y_local = np.asarray(f['hu_obs'][i0:i1, km], dtype=np.float64)  # (local_nd,)

    #             # Ensemble mean local slice
    #             mean_file_path = f"{self.base_path}/{self.file_prefix}_mean.h5"
    #             with h5py.File(mean_file_path, 'r', driver='mpio', comm=comm) as f:
    #                 # _k = k + 1 if k < self.nt - 1 else k
    #                 ensemble_mean_local = np.asarray(f['mean'][i0:i1, k], dtype=np.float64)  # (local_nd,)

    #             # ----------------------------------------------------------------------
    #             # Compute d = H * y_obs   (GEMV + Allreduce)
    #             # ----------------------------------------------------------------------
    #             d_local = H_local @ y_local                    # (m,)
    #             d_global = np.empty_like(d_local)
    #             comm.Allreduce([d_local, mpi_type], [d_global, mpi_type], op=MPI.SUM)

    #             # ----------------------------------------------------------------------
    #             # Compute Hbar = H * ensemble_mean   (GEMV + Allreduce)
    #             # ----------------------------------------------------------------------
    #             Hbar_local = H_local @ ensemble_mean_local     # (m,)
    #             Hbar_global = np.empty_like(Hbar_local)
    #             comm.Allreduce([Hbar_local, mpi_type], [Hbar_global, mpi_type], op=MPI.SUM)

    #             # ----------------------------------------------------------------------
    #             # Compute Eta = H * (ensemble_vec - ensemble_mean)
    #             # ----------------------------------------------------------------------
    #             # Each rank reads its local portion of ensemble states
    #             States_local = np.empty((local_nd, Nens), dtype=np.float64, order='C')
    #             for j in range(Nens):
    #                 States_local[:, j] = np.asarray(self.read_analysis(k, j), dtype=np.float64)

    #             # Compute local ensemble perturbations
    #             perturb_local = States_local - ensemble_mean_local[:, None]  # (local_nd, Nens)

    #             # Project ensemble perturbations into observation space
    #             Eta_local = H_local @ perturb_local                          # (m, Nens)
    #             Eta = Eta_local.copy(order='C')
    #             comm.Allreduce(MPI.IN_PLACE, [Eta, mpi_type], op=MPI.SUM)    # (m, Nens)

    #             # ----------------------------------------------------------------------
    #             # Compute Dprime and HAprime (final outputs)
    #             # ----------------------------------------------------------------------
    #             Dprime  = (d_global - Hbar_global)          # (m,1)
    #             HAprime = Eta                                                # (m,Nens)

    #             # Optional: diagnostics
    #             if rank == 0:
    #                 print(f"[rank {rank}] Shapes -> Dprime: {Dprime.shape}, HAprime: {HAprime.shape}")

    #             return Dprime, Eta, HAprime, kwargs

    #     except Exception as e:
    #         print(f"[rank {rank}] Error in compute_X5_utils: {e}")
    #         tb_str = "".join(traceback.format_exception(*sys.exc_info()))
    #         print(f"Traceback details:\n{tb_str}")
    #         comm.Abort(1)


    