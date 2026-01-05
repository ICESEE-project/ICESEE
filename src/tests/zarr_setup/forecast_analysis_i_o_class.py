import h5py
import numpy as np
from mpi4py import MPI
import os
import glob
import gc
import zarr
import traceback
import sys
import shutil
from numcodecs import blosc
import time
import functools

blosc.use_threads = False

from typing import Callable, Optional, TypeVar, Any

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

class EnKFIO:
    def __init__(self, file_prefix, nd, nens, nt, subcomm, mpi_comm, params, serial_file_creation=True, base_path="enkf_data", batch_size=50):
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

            # # Divide nd among ranks
            # nd_local_base = nd // self.size
            # remainder = nd % self.size
            # if self.rank < remainder:
            #     self.nd_local = nd_local_base + 1
            #     self.nd_start = self.rank * (nd_local_base + 1)
            # else:
            #     self.nd_local = nd_local_base
            #     self.nd_start = remainder * (nd_local_base + 1) + (self.rank - remainder) * nd_local_base
            # self.nd_end = self.nd_start + self.nd_local

            # # foranalysis lets use mpi_com
            # size_world = mpi_comm.Get_size()
            # rank_world = mpi_comm.Get_rank()
            # nd_local_base_world = nd // size_world
            # remainder_world = nd % size_world
            # if rank_world < remainder_world:
            #     self.nd_local_world = nd_local_base_world + 1
            #     self.nd_start_world = rank_world * (nd_local_base_world + 1)
            # else:
            #     self.nd_local_world = nd_local_base_world
            #     self.nd_start_world = remainder_world * (nd_local_base_world + 1) + (rank_world - remainder_world) * nd_local_base_world
            # self.nd_end_world = self.nd_start_world + self.nd_local_world

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
            self.comm.Barrier()

            # Initialize file and dataset lists
            self.files = []
            self.datasets = []
            self.current_batch_start = -1

            # Create initial batch
            if self.serial_file_creation:
                self._create_batch_serial(0)
            else:
                self._create_batch_parallel(0)
        except Exception as e:
            print(f"Error occurred in __init__: {e}")
            tb_str = "".join(traceback.format_exception(*sys.exc_info()))
            print(f"Traceback details:\n{tb_str}")
            self.mpi_comm.Abort(1)

    def _create_batch_serial(self, t_start):
        try:
            start = MPI.Wtime()
            self._close_batch()
            self.files = []
            self.datasets = []
            self.current_batch_start = t_start
            nfiles = min(self.batch_size, self.nt - t_start)

            if MPI.COMM_WORLD.Get_rank() == 0:
                for t in range(t_start, t_start + nfiles):
                    fname = f"{self.base_path}/{self.file_prefix}_{t:04d}.h5"
                    with h5py.File(fname, 'w') as f:
                        row_chunk = min(1024, self.nd)
                        col_chunk = 1
                        f.create_dataset(
                            'states', (self.nd, self.nens),
                            chunks=(row_chunk, col_chunk),
                            # compression="gzip", compression_opts=4,
                            compression="lzf",
                            dtype='f8'
                        )
                        # f.create_dataset(
                        #     'states', (self.nd, self.nens),
                        #     chunks=(row_chunk, col_chunk),
                        #     dtype='f8'
                        # )
                        
            self.mpi_comm.Barrier()

            for t in range(t_start, t_start + nfiles):
                fname = f"{self.base_path}/{self.file_prefix}_{t:04d}.h5"
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
                row_chunk = min(1024, self.nd)
                col_chunk = 1
                dset = f.create_dataset(
                    'states', (self.nd, self.nens),
                    chunks=(row_chunk, col_chunk),
                    compression="gzip", compression_opts=4,
                    dtype='f8'
                )
                self.files.append(f)
                self.datasets.append(dset)
        except Exception as e:
            print(f"Error occurred in _create_batch_parallel: {e}")
            tb_str = "".join(traceback.format_exception(*sys.exc_info()))
            print(f"Traceback details:\n{tb_str}")
            self.mpi_comm.Abort(1)

    def _close_batch(self):
        try:
            for f in self.files:
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
                if self.serial_file_creation:
                    self._create_batch_serial(batch_start)
                else:
                    self._create_batch_parallel(batch_start)
        except Exception as e:
            print(f"Error occurred in _ensure_batch: {e}")
            tb_str = "".join(traceback.format_exception(*sys.exc_info()))
            print(f"Traceback details:\n{tb_str}")
            self.mpi_comm.Abort(1)

    def read_forecast(self, t, ens_idx):
        try:
            self._ensure_batch(t)
            batch_idx = t - self.current_batch_start
            start = MPI.Wtime()
            data = self.datasets[batch_idx][self.nd_start:self.nd_end, ens_idx]
            read_time = MPI.Wtime() - start
            return data
        except Exception as e:
            print(f"Error occurred in read_forecast: {e}")
            tb_str = "".join(traceback.format_exception(*sys.exc_info()))
            print(f"Traceback details:\n{tb_str}")
            self.mpi_comm.Abort(1)

    def write_forecast(self, t, data, ens_idx):
        try:
            self._ensure_batch(t)
            batch_idx = t - self.current_batch_start
            start = MPI.Wtime()
            self.datasets[batch_idx][self.nd_start:self.nd_end, ens_idx] = data
            write_time = MPI.Wtime() - start
        except Exception as e:
            print(f"Error occurred in write_forecast: {e}")
            tb_str = "".join(traceback.format_exception(*sys.exc_info()))
            print(f"Traceback details:\n{tb_str}")
            self.mpi_comm.Abort(1)

    def read_analysis(self, t, ens_idx):
        try:
            self._ensure_batch(t)
            batch_idx = t - self.current_batch_start
            start = MPI.Wtime()
            data = self.datasets[batch_idx][self.nd_start_world:self.nd_end_world, ens_idx]
            read_time = MPI.Wtime() - start
            return data
        except Exception as e:
            print(f"Error occurred in read_analysis: {e}")
            tb_str = "".join(traceback.format_exception(*sys.exc_info()))
            print(f"Traceback details:\n{tb_str}")
            self.mpi_comm.Abort(1)

    def write_analysis(self, t, data, ens_idx):
        try:
            self._ensure_batch(t)
            batch_idx = t - self.current_batch_start
            start = MPI.Wtime()
            # self.datasets[batch_idx][self.nd_start:self.nd_end, ens_idx] = data
            self.datasets[batch_idx][self.nd_start_world:self.nd_end_world, ens_idx] = data
            write_time = MPI.Wtime() - start
        except Exception as e:
            print(f"Error occurred in write_analysis: {e}")
            tb_str = "".join(traceback.format_exception(*sys.exc_info()))
            print(f"Traceback details:\n{tb_str}")
            self.mpi_comm.Abort(1)

    # def read_analysis(self, t, ens_idx, start=None, end=None):
    #     try:
    #         self._ensure_batch(t)
    #         batch_idx = t - self.current_batch_start
    #         start = MPI.Wtime()
    #         if start is None or end is None:
    #             data = self.datasets[batch_idx][self.nd_start_world:self.nd_end_world, ens_idx]
    #         else:
    #             data = self.datasets[batch_idx][self.nd_start_world + start:self.nd_start_world + end, ens_idx]
    #         read_time = MPI.Wtime() - start
    #         return data
    #     except Exception as e:
    #         print(f"Error occurred in read_analysis: {e}")
    #         tb_str = "".join(traceback.format_exception(*sys.exc_info()))
    #         print(f"Traceback details:\n{tb_str}")
    #         self.mpi_comm.Abort(1)

    # def write_analysis(self, t, data, ens_idx, offset=None):
    #     try:
    #         self._ensure_batch(t)
    #         batch_idx = t - self.current_batch_start
    #         start = MPI.Wtime()
    #         if offset is None:
    #             self.datasets[batch_idx][self.nd_start:self.nd_end, ens_idx] = data
    #         else:
    #             self.datasets[batch_idx][self.nd_start + offset:self.nd_start + offset + len(data), ens_idx] = data
    #         write_time = MPI.Wtime() - start
    #     except Exception as e:
    #         print(f"Error occurred in write_analysis: {e}")
    #         tb_str = "".join(traceback.format_exception(*sys.exc_info()))
    #         print(f"Traceback details:\n{tb_str}")
    #         self.mpi_comm.Abort(1)

    def write_matrix(self, t, dataset_name, data, ens_idx):
        try:
            self._ensure_batch(t)
            batch_idx = t - self.current_batch_start
            start = MPI.Wtime()
            self.files[batch_idx][dataset_name][self.nd_start:self.nd_end, ens_idx] = data
            write_time = MPI.Wtime() - start
        except Exception as e:
            print(f"Error occurred in write_matrix: {e}")
            tb_str = "".join(traceback.format_exception(*sys.exc_info()))
            print(f"Traceback details:\n{tb_str}")
            self.mpi_comm.Abort(1)

    def read_matrix(self, t, dataset_name, ens_idx):
        try:
            self._ensure_batch(t)
            batch_idx = t - self.current_batch_start
            start = MPI.Wtime()
            data = self.files[batch_idx][dataset_name][self.nd_start:self.nd_end, ens_idx]
            read_time = MPI.Wtime() - start
            return data
        except Exception as e:
            print(f"Error occurred in read_matrix: {e}")
            tb_str = "".join(traceback.format_exception(*sys.exc_info()))
            print(f"Traceback details:\n{tb_str}")
            self.mpi_comm.Abort(1)

    def gather_matrix(self, t, dataset_name):
        try:
            self._ensure_batch(t)
            batch_idx = t - self.current_batch_start
            start = MPI.Wtime()
            local_data = self.files[batch_idx][dataset_name][self.nd_start:self.nd_end, :]
            counts = self.mpi_comm.allgather(self.nd_local)
            displacements = self.mpi_comm.allgather(self.nd_start)
            global_data = np.zeros((self.nd, self.nens), dtype='f8')
            self.mpi_comm.Allgatherv(local_data, [global_data, counts, displacements, MPI.DOUBLE])
            gather_time = MPI.Wtime() - start
            return global_data
        except Exception as e:
            print(f"Error occurred in gather_matrix: {e}")
            tb_str = "".join(traceback.format_exception(*sys.exc_info()))
            print(f"Traceback details:\n{tb_str}")
            self.mpi_comm.Abort(1)

    def compute_forecast_mean(self, t):
        try:
            self._ensure_batch(t)
            comm = self.mpi_comm
            rank = comm.Get_rank()
            size = comm.Get_size()

            batch_idx = t - self.current_batch_start
            start = MPI.Wtime()

            self.mpi_comm.barrier()
            local_data = self.datasets[batch_idx][self.nd_start_world:self.nd_end_world, :]
            local_mean = np.mean(local_data, axis=1).astype('f8', copy=False)

            local_row_count = np.array([local_mean.shape[0]], dtype='i8')
            global_row_counts = np.zeros(size, dtype='i8') if rank == 0 else None
            comm.Gather(local_row_count, global_row_counts, root=0)

            if rank == 0:
                displacements = np.zeros(size, dtype='i8')
                displacements[1:] = np.cumsum(global_row_counts[:-1])
                total_rows = np.sum(global_row_counts)
                global_mean = np.zeros(total_rows, dtype='f8')
            else:
                displacements = None
                global_mean = None

            comm.Gatherv(local_mean, [global_mean, global_row_counts, displacements, MPI.DOUBLE], root=0)

            if rank == 0:
                result = global_mean
                file_path = f"{self.base_path}/{self.file_prefix}_mean.h5"
                with h5py.File(file_path, 'a') as f:
                    if 'mean' not in f:
                        f.create_dataset(
                            'mean', (self.nd, self.nt),
                            chunks=(min(self.nd,1000),1),
                            dtype='f8'
                        )
                    f['mean'][:, t] = result
            self.mpi_comm.barrier()
        except Exception as e:
            print(f"Error occurred in compute_forecast_mean: {e}")
            tb_str = "".join(traceback.format_exception(*sys.exc_info()))
            print(f"Traceback details:\n{tb_str}")
            self.mpi_comm.Abort(1)

    @retry_on_failure(max_attempts=5, delay=1.0, mpi_comm=MPI.COMM_WORLD)
    def compute_forecast_mean_chunked_gather(self, t, ens_chunk_size=1):
        try:
            self._ensure_batch(t)
            comm = self.mpi_comm
            rank = comm.Get_rank()
            size = comm.Get_size()

            batch_idx = t - self.current_batch_start
            start = MPI.Wtime()

            local_rows = self.nd_end_world - self.nd_start_world
            local_sum = np.zeros(local_rows, dtype='f8')

            for start_ens in range(0, self.nens, ens_chunk_size):
                end_ens = min(start_ens + ens_chunk_size, self.nens)
                local_data = self.datasets[batch_idx][self.nd_start_world:self.nd_end_world, start_ens:end_ens]
                partial_sum = np.sum(local_data, axis=1).astype('f8', copy=False)
                local_sum += partial_sum

            local_mean = local_sum / self.nens

            local_row_count = np.array([local_mean.shape[0]], dtype='i8')
            global_row_counts = np.zeros(size, dtype='i8') if rank == 0 else None
            comm.Gather(local_row_count, global_row_counts, root=0)

            if rank == 0:
                displacements = np.zeros(size, dtype='i8')
                displacements[1:] = np.cumsum(global_row_counts[:-1])
                total_rows = np.sum(global_row_counts)
                global_mean = np.zeros(total_rows, dtype='f8')
            else:
                displacements = None
                global_mean = None

            comm.Gatherv(local_mean, [global_mean, global_row_counts, displacements, MPI.DOUBLE], root=0)

            if rank == 0:
                result = global_mean
                file_path = f"{self.base_path}/{self.file_prefix}_mean.h5"
                with h5py.File(file_path, 'a') as f:
                    if 'mean' not in f:
                        f.create_dataset(
                            'mean', (self.nd, self.nt),
                            chunks=(min(self.nd,1000),1),
                            dtype='f8'
                        )
                    f['mean'][:, t] = result
        except Exception as e:
            print(f"Error occurred in compute_forecast_mean_chunked_gather: {e}")
            tb_str = "".join(traceback.format_exception(*sys.exc_info()))
            print(f"Traceback details:\n{tb_str}")
            self.mpi_comm.Abort(1)
    
    def compute_forecast_mean_chunked_v2(self, k, **kwargs):
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
            for ens_idx in range(self.nens):
                if local_rows > 0:
                    v = self.read_analysis(k, ens_idx)
                    v = np.asarray(v, dtype='f8')
                    if v.ndim != 1 or v.size != local_rows:
                        v = v.reshape(-1)
                        assert v.size == local_rows, "read_analysis must return (local_rows,)"
                    local_sum += v

            local_mean = (local_sum / float(self.nens)) if local_rows > 0 else np.empty((0,), dtype='f8')

            # ---- Parallel HDF5: collective create + collective write -------------
            file_path = f"{self.base_path}/{self.file_prefix}_mean.h5"

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

    # @retry_on_failure(max_attempts=5, delay=1.0, mpi_comm=MPI.COMM_WORLD)
    def compute_forecast_mean_chunked(self, t, ens_chunk_size=None, use_collective_io=False, max_ranks=4):
        self._ensure_batch(t)
        comm = self.mpi_comm
        rank = comm.Get_rank()
        size = comm.Get_size()

        import numpy as np
        import h5py.h5p as h5p
        import h5py.h5s as h5s
        import h5py.h5fd as h5fd
        try:
            if size == 2:
                max_ranks = 1
                active_ranks = 1
            else:
                max_ranks = min(8, (size // 2) + 1)
                active_ranks = min(max_ranks, size)

            # Exclude ranks with no work
            color = 0 if (rank < active_ranks) else MPI.UNDEFINED
            sub_comm = comm.Split(color, rank)
            active = sub_comm != MPI.COMM_NULL

            if active:
                sub_rank = sub_comm.Get_rank()
                sub_size = sub_comm.Get_size()

                local_rows = (self.nd_end_world - self.nd_start_world) if rank == 0 else 0
                local_rows = sub_comm.bcast(local_rows, root=0)
                rows_per_rank = local_rows // sub_size
                extra_rows = local_rows % sub_size
                nd_start = sub_rank * rows_per_rank + min(sub_rank, extra_rows)
                nd_end = nd_start + rows_per_rank + (1 if sub_rank < extra_rows else 0)

                batch_idx = t - self.current_batch_start
                ds = self.datasets[batch_idx]

                if ens_chunk_size is None:
                    bytes_per_element = 8
                    target_memory = 1e9
                    ens_chunk_size = max(1, int(target_memory / (nd_end - nd_start) / bytes_per_element))
                    ens_chunk_size = min(ens_chunk_size, self.nens)
                    if sub_rank == 0:
                        print(f"Dynamic ens_chunk_size: {ens_chunk_size}")

                local_sum = np.zeros(nd_end - nd_start, dtype='f8')

                dxpl = h5p.create(h5p.DATASET_XFER)
                if use_collective_io:
                    dxpl.set_dxpl_mpio(h5fd.MPIO_COLLECTIVE)
                else:
                    dxpl.set_dxpl_mpio(h5fd.MPIO_INDEPENDENT)

                t_io = 0.0
                t_comp = 0.0
                start_total = MPI.Wtime()

                buffers = [np.empty((nd_end - nd_start, ens_chunk_size), dtype='f8') for _ in range(2)]
                current_buffer = 0

                for start_ens in range(0, self.nens, ens_chunk_size):
                    end_ens = min(start_ens + ens_chunk_size, self.nens)
                    chunk_cols = end_ens - start_ens

                    if chunk_cols < ens_chunk_size:
                        buffers[current_buffer] = np.empty((nd_end - nd_start, chunk_cols), dtype='f8')

                    t_start_io = MPI.Wtime()
                    file_space = ds.id.get_space()
                    file_space.select_hyperslab((self.nd_start_world + nd_start, start_ens), (nd_end - nd_start, chunk_cols))
                    mem_space = h5s.create_simple((nd_end - nd_start, chunk_cols))
                    ds.id.read(mem_space, file_space, buffers[current_buffer], dxpl=dxpl)
                    t_io += MPI.Wtime() - t_start_io

                    if start_ens > 0:
                        t_start_comp = MPI.Wtime()
                        local_sum += np.sum(buffers[1 - current_buffer], axis=1)
                        t_comp += MPI.Wtime() - t_start_comp

                    current_buffer = 1 - current_buffer

                t_start_comp = MPI.Wtime()
                local_sum += np.sum(buffers[current_buffer], axis=1)
                t_comp += MPI.Wtime() - t_start_comp

                local_mean = local_sum / self.nens

                file_path = f"{self.base_path}/{self.file_prefix}_mean.h5"
                sub_comm.Barrier()
                t_start_io = MPI.Wtime()
                f = h5py.File(file_path, 'a', driver='mpio', comm=sub_comm)

                if 'mean' not in f:
                    chunk_rows = min(self.nd, 1000)
                    f.create_dataset(
                        'mean', (self.nd, self.nt),
                        chunks=(chunk_rows, 1),
                        dtype='f8'
                    )

                out_ds = f['mean']
                file_space = out_ds.id.get_space()
                file_space.select_hyperslab((self.nd_start_world + nd_start, t), (nd_end - nd_start, 1))
                mem_space = h5s.create_simple((nd_end - nd_start,))
                out_ds.id.write(mem_space, file_space, local_mean, dxpl=dxpl)

                f.close()
                t_io += MPI.Wtime() - t_start_io
                sub_comm.Barrier()

                if sub_rank == 0:
                    print(f"Total time: {MPI.Wtime() - start_total:.2f}s, I/O: {t_io:.2f}s, Compute: {t_comp:.2f}s")

            comm.Barrier()
        except Exception as e:
            print(f"Error occurred in compute_forecast_mean_chunked: {e}")
            tb_str = "".join(traceback.format_exception(*sys.exc_info()))
            print(f"Traceback details:\n{tb_str}")
            self.mpi_comm.Abort(1)

    def _compute_forecast_mean_chunked(self, t, ens_chunk_size=1):
        try:
            self._ensure_batch(t)
            comm = self.mpi_comm
            rank = comm.Get_rank()
            size = comm.Get_size()

            batch_idx = t - self.current_batch_start
            start = MPI.Wtime()

            local_rows = self.nd_end_world - self.nd_start_world
            local_sum = np.zeros(local_rows, dtype='f8')

            for start_ens in range(0, self.nens, ens_chunk_size):
                end_ens = min(start_ens + ens_chunk_size, self.nens)
                local_data = self.datasets[batch_idx][self.nd_start_world:self.nd_end_world, start_ens:end_ens]
                partial_sum = np.sum(local_data, axis=1).astype('f8', copy=False)
                local_sum += partial_sum

            local_mean = local_sum / self.nens

            file_path = f"{self.base_path}/{self.file_prefix}_mean.h5"
            with h5py.File(file_path, 'a', driver='mpio', comm=comm) as f:
                if 'mean' not in f:
                    if rank == 0:
                        f.create_dataset(
                            'mean', (self.nd, self.nt),
                            chunks=(min(self.nd, 1000), 1),
                            dtype='f8'
                        )
                comm.Barrier()
        except Exception as e:
            print(f"Error occurred in _compute_forecast_mean_chunked: {e}")
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

    def icesee_get_index(self, **kwargs):
        try:
            vec_inputs = kwargs.get("vec_inputs", None)
            params = kwargs.get("params", None)
            nd = kwargs.get("nd", params.get("nd", None))

            if params["default_run"]:
                comm = kwargs.get("subcomm", None)
            else:
                comm = kwargs.get("comm_world", None)
            
            len_vec = params["total_state_param_vars"]
            dim_list_param = np.array(kwargs.get('dim_list', None)) // len_vec
            hdim = nd // len_vec

            if comm is None:
                rank = 0
                dim = dim_list_param[rank]
                offsets = [0]
            else:
                if params["even_distribution"]:
                    rank = 0
                    dim = dim_list_param[rank]
                    offsets = [0]
                else:
                    rank = comm.Get_rank()
                    dim = dim_list_param[rank]
                    offsets = np.cumsum(np.insert(dim_list_param, 0, 0))

            start_idx = offsets[rank]
            index_map = {}
            var_start = 0

            for var in vec_inputs:
                start = var_start + start_idx
                end = start + dim
                index_map[var] = np.arange(start, end)
                var_start += hdim

            local_size_per_rank = kwargs.get('dim_list', None)
            return index_map, local_size_per_rank[rank]
        except Exception as e:
            print(f"Error occurred in icesee_get_index: {e}")
            tb_str = "".join(traceback.format_exception(*sys.exc_info()))
            print(f"Traceback details:\n{tb_str}")
            self.mpi_comm.Abort(1)

    @retry_on_failure(max_attempts=5, delay=1.0, mpi_comm=MPI.COMM_WORLD)
    def generate_observation_schedule(self, **kwargs):
        try:
            t = np.array(kwargs["t"])
            freq_obs = self.params["freq_obs"]
            obs_start_time = self.params["obs_start_time"]
            obs_max_time = self.params["obs_max_time"]

            max_t = np.max(t)
            obs_max_time = min(obs_max_time, max_t)

            obs_t = np.arange(obs_start_time, obs_max_time + freq_obs, freq_obs)
            obs_t = obs_t[obs_t <= obs_max_time]

            obs_idx = []
            for time in obs_t:
                idx = np.argmin(np.abs(t - time))
                obs_idx.append(idx)
            obs_idx = np.array(obs_idx, dtype=int)

            num_observations = len(obs_idx)
            return obs_t, obs_idx, num_observations
        except Exception as e:
            print(f"Error occurred in generate_observation_schedule: {e}")
            tb_str = "".join(traceback.format_exception(*sys.exc_info()))
            print(f"Traceback details:\n{tb_str}")
            self.mpi_comm.Abort(1)

    @retry_on_failure(max_attempts=5, delay=1.0, mpi_comm=MPI.COMM_WORLD)
    def _create_synthetic_observations(self, **kwargs):
        try:
            synthetic_obs_zarr_path = kwargs.get('synthetic_obs_zarr_path')
            error_R_zarr_path = kwargs.get('error_R_zarr_path')
            nd = kwargs.get('nd')
            nt = kwargs.get('nt')

            obs_t, ind_m, m_obs = self.generate_observation_schedule(**kwargs)
            m = m_obs
            m_R = m_obs*2 +1

            # print(f"\n m={m}, m_R={m_R} \n")

            rank = self.mpi_comm.Get_rank()
            size = self.mpi_comm.Get_size()

            if rank == 0:
                if os.path.exists(synthetic_obs_zarr_path):
                    shutil.rmtree(synthetic_obs_zarr_path)
                if os.path.exists(error_R_zarr_path):
                    shutil.rmtree(error_R_zarr_path)
            self.mpi_comm.Barrier()

            if rank == 0:
                hu_obs = zarr.create_array(store=synthetic_obs_zarr_path, shape=(nd, m), chunks=(min(1000, nd), min(50, m)), dtype='f8', overwrite=True)
                error_R = zarr.create_array(store=error_R_zarr_path, shape=(nd, m_R), chunks=(min(1000, nd), min(50, m_R)), dtype='f8', overwrite=True)
            
            self.mpi_comm.Barrier()
            hu_obs = zarr.open_array(store=synthetic_obs_zarr_path, mode='r+')
            error_R = zarr.open_array(store=error_R_zarr_path, mode='r+')
            self.mpi_comm.Barrier()

            if kwargs.get('joint_estimation', False) or self.params.get('localization_flag', False):
                hdim = nd // self.params["total_state_param_vars"]
            else:
                hdim = nd // self.params["total_state_param_vars"]

            if rank == 0:
                for i, sig in enumerate(self.params["sig_obs"]):
                    start_idx = i*hdim
                    end_idx = start_idx + hdim
                    error_R[start_idx:end_idx,:] = np.ones((hdim,1)) * sig
            self.mpi_comm.Barrier()

            statevec_true = zarr.open_array(store="output/statevec_true.zarr", mode='r+')
            indx_map, _ = self.icesee_get_index(**kwargs)
            if self.nd < 10000:
                if rank==0:
                    km = 0
                    for step in range(nt):
                        if (km<m_obs) and (step+1 == ind_m[km]):
                            for key in kwargs['vec_inputs']:
                                # hu_obs[indx_map[key],km] = statevec_true[indx_map[key],step+1]
                                hu_obs[indx_map[key],km] = statevec_true[indx_map[key],step+1] + np.random.normal(0,error_R[indx_map[key],km],len(indx_map[key]))

                            km += 1
                self.mpi_comm.Barrier()
            else:
                if size >= m_obs:
                    obs_per_process = m_obs // size
                    remainder = m_obs % size
                    start_obs = rank * obs_per_process + min(rank, remainder)
                    num_obs = obs_per_process + 1 if rank < remainder else obs_per_process

                    rows_per_process = hdim // size
                    row_remainder = hdim % size
                    row_start = rank * rows_per_process + min(rank, row_remainder)
                    row_end = row_start + (rows_per_process + 1 if rank < row_remainder else rows_per_process)

                    for km in range(start_obs, start_obs + num_obs):
                        if km < m_obs:
                            step = ind_m[km] - 1
                            if 0 <= step < nt:
                                for key in kwargs['vec_inputs']:
                                    indices = indx_map[key]
                                    local_indices = indices[(indices >= row_start) & (indices < row_end)]
                                    if len(local_indices) > 0:
                                        state_data = statevec_true[local_indices, step]
                                        error_data = error_R[local_indices, km]
                                        result = state_data + np.random.normal(0, error_data, len(local_indices))
                                        if result.shape != (len(local_indices),):
                                            raise ValueError(f"Rank {rank}: Shape mismatch at km={km}: expected {len(local_indices)}, got {result.shape}")
                                        hu_obs[local_indices, km] = result
                    self.mpi_comm.Barrier()
                else:
                    obs_per_process = m_obs // size
                    remainder = m_obs % size
                    start_obs = rank * obs_per_process + min(rank, remainder)
                    num_obs = obs_per_process + 1 if rank < remainder else obs_per_process

                    rows_per_process = nd // size
                    row_remainder = nd % size
                    row_start = rank * rows_per_process + min(rank, row_remainder)
                    row_end = min(row_start + (rows_per_process + 1 if rank < row_remainder else rows_per_process), nd)

                    for km in range(start_obs, start_obs + num_obs):
                        if km < m_obs:
                            step = ind_m[km] - 1
                            if 0 <= step < nt:
                                for key in kwargs['vec_inputs']:
                                    indices = indx_map[key]
                                    local_indices = indices[(indices >= row_start) & (indices < row_end)]
                                    if len(local_indices) > 0:
                                        state_data = statevec_true[local_indices, step]
                                        error_data = error_R[local_indices, km]
                                        result = state_data + np.random.normal(0, error_data, len(local_indices))
                                        if result.shape != (len(local_indices),):
                                            raise ValueError(f"Rank {rank}: Shape mismatch at km={km}: expected {len(local_indices)}, got {result.shape}")
                                        hu_obs[local_indices, km] = result
                    self.mpi_comm.Barrier()

            return obs_t, m_obs
        except Exception as e:
            print(f"Error in _create_synthetic_observations: {e}")
            tb_str = "".join(traceback.format_exception(*sys.exc_info()))
            print(f"Traceback details:\n{tb_str}")
            self.mpi_comm.Abort(1)

    def H_matrix(self, **kwargs):
        try:
            zarr_path = kwargs.get('H_matrix_zarr_path')
            nd = kwargs.get('nd')
            m_obs = kwargs.get('m_obs')
            m = m_obs * 2 + 1
            di = int((nd - 2) / (2 * m_obs))

            H_matrix_file = zarr.create_array(store=zarr_path, shape=(m, nd), chunks=(min(50, m), min(1000, nd)), dtype='f8', overwrite=True)
            for i in range(1, m_obs + 1):
                H_matrix_file[i - 1, i * di - 1] = 1
                H_matrix_file[m_obs + i - 1, int((nd - 2) / 2) + i * di - 1] = 1

            H_matrix_file[m_obs * 2, nd - 2] = 1

            if self.params.get('joint_estimation', False):
                ndim = nd // self.params["total_state_param_vars"]
                state_variables_size = ndim * self.params["num_state_vars"]
                H_matrix_file[:, state_variables_size:] = 0
        except Exception as e:
            print(f"Error in H_matrix: {e}")
            tb_str = "".join(traceback.format_exception(*sys.exc_info()))
            print(f"Traceback details:\n{tb_str}")
            self.mpi_comm.Abort(1)

    def Eta_matrix(self, k, ens_idx, **kwargs):
        try:
            H_matrix_zarr_path = kwargs.get('H_matrix_zarr_path', "output/H_matrix.zarr")
            # Eta_matrix_zarr_path = kwargs.get('Eta_matrix_zarr_path', "output/Eta_matrix.zarr")

            H_matrix_file = zarr.open_array(H_matrix_zarr_path, mode='r')
            H_local = H_matrix_file[:, self.nd_start_world:self.nd_end_world]

            mean_file_path = f"{self.base_path}/{self.file_prefix}_mean.h5"
            with h5py.File(mean_file_path, 'r', driver='mpio', comm=self.mpi_comm) as f:
                ens_mean = f['mean'][self.nd_start_world:self.nd_end_world, k]
                
            state = self.read_analysis(k, ens_idx)
            ens_pertubations = state - ens_mean

            Eta_local = np.dot(H_local, ens_pertubations)

            # compute Eta_global from all ranks 
            Eta_global = np.empty_like(Eta_local)
            self.mpi_comm.Allreduce(Eta_local, Eta_global, op=MPI.SUM)

            return Eta_global
        except Exception as e:
            print(f"Error in Eta_matrix: {e}")
            tb_str = "".join(traceback.format_exception(*sys.exc_info()))
            print(f"Traceback details:\n{tb_str}")
            self.mpi_comm.Abort(1)
    
    def Eta_matrix_root(self, k, ens_idx, **kwargs):
        try:
            H_matrix_zarr_path = kwargs.get('H_matrix_zarr_path', "output/H_matrix.zarr")
            # Eta_matrix_zarr_path = kwargs.get('Eta_matrix_zarr_path', "output/Eta_matrix.zarr")

            H_matrix_file = zarr.open_array(H_matrix_zarr_path, mode='r')
            H_local = H_matrix_file[:, self.nd_start_world:self.nd_end_world]

            mean_file_path = f"{self.base_path}/{self.file_prefix}_mean.h5"
            with h5py.File(mean_file_path, 'r', driver='mpio', comm=self.mpi_comm) as f:
                ens_mean = f['mean'][self.nd_start_world:self.nd_end_world, k]
                
            state = self.read_analysis(k, ens_idx)
            ens_pertubations = state - ens_mean

            Eta_local = np.dot(H_local, ens_pertubations)

            # compute Eta_global from all ranks 
            Eta_global = np.empty_like(Eta_local)
            self.mpi_comm.Allreduce(Eta_local, Eta_global, op=MPI.SUM)

            return Eta_global
        except Exception as e:
            print(f"Error in Eta_matrix: {e}")
            tb_str = "".join(traceback.format_exception(*sys.exc_info()))
            print(f"Traceback details:\n{tb_str}")
            self.mpi_comm.Abort(1)

    def compute_X5_root(self, km, **kwargs):
        try:
            H_matrix_zarr_path = kwargs.get('H_matrix_zarr_path', "output/H_matrix.zarr")
            synthetic_obs_zarr_path = kwargs.get('synthetic_obs_zarr_path', "output/synthetic_obs.zarr")
            # d_matrix_zarr_path = kwargs.get('d_matrix_zarr_path', "output/d_matrix.zarr")

            H_matrix = zarr.open_array(H_matrix_zarr_path, mode='r')
            H_local = H_matrix[:, self.nd_start_world:self.nd_end_world]

            #  get the synthetic observations
            synthetic_obs = zarr.open_array(synthetic_obs_zarr_path, mode='r')
            synthetic_obs_local = synthetic_obs[self.nd_start_world:self.nd_end_world, km]
            # print(f"\n synthetic_obs shape: {synthetic_obs.shape} synthetic_obs_local shape: {synthetic_obs_local.shape}\n")

            rank = self.mpi_comm.Get_rank()

            #  compute the d matrix
            # d_local = zarr.create_array(store=f"output/d_matrix_{rank}.zarr", shape=(H_matrix.shape[0],), chunks=(min(50, H_matrix.shape[0]),), dtype='f8', overwrite=True)
            d_local = np.dot(H_local, synthetic_obs_local)

            # compute global d from all ranks
            d_global = np.empty_like(d_local)
            # d_global = za
            self.mpi_comm.Allreduce(d_local, d_global, op=MPI.SUM)

            # read all ensemble members locally into a zarr array
            state_local_matrix = zarr.create_array(store=f"output/state_local_matrix_{rank}.zarr", shape=(self.nd_end_world - self.nd_start_world, self.nens), chunks=(min(1000, self.nd_end_world - self.nd_start_world), min(10, self.nens)), dtype='f8', overwrite=True)
            for ens_idx in range(self.nens):
                state_local = self.read_analysis(km, ens_idx)
                state_local_matrix[:, ens_idx] = state_local

            # -- compute innovations D - HA for all ensemble members
            HA_local = np.dot(H_local, state_local_matrix)  # (m, nens)
            HA = np.empty_like(HA_local)
            self.mpi_comm.Allreduce(HA_local, HA, op=MPI.SUM)

            # compute Eta and D for all ensemble members
            # get the ensemble mean
            mean_file_path = f"{self.base_path}/{self.file_prefix}_mean.h5"
            with h5py.File(mean_file_path, 'r', driver='mpio', comm=self.mpi_comm) as f:
                ens_mean = f['mean'][self.nd_start_world:self.nd_end_world, km]
            ens_pertubations = state_local_matrix - ens_mean[:, np.newaxis]  # (nd_local, nens)
            Eta_local = np.dot(H_local, ens_pertubations)  # (m, nens)
            Eta = np.empty_like(Eta_local)
            self.mpi_comm.Allreduce(Eta_local, Eta, op=MPI.SUM)

            D_global = d_global + Eta  # (m, nens)
            Dprime = D_global - HA  # (m, nens)

            # -- create a zarr file for each processor to write to  
            m = kwargs.get('m_obs')*2+1
            Nens = self.nens
        
            X5 = np.empty((Nens, Nens))
            if self.mpi_comm.Get_rank() == 0:
            # if False:
                # compute the HAbar
                # HAbar = np.mean(HA, axis=1)
                # HAprime_local = HA - HAbar[:, np.newaxis]
                one_N = np.ones((Nens,Nens))/Nens
                HAprime= HA@(np.eye(Nens) - one_N) # mxNens

                # compute HAprime + Eta
                HAprime_Eta = HAprime + Eta
                print(f"[Rank {self.mpi_comm.Get_rank()}] HAprime_Eta_local norm: {np.linalg.norm(HAprime_Eta)}, shape: {HAprime_Eta.shape}")
                # print(f"\n [Rank {self.mpi_comm.Get_rank()}] Dprime_local shape: {Dprime_local.shape} HAprime_local shape: {HAprime_local.shape} HAprime_Eta_local shape: {HAprime_Eta_local.shape}\n ")

                # compute SVD of HAprime_Eta
                U, sig, Vt = np.linalg.svd(HAprime_Eta, full_matrices=False)

                # get the min (m Nens)
                nrmin = min(Nens, m)
                
                # convert S to eigenvalues
                sig = sig**2

                # compute the number of significant eigenvalues
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

                #  compute X3 = U*X2 # m_obs x Nens
                X3 = np.dot(U, X2)

                # print(f"[ICESEE] Rank: {rank_world} X3 shape: {X3.shape}")
                # compute X4 = (HAprime.T)*X3 # Nens x Nens
                X4 = np.dot(HAprime.T, X3)
                del X2, X3, U, HAprime, HA, Eta, Dprime
                gc.collect()

                # compute X5 = X4 + I
                X5 = X4 + np.eye(Nens)
                # sum of each column of X5 should be 1
                if np.sum(X5, axis=0).all() != 1.0:
                    print(f"[ICESEE] Sum of each X5 column is not 1.0: {np.sum(X5, axis=0)}")
                # print(f"[ICESEE] Rank: {comm_world.Get_rank()} X5 sum: {np.sum(X5, axis=0)}")
                del X4; gc.collect()
            
            # Broadcast X5 to all ranks
            self.mpi_comm.Bcast(X5, root=0)

            return X5
        
        except Exception as e:
            print(f"Error in compute_X5_analysis_mean: {e}")
            tb_str = "".join(traceback.format_exception(*sys.exc_info()))
            print(f"Traceback details:\n{tb_str}")
            self.mpi_comm.Abort(1)

    @retry_on_failure(max_attempts=5, delay=1.0, mpi_comm=MPI.COMM_WORLD)
    def compute_X5_(self, km, **kwargs):
        # innovation (Dprime) = D - HA 
        try:
            H_matrix_zarr_path = kwargs.get('H_matrix_zarr_path', "output/H_matrix.zarr")
            synthetic_obs_zarr_path = kwargs.get('synthetic_obs_zarr_path', "output/synthetic_obs.zarr")
            # d_matrix_zarr_path = kwargs.get('d_matrix_zarr_path', "output/d_matrix.zarr")

            H_matrix = zarr.open_array(H_matrix_zarr_path, mode='r')
            H_local = H_matrix[:, self.nd_start_world:self.nd_end_world]

            #  get the synthetic observations
            synthetic_obs = zarr.open_array(synthetic_obs_zarr_path, mode='r')
            synthetic_obs_local = synthetic_obs[self.nd_start_world:self.nd_end_world, km]
            # print(f"\n synthetic_obs shape: {synthetic_obs.shape} synthetic_obs_local shape: {synthetic_obs_local.shape}\n")

            #  compute the d matrix
            d_local = np.dot(H_local, synthetic_obs_local)

            # compute global d from all ranks
            d_global = np.empty_like(d_local)
            self.mpi_comm.Allreduce(d_local, d_global, op=MPI.SUM)

            # wrapper function for all ensemble members
            def D_HA(ens_idx):
                # call Eta_matrix
                Eta_global = self.Eta_matrix(km, ens_idx, **kwargs)

                # compute D
                D_global = d_global + Eta_global

                # compute HA
                state_local = self.read_analysis(km, ens_idx)
                HA_local = np.dot(H_local, state_local)
                # compute global HA from all ranks
                HA_global = np.empty_like(HA_local)
                self.mpi_comm.Allreduce(HA_local, HA_global, op=MPI.SUM)
                return Eta_global, D_global, HA_global

            # compute innovations for all ens members
            # -- create a zarr file for each processor to write to  
            m = kwargs.get('m_obs')*2+1
            Nens = self.nens
            Dprime = np.zeros((m,self.nens))
            HA = np.zeros((m, self.nens))
            Eta = np.zeros((m, self.nens))
            for ens_idx in range(self.nens):
                Eta_global, D_global, HA_global = D_HA(ens_idx)
                Dprime[:, ens_idx] = D_global - HA_global
                HA[:, ens_idx] = HA_global
                Eta[:, ens_idx] = Eta_global

            # compute the HAbar
            # HAbar = np.mean(HA, axis=1)
            # HAprime_local = HA - HAbar[:, np.newaxis]
            one_N = np.ones((Nens,Nens))/Nens
            HAprime= HA@(np.eye(Nens) - one_N) # mxNens

            # compute HAprime + Eta
            HAprime_Eta = HAprime + Eta
            print(f"[Rank {self.mpi_comm.Get_rank()}] HAprime_Eta_local norm: {np.linalg.norm(HAprime_Eta)}, shape: {HAprime_Eta.shape}")
            # print(f"\n [Rank {self.mpi_comm.Get_rank()}] Dprime_local shape: {Dprime_local.shape} HAprime_local shape: {HAprime_local.shape} HAprime_Eta_local shape: {HAprime_Eta_local.shape}\n ")

            # compute SVD of HAprime_Eta
            U, sig, Vt = np.linalg.svd(HAprime_Eta, full_matrices=False)

            # get the min (m Nens)
            nrmin = min(Nens, m)
            
            # convert S to eigenvalues
            sig = sig**2

            # compute the number of significant eigenvalues
            # sigsum = np.sum(sig[:nrmin]) #computes the total sum of the first nrmin eigenvalues
            # print(f"[Rank {self.mpi_comm.Get_rank()}] sigsum: {sigsum}, sig: {sig[:nrmin]}")
            # sigsum1 = 0.0
            # nrsigma = 0

            # for i in range(nrmin):
            #     if sigsum1 / sigsum < 0.999:
            #         nrsigma += 1
            #         sigsum1 += sig[i]
            #         sig[i] = 1.0 / sig[i]  # Inverse of eigenvalue
            #     else:
            #         sig[i:nrmin] = 0.0  # Set remaining eigenvalues to 0
            #         break  # Exit the loop

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

            #  compute X3 = U*X2 # m_obs x Nens
            X3 = np.dot(U, X2)

            # print(f"[ICESEE] Rank: {rank_world} X3 shape: {X3.shape}")
            # compute X4 = (HAprime.T)*X3 # Nens x Nens
            X4 = np.dot(HAprime.T, X3)
            del X2, X3, U, HAprime, HA, Eta, Dprime
            gc.collect()

            # compute X5 = X4 + I
            X5 = X4 + np.eye(Nens)
            # sum of each column of X5 should be 1
            if np.sum(X5, axis=0).all() != 1.0:
                print(f"[ICESEE] Sum of each X5 column is not 1.0: {np.sum(X5, axis=0)}")
            # print(f"[ICESEE] Rank: {comm_world.Get_rank()} X5 sum: {np.sum(X5, axis=0)}")
            del X4; gc.collect()

            return X5
        
        except Exception as e:
            print(f"Error in compute_X5_analysis_mean: {e}")
            tb_str = "".join(traceback.format_exception(*sys.exc_info()))
            print(f"Traceback details:\n{tb_str}")
            self.mpi_comm.Abort(1)
    
    def compute_X5_utils_batch(self, km, **kwargs):
        """
        Returns:
        Eta   : (m, Nens)
        Dprime: (m, Nens)   # constant across Nens (each column equal)
        HA    : (m, Nens)
        """
        try:
            comm = self.mpi_comm
            rank = comm.Get_rank()

            # ---- Config / inputs
            H_matrix_zarr_path = kwargs.get('H_matrix_zarr_path', "output/H_matrix.zarr")
            synthetic_obs_zarr_path = kwargs.get('synthetic_obs_zarr_path', "output/synthetic_obs.zarr")
            m = kwargs.get('m_obs') * 2 + 1
            Nens = int(self.nens)
            block_size = int(kwargs.get('block_size', max(16, min(64, Nens))))  # tuneable batch size

            # ---- Open H once and slice local columns
            H_matrix = zarr.open_array(H_matrix_zarr_path, mode='r')  # shape (m_total, nd_total) or (m, nd_total) as per your layout
            H_local = H_matrix[:, self.nd_start_world:self.nd_end_world]  # shape (m, local_nd)
            H_local = np.ascontiguousarray(H_local, dtype=np.float64)
            m_local, local_nd = H_local.shape  # expect m_local == m

            # ---- Read local ensemble mean once
            mean_file_path = f"{self.base_path}/{self.file_prefix}_mean.h5"
            with h5py.File(mean_file_path, 'r', driver='mpio', comm=comm) as f:
                ens_mean_local = f['mean'][self.nd_start_world:self.nd_end_world, km]
            ens_mean_local = np.ascontiguousarray(ens_mean_local, dtype=np.float64)  # (local_nd,)

            # ---- Synthetic observations (local slice)
            synthetic_obs = zarr.open_array(synthetic_obs_zarr_path, mode='r')
            synthetic_obs_local = synthetic_obs[self.nd_start_world:self.nd_end_world, km]
            synthetic_obs_local = np.ascontiguousarray(synthetic_obs_local, dtype=np.float64)  # (local_nd,)

            # ---- Fuse d and Hmean into a single GEMM + single Allreduce
            # Build a 2-column local matrix [y_obs, ens_mean]
            V_local = np.empty((local_nd, 2), dtype=np.float64, order='C')
            V_local[:, 0] = synthetic_obs_local
            V_local[:, 1] = ens_mean_local

            Y_local = H_local @ V_local                   # (m, 2)
            Y_global = np.empty_like(Y_local, order='C')  # (m, 2)

            # Single collective for both vectors
            comm.Allreduce([Y_local, MPI.DOUBLE], [Y_global, MPI.DOUBLE], op=MPI.SUM)
            d_global     = Y_global[:, 0]                 # (m,)
            Hmean_global = Y_global[:, 1]                 # (m,)

            # ---- Compute HA for all ensemble members in batches
            HA_global = np.empty((m_local, Nens), dtype=np.float64, order='C')

            for j0 in range(0, Nens, block_size):
                j1 = min(j0 + block_size, Nens)
                B = j1 - j0

                # Load a contiguous local block of states: shape (local_nd, B)
                States_local_blk = np.empty((local_nd, B), dtype=np.float64, order='C')
                for jj, ens_idx in enumerate(range(j0, j1)):
                    States_local_blk[:, jj] = self.read_analysis(km, ens_idx)  # must return local slice (local_nd,)

                # Local GEMM then one Allreduce for this batch
                HA_local_blk = H_local @ States_local_blk             # (m, B)
                HA_global_blk = np.empty_like(HA_local_blk, order='C')
                comm.Allreduce([HA_local_blk, MPI.DOUBLE], [HA_global_blk, MPI.DOUBLE], op=MPI.SUM)

                HA_global[:, j0:j1] = HA_global_blk

            # ---- Eta and D'
            # Eta = HA - Hmean[:, None]
            Eta = HA_global - Hmean_global[:, None]                   # (m, Nens)

            # D' = (d - Hmean) broadcast across columns
            d_minus_Hmean = (d_global - Hmean_global)                 # (m,)
            # Make every column identical, no extra collectives
            Dprime = np.broadcast_to(d_minus_Hmean[:, None], (m_local, Nens)).copy()

            return Eta, Dprime, HA_global

        except Exception as e:
            print(f"Error in compute_X5_utils: {e}")
            tb_str = "".join(traceback.format_exception(*sys.exc_info()))
            print(f"Traceback details:\n{tb_str}")
            self.mpi_comm.Abort(1)

    def compute_X5_utils(self, km, **kwargs):
        # Eta = HA-Hmean where HA = H*state and Hmean = H*mean(state)
        # Dprime[:ens_idx] = d - Hmean
        try:
            H_matrix_zarr_path = kwargs.get('H_matrix_zarr_path', "output/H_matrix.zarr")
            synthetic_obs_zarr_path = kwargs.get('synthetic_obs_zarr_path', "output/synthetic_obs.zarr")
            m = kwargs.get('m_obs') * 2 + 1
            Nens = self.nens

            # --- Open H (read-only) once and slice local columns
            H_matrix = zarr.open_array(H_matrix_zarr_path, mode='r')
            # H has shape (m, nd_total); we take our local column block
            H_local = H_matrix[:, self.nd_start_world:self.nd_end_world]  # (m, local_nd)

            local_nd = self.nd_end_world - self.nd_start_world

            # --- Read ensemble mean ONCE (parallel)
            mean_file_path = f"{self.base_path}/{self.file_prefix}_mean.h5"
            with h5py.File(mean_file_path, 'r', driver='mpio', comm=self.mpi_comm) as f:
                ens_mean_local = f['mean'][self.nd_start_world:self.nd_end_world, km]  # (local_nd,)

            # --- Synthetic obs (once)
            synthetic_obs = zarr.open_array(synthetic_obs_zarr_path, mode='r')
            synthetic_obs_local = synthetic_obs[self.nd_start_world:self.nd_end_world, km]  # (local_nd,)

            # --- d = H * y_obs (one GEMV + one Allreduce)
            d_local = H_local @ synthetic_obs_local              # (m,)
            d_global = np.empty_like(d_local)
            self.mpi_comm.Allreduce(d_local, d_global, op=MPI.SUM)  # 1st collective

            # --- Hmean = H * ens_mean (one GEMV + (no extra open calls))
            Hmean_local = H_local @ ens_mean_local               # (m,)
            Hmean_global = np.empty_like(Hmean_local)
            self.mpi_comm.Allreduce(Hmean_local, Hmean_global, op=MPI.SUM)  # 2nd collective

            # --- Read all ensemble states locally and batch into a matrix
            # Shape: (local_nd, Nens)
            States_local = np.empty((local_nd, Nens), dtype=H_local.dtype, order='C')
            for j in range(Nens):
                States_local[:, j] = self.read_analysis(km, j)  # each returns local slice (local_nd,)

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

            return Dprime, Eta, HA
        except Exception as e:
            print(f"Error in compute_X5_modified: {e}")
            tb_str = "".join(traceback.format_exception(*sys.exc_info()))
            print(f"Traceback details:\n{tb_str}")
            self.mpi_comm.Abort(1)


    def compute_X5_modified(self, km, **kwargs):
        # Eta = HA-Hmean where HA = H*state and Hmean = H*mean(state)
        # Dprime[:ens_idx] = d - Hmean
        try:
            m = kwargs.get('m_obs') * 2 + 1
            Nens = self.nens

            Dprime, Eta, HA = self.compute_X5_utils(km, **kwargs)
            # Dprime, Eta, HA = self.compute_X5_utils_batch(km, **kwargs)

            # compute the HAbar
            # HAbar = np.mean(HA, axis=1)
            # HAprime_local = HA - HAbar[:, np.newaxis]
            one_N = np.ones((Nens,Nens))/Nens
            HAprime= HA@(np.eye(Nens) - one_N) # mxNens

            # compute HAprime + Eta
            HAprime_Eta = HAprime + Eta
            print(f"[Rank {self.mpi_comm.Get_rank()}] HAprime_Eta_local norm: {np.linalg.norm(HAprime_Eta)}, shape: {HAprime_Eta.shape}")
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

            #  compute X3 = U*X2 # m_obs x Nens
            X3 = np.dot(U, X2)

            # print(f"[ICESEE] Rank: {rank_world} X3 shape: {X3.shape}")
            # compute X4 = (HAprime.T)*X3 # Nens x Nens
            X4 = np.dot(HAprime.T, X3)
            del X2, X3, U, HAprime, HA, Eta, Dprime
            gc.collect()

            # compute X5 = X4 + I
            X5 = X4 + np.eye(Nens)
            # sum of each column of X5 should be 1
            if np.sum(X5, axis=0).all() != 1.0:
                print(f"[ICESEE] Sum of each X5 column is not 1.0: {np.sum(X5, axis=0)}")
            # print(f"[ICESEE] Rank: {comm_world.Get_rank()} X5 sum: {np.sum(X5, axis=0)}")
            del X4; gc.collect()

            return X5

        except Exception as e:
            print(f"Error in compute_X5_root_optimized: {e}")
            tb_str = "".join(traceback.format_exception(*sys.exc_info()))
            print(f"Traceback details:\n{tb_str}")
            self.mpi_comm.Abort(1)


    # compute analysis mean
    def compute_analysis_update(self, km, **kwargs):
        # Compute the analysis update for each rank
        try:
            self._ensure_batch(km)
            comm = self.mpi_comm
            rank = comm.Get_rank()
            size = comm.Get_size()

            batch_idx = km - self.current_batch_start

            start = MPI.Wtime()

            # call the compute X5 function
            # X5 = self.compute_X5_(km, **kwargs)
            X5 = self.compute_X5_modified(km, **kwargs)
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

            for i in range(Nens):
                all_states_zarr[:, i] = self.read_analysis(km, i)

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
            min_thickness = 1e-3
            vec_inputs = kwargs.get("vec_inputs", None)
         
            for i, input_var in enumerate(vec_inputs or []):
                if input_var in THICKNESS_VARS:
                    start = i * ndim
                    end = start + ndim
                    analysis_updates[start:end, :] = np.maximum(analysis_updates[start:end, :], min_thickness)

            # Write back all analysis updates
            for j in range(Nens):
                self.write_analysis(km, analysis_updates[:, j], j)

            # compute the anlysis mean and write to h5 file
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
                f['mean'][self.nd_start_world:self.nd_end_world, km] = analysis_mean
            
            # clean up zarr file
            # if os.path.exists(zarr_path):
            #     shutil.rmtree(zarr_path)
            for path in [allstates_sate_zarr_path, mean_params_zarr_path, pertubations_zarr_path, analysis_updates_zarr_path]:
                if os.path.exists(path):
                    shutil.rmtree(path)
            
            self.mpi_comm.Barrier()
            # del all_states_zarr
            gc.collect()
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