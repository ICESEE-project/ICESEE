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
            self._create_batch(0)
            # if self.serial_file_creation:
            #     self._create_batch_serial(0)
            # else:
            #     self._create_batch_parallel(0)
        except Exception as e:
            print(f"Error occurred in __init__: {e}")
            tb_str = "".join(traceback.format_exception(*sys.exc_info()))
            print(f"Traceback details:\n{tb_str}")
            self.mpi_comm.Abort(1)

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
                try:
                    f.flush()
                except Exception:
                    pass
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
                self._create_batch(batch_start)
                # if self.serial_file_creation:
                #     self._create_batch_serial(batch_start)
                # else:
                #     self._create_batch_parallel(batch_start)
        except Exception as e:
            print(f"Error occurred in _ensure_batch: {e}")
            tb_str = "".join(traceback.format_exception(*sys.exc_info()))
            print(f"Traceback details:\n{tb_str}")
            self.mpi_comm.Abort(1)

    @retry_on_failure(max_attempts=3, delay=0.5, mpi_comm=MPI.COMM_WORLD)
    def read_forecast(self, t, ens):
        self._ensure_batch(t)
        batch_idx = t - self.current_batch_start
        try:
            data = self.datasets[batch_idx][self.nd_start:self.nd_end, ens]
            print(f"[ICESEE] Finished reading ensemble {ens} ensemble shape: {data.shape} norm {np.linalg.norm(data)}")
            return data
        except Exception as e:
            print(f"Error occurred in read_forecast: {e}")
            tb_str = "".join(traceback.format_exception(*sys.exc_info()))
            print(f"Traceback details:\n{tb_str}")
            self.mpi_comm.Abort(1)


    @retry_on_failure(max_attempts=3, delay=0.5, mpi_comm=MPI.COMM_WORLD)
    def write_forecast(self, t, data_batch, ens_indices):
        self._ensure_batch(t)
        batch_idx = t - self.current_batch_start
        start = MPI.Wtime()

        # data_batch: (nd_local, len(ens_indices))
        # for i, ens_idx in enumerate(ens_indices):
        #     self.datasets[batch_idx][self.nd_start:self.nd_end, ens_idx] = data_batch[:, i]
        self.datasets[batch_idx][self.nd_start:self.nd_end, ens_indices] = data_batch
        # Flush once per batch/call
        # self.files[batch_idx].flush()
        write_time = MPI.Wtime() - start
      

    def close(self):
        try:
            self._close_batch()
        except Exception as e:
            print(f"Error occurred in close: {e}")
            tb_str = "".join(traceback.format_exception(*sys.exc_info()))
            print(f"Traceback details:\n{tb_str}")
            self.mpi_comm.Abort(1)

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
        # try:
        synthetic_obs_zarr_path = kwargs.get('synthetic_obs_zarr_path')
        error_R_zarr_path = kwargs.get('error_R_zarr_path')
        nd = self.nd
        nt = self.nt

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

        # if rank == 0:
        #     for i, sig in enumerate(self.params["sig_obs"]):
        #         start_idx = i*hdim
        #         end_idx = start_idx + hdim
        #         error_R[start_idx:end_idx,:] = np.ones((hdim,1)) * sig
        # self.mpi_comm.Barrier()

        statevec_true = zarr.open_array(store=f"{self.base_path}/statevec_true.zarr", mode='r+')
        _, indx_map, _ = icesee_get_index(**kwargs)
        if self.nd < 10000:
            if rank==0:
                print("[ICESEE] Generating synthetic observations ...")
                km = 0
                for step in range(nt):
                    if (km<m_obs) and (step+1 == ind_m[km]):
                        for key in kwargs['vec_inputs']:
                            # hu_obs[indx_map[key],km] = statevec_true[indx_map[key],step+1]
                            hu_obs[indx_map[key],km] = statevec_true[indx_map[key],step+1] + np.random.normal(0,error_R[indx_map[key],km],len(indx_map[key]))

                        km += 1
                # print(f"\n nd = {nd}, Nens = {self.nens}, nt = {nt}\n")        
            self.mpi_comm.Barrier()
        else:
            if rank == 0:
                print("[ICESEE] Generating synthetic observations in parallel ...")
                # print(f"\n nd = {nd}, Nens = {self.nens}, nt = {nt}\n")
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
        # except Exception as e:
        #     print(f"Error in _create_synthetic_observations: {e}")
        #     tb_str = "".join(traceback.format_exception(*sys.exc_info()))
        #     print(f"Traceback details:\n{tb_str}")
        #     self.mpi_comm.Abort(1)

    def H_matrix(self, **kwargs):
        try:
            zarr_path = kwargs.get('H_matrix_zarr_path')
            nd = self.nd
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