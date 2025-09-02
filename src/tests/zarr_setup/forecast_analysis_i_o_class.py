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

# import utility functions
# from ICESEE.src.utils.tools import icesee_get_index

class EnKFIO:
    def __init__(self, file_prefix, nd, nens, nt, subcomm, mpi_comm, params, serial_file_creation=True, base_path="enkf_data", batch_size=50):
        """
        Initialize EnKF I/O manager for (nd, nens) data with dynamic batch file creation.
        
        Args:
            file_prefix (str): Prefix for HDF5 files (e.g., 'enkf' -> 'enkf_0000.h5').
            nd (int): State dimension.
            nens (int): Number of ensemble members.
            nt (int): Total number of time steps.
            tobserve (list): List of observation time steps (0-based indices).
            mpi_comm: MPI communicator (e.g., MPI.COMM_WORLD).
            base_path (str): Directory for HDF5 files.
            batch_size (int): Number of files per batch for large nt and nd.
        """
        self.nd = nd
        self.nens = nens
        self.nt = nt
        # self.tobserve = tobserve
        self.params = params
        # self.nt_obs = len(tobserve)
        self.base_path = base_path
        self.file_prefix = file_prefix
        self.batch_size = batch_size
        self.mpi_comm = mpi_comm
        self.comm = subcomm if nens >= mpi_comm.Get_size() else mpi_comm
        self.rank = self.comm.Get_rank() if nens >= mpi_comm.Get_size() else mpi_comm.Get_rank()
        self.size = self.comm.Get_size() if nens >= mpi_comm.Get_size() else mpi_comm.Get_size()
        # self.rank = self.mpi_comm.Get_rank()  # Use mpi_comm for domain decomposition
        # self.size = self.mpi_comm.Get_size()
        self.subcomm = subcomm
        self.serial_file_creation = serial_file_creation

        # Divide nd among ranks
        nd_local_base = nd // self.size
        remainder = nd % self.size
        if self.rank < remainder:
            self.nd_local = nd_local_base + 1
            self.nd_start = self.rank * (nd_local_base + 1)
        else:
            self.nd_local = nd_local_base
            self.nd_start = remainder * (nd_local_base + 1) + (self.rank - remainder) * nd_local_base
        self.nd_end = self.nd_start + self.nd_local

        # print(f"Rank {self.rank}/{self.size}: nd_local={self.nd_local}, nd_start={self.nd_start}, nd_end={self.nd_end}")

        # foranalysis lets use mpi_com
        size_world = mpi_comm.Get_size()
        rank_world = mpi_comm.Get_rank()
        nd_local_base_world = nd // size_world
        remainder_world = nd % size_world
        if rank_world < remainder_world:
            self.nd_local_world = nd_local_base_world + 1
            self.nd_start_world = rank_world * (nd_local_base_world + 1)
        else:
            self.nd_local_world = nd_local_base_world
            self.nd_start_world = remainder_world * (nd_local_base_world + 1) + (rank_world - remainder_world) * nd_local_base_world
        self.nd_end_world = self.nd_start_world + self.nd_local_world

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
    
    def _create_batch_serial(self, t_start):
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
                    # row_chunk = self.nd_end_world - self.nd_start_world
                    col_chunk = 1
                    f.create_dataset(
                        'states', (self.nd, self.nens),
                        chunks=(row_chunk, col_chunk),
                        compression="gzip", compression_opts=4,
                        dtype='f8'
                    )
                    
        self.mpi_comm.Barrier()

        for t in range(t_start, t_start + nfiles):
            fname = f"{self.base_path}/{self.file_prefix}_{t:04d}.h5"
            f = h5py.File(fname, 'a', driver='mpio', comm=self.comm)
            dset = f['states']
            self.files.append(f)
            self.datasets.append(dset)

        # if self.rank == 0:
        #     print(f"Batch creation time (t_start={t_start}): {MPI.Wtime() - start:.2f} seconds")

    def _create_batch_parallel(self, t_start):
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
            # row_chunk = self.nd_end_world - self.nd_start_world
            col_chunk = 1
            dset = f.create_dataset(
                'states', (self.nd, self.nens),
                chunks=(row_chunk, col_chunk),
                compression="gzip", compression_opts=4,
                dtype='f8'
            )
            self.files.append(f)
            self.datasets.append(dset)
        # if self.rank == 0:
        #     print(f"Batch creation time (t_start={t_start}): {MPI.Wtime() - start:.2f} seconds")

    def _close_batch(self):
        for f in self.files:
            f.close()
        self.files = []
        self.datasets = []

    def _ensure_batch(self, t):
        batch_start = (t // self.batch_size) * self.batch_size
        if batch_start != self.current_batch_start:
            if self.serial_file_creation:
                self._create_batch_serial(batch_start)
            else:
                self._create_batch_parallel(batch_start)

    def read_forecast(self, t, ens_idx):
        self._ensure_batch(t)
        batch_idx = t - self.current_batch_start
        start = MPI.Wtime()
        data = self.datasets[batch_idx][self.nd_start:self.nd_end, ens_idx]
        read_time = MPI.Wtime() - start
        return data

    def write_forecast(self, t, data, ens_idx):
        self._ensure_batch(t)
        batch_idx = t - self.current_batch_start
        start = MPI.Wtime()
        self.datasets[batch_idx][self.nd_start:self.nd_end, ens_idx] = data
        write_time = MPI.Wtime() - start

    def read_analysis(self, t, ens_idx):
        # if t not in self.tobserve:
        #     raise ValueError(f"Time step {t} is not an observation time")
        self._ensure_batch(t)
        batch_idx = t - self.current_batch_start
        start = MPI.Wtime()
        # data = self.datasets[batch_idx][self.nd_start:self.nd_end, ens_idx]
        data = self.datasets[batch_idx][self.nd_start_world:self.nd_end_world, ens_idx]
        read_time = MPI.Wtime() - start
        return data

    def write_analysis(self, t, data, ens_idx):
        # if t not in self.tobserve:
        #     raise ValueError(f"Time step {t} is not an observation time")
        self._ensure_batch(t)
        batch_idx = t - self.current_batch_start
        start = MPI.Wtime()
        self.datasets[batch_idx][self.nd_start:self.nd_end, ens_idx] = data
        write_time = MPI.Wtime() - start

    def write_matrix(self, t, dataset_name, data, ens_idx):
        # if t not in self.tobserve:
        #     raise ValueError(f"Time step {t} is not an observation time")
        self._ensure_batch(t)
        batch_idx = t - self.current_batch_start
        start = MPI.Wtime()
        self.files[batch_idx][dataset_name][self.nd_start:self.nd_end, ens_idx] = data
        write_time = MPI.Wtime() - start

    def read_matrix(self, t, dataset_name, ens_idx):
        # if t not in self.tobserve:
        #     raise ValueError(f"Time step {t} is not an observation time")
        self._ensure_batch(t)
        batch_idx = t - self.current_batch_start
        start = MPI.Wtime()
        data = self.files[batch_idx][dataset_name][self.nd_start:self.nd_end, ens_idx]
        read_time = MPI.Wtime() - start
        return data

    def gather_matrix(self, t, dataset_name):
        # if t not in self.tobserve:
        #     raise ValueError(f"Time step {t} is not an observation time")
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

    def compute_forecast_mean(self, t):
        self._ensure_batch(t)
        comm = self.mpi_comm
        rank = comm.Get_rank()
        size = comm.Get_size()

        batch_idx = t - self.current_batch_start
        start = MPI.Wtime()

        self.mpi_comm.barrier()
        # Each rank holds a disjoint block of rows in [nd_start_world:nd_end_world)
        local_data = self.datasets[batch_idx][self.nd_start_world:self.nd_end_world, :]
        # Per-row mean across columns for the local rows
        local_mean = np.mean(local_data, axis=1).astype('f8', copy=False)  # shape: (nd_local,)

        # Gather the number of rows each rank is processing
        local_row_count = np.array([local_mean.shape[0]], dtype='i8')
        global_row_counts = np.zeros(size, dtype='i8') if rank == 0 else None
        comm.Gather(local_row_count, global_row_counts, root=0)

        # Compute displacements for Gatherv
        if rank == 0:
            displacements = np.zeros(size, dtype='i8')
            displacements[1:] = np.cumsum(global_row_counts[:-1])
            total_rows = np.sum(global_row_counts)
            global_mean = np.zeros(total_rows, dtype='f8')
        else:
            displacements = None
            global_mean = None

        # Use Gatherv to collect local means into global_mean on root
        comm.Gatherv(local_mean, [global_mean, global_row_counts, displacements, MPI.DOUBLE], root=0)

        # Broadcast the global mean to all ranks
        if rank == 0:
            result = global_mean
            
            # write result h5py file for each time step
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

    def compute_forecast_mean_chunked_gather(self, t, ens_chunk_size=1):
        self._ensure_batch(t)
        comm = self.mpi_comm
        rank = comm.Get_rank()
        size = comm.Get_size()

        batch_idx = t - self.current_batch_start
        start = MPI.Wtime()

        # comm.barrier()
        # Each rank holds a disjoint block of rows in [nd_start_world:nd_end_world)
        # Initialize local sum for accumulating partial sums across ensemble chunks
        local_rows = self.nd_end_world - self.nd_start_world
        local_sum = np.zeros(local_rows, dtype='f8')

        # Process ensembles in contiguous chunks to optimize I/O for chunked HDF5 datasets
        for start_ens in range(0, self.nens, ens_chunk_size):
            end_ens = min(start_ens + ens_chunk_size, self.nens)
            # Load a chunk of columns (ensembles) for local rows
            local_data = self.datasets[batch_idx][self.nd_start_world:self.nd_end_world, start_ens:end_ens]
            # Compute sum across the chunk's columns
            partial_sum = np.sum(local_data, axis=1).astype('f8', copy=False)
            local_sum += partial_sum

        # Compute the local mean after processing all chunks
        local_mean = local_sum / self.nens

        # Gather the number of rows each rank is processing
        local_row_count = np.array([local_mean.shape[0]], dtype='i8')
        global_row_counts = np.zeros(size, dtype='i8') if rank == 0 else None
        comm.Gather(local_row_count, global_row_counts, root=0)

        # Compute displacements for Gatherv
        if rank == 0:
            displacements = np.zeros(size, dtype='i8')
            displacements[1:] = np.cumsum(global_row_counts[:-1])
            total_rows = np.sum(global_row_counts)
            global_mean = np.zeros(total_rows, dtype='f8')
        else:
            displacements = None
            global_mean = None

        # Use Gatherv to collect local means into global_mean on root
        comm.Gatherv(local_mean, [global_mean, global_row_counts, displacements, MPI.DOUBLE], root=0)

        # Write result to h5py file for each time step on root
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
        # comm.barrier()


    def compute_forecast_mean_chunked(self, t, ens_chunk_size=None, use_collective_io=False, max_ranks=4):
        self._ensure_batch(t)
        comm = self.mpi_comm
        rank = comm.Get_rank()
        size = comm.Get_size()

        import numpy as np
        import h5py.h5p as h5p
        import h5py.h5s as h5s
        import h5py.h5fd as h5fd

        # Create sub-communicator for active ranks
        active_ranks = min(max_ranks, size)  # Use at most max_ranks processes
        color = 0 if rank < active_ranks else MPI.UNDEFINED
        sub_comm = comm.Split(color, rank)
        active = sub_comm != MPI.COMM_NULL

        if active:
            sub_rank = sub_comm.Get_rank()
            sub_size = sub_comm.Get_size()

            # Distribute rows among active ranks
            local_rows = (self.nd_end_world - self.nd_start_world) if rank == 0 else 0
            local_rows = sub_comm.bcast(local_rows, root=0)
            rows_per_rank = local_rows // sub_size
            extra_rows = local_rows % sub_size
            nd_start = sub_rank * rows_per_rank + min(sub_rank, extra_rows)
            nd_end = nd_start + rows_per_rank + (1 if sub_rank < extra_rows else 0)

            batch_idx = t - self.current_batch_start
            ds = self.datasets[batch_idx]  # Input dataset (nd, nens)

            # Dynamic chunk size
            if ens_chunk_size is None:
                bytes_per_element = 8  # float64
                target_memory = 1e9  # 1 GB
                ens_chunk_size = max(1, int(target_memory / (nd_end - nd_start) / bytes_per_element))
                ens_chunk_size = min(ens_chunk_size, self.nens)
                if sub_rank == 0:
                    print(f"Dynamic ens_chunk_size: {ens_chunk_size}")

            local_sum = np.zeros(nd_end - nd_start, dtype='f8')

            # I/O property list
            dxpl = h5p.create(h5p.DATASET_XFER)
            if use_collective_io:
                dxpl.set_dxpl_mpio(h5fd.MPIO_COLLECTIVE)
            else:
                dxpl.set_dxpl_mpio(h5fd.MPIO_INDEPENDENT)

            # Timing for diagnostics
            t_io = 0.0
            t_comp = 0.0
            start_total = MPI.Wtime()

            # Double buffering
            buffers = [np.empty((nd_end - nd_start, ens_chunk_size), dtype='f8') for _ in range(2)]
            current_buffer = 0

            for start_ens in range(0, self.nens, ens_chunk_size):
                end_ens = min(start_ens + ens_chunk_size, self.nens)
                chunk_cols = end_ens - start_ens

                if chunk_cols < ens_chunk_size:
                    buffers[current_buffer] = np.empty((nd_end - nd_start, chunk_cols), dtype='f8')

                # Read data
                t_start_io = MPI.Wtime()
                file_space = ds.id.get_space()
                file_space.select_hyperslab((self.nd_start_world + nd_start, start_ens), (nd_end - nd_start, chunk_cols))
                mem_space = h5s.create_simple((nd_end - nd_start, chunk_cols))
                ds.id.read(mem_space, file_space, buffers[current_buffer], dxpl=dxpl)
                t_io += MPI.Wtime() - t_start_io

                # Compute sum for previous buffer
                if start_ens > 0:
                    t_start_comp = MPI.Wtime()
                    local_sum += np.sum(buffers[1 - current_buffer], axis=1)
                    t_comp += MPI.Wtime() - t_start_comp

                current_buffer = 1 - current_buffer

            # Compute sum for last buffer
            t_start_comp = MPI.Wtime()
            local_sum += np.sum(buffers[current_buffer], axis=1)
            t_comp += MPI.Wtime() - t_start_comp

            # Compute local mean
            local_mean = local_sum / self.nens

            # Output file handling
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

            # Diagnostics
            if sub_rank == 0:
                print(f"Total time: {MPI.Wtime() - start_total:.2f}s, I/O: {t_io:.2f}s, Compute: {t_comp:.2f}s")

        # Ensure all ranks synchronize before returning
        comm.Barrier()
        # sub_comm.Free()

    def _compute_forecast_mean_chunked(self, t, ens_chunk_size=1):
        self._ensure_batch(t)
        comm = self.mpi_comm
        rank = comm.Get_rank()
        size = comm.Get_size()

        batch_idx = t - self.current_batch_start
        start = MPI.Wtime()

        # Each rank holds a disjoint block of rows in [nd_start_world:nd_end_world)
        # Initialize local sum for accumulating partial sums across ensemble chunks
        local_rows = self.nd_end_world - self.nd_start_world
        local_sum = np.zeros(local_rows, dtype='f8')

        # Process ensembles in contiguous chunks to optimize I/O for chunked HDF5 datasets
        for start_ens in range(0, self.nens, ens_chunk_size):
            end_ens = min(start_ens + ens_chunk_size, self.nens)
            # Load a chunk of columns (ensembles) for local rows
            local_data = self.datasets[batch_idx][self.nd_start_world:self.nd_end_world, start_ens:end_ens]
            # Compute sum across the chunk's columns
            partial_sum = np.sum(local_data, axis=1).astype('f8', copy=False)
            local_sum += partial_sum

        # Compute the local mean after processing all chunks
        local_mean = local_sum / self.nens

        # Each rank writes its portion of the mean to the HDF5 file
        file_path = f"{self.base_path}/{self.file_prefix}_mean.h5"
        with h5py.File(file_path, 'a', driver='mpio', comm=comm) as f:
            if 'mean' not in f:
                # Only one rank (e.g., rank 0) creates the dataset to avoid race conditions
                if rank == 0:
                    f.create_dataset(
                        'mean', (self.nd, self.nt),
                        chunks=(min(self.nd, 1000), 1),
                        dtype='f8'
                    )
            # Ensure dataset is created before all ranks write
            comm.Barrier()

            # Each rank writes its local_mean to its portion of the dataset
            # f['mean'][self.nd_start_world:self.nd_end_world, t] = local_mean

        return
    
    def compare_forecast_means(self, t, ens_chunk_size=1, rtol=1e-5, atol=1e-8):
        """
        Compare the means computed by compute_forecast_mean and compute_forecast_mean_chunked.
        
        Parameters:
        - t: Time step to compute the mean for.
        - ens_chunk_size: Chunk size for the chunked method.
        - rtol: Relative tolerance for comparison.
        - atol: Absolute tolerance for comparison.
        
        Returns:
        - bool: True if means are equivalent, False otherwise (only on root rank).
        """
        comm = self.mpi_comm
        rank = comm.Get_rank()

        # Compute mean using the original method
        # self.compute_forecast_mean(t)
        mean_original = None
        if rank == 0:
            file_path = f"{self.base_path}/{self.file_prefix}_mean.h5"
            # with h5py.File(file_path, 'r') as f:
            #     mean_original = f['mean'][:, t].copy()
            self._ensure_batch(t)
            batch_idx = t - self.current_batch_start
            ens_mean = self.datasets[batch_idx][:,:]
            mean_original = np.mean(ens_mean, axis=1)

        # Compute mean using the chunked method
        self.compute_forecast_mean_chunked(t, ens_chunk_size=ens_chunk_size)
        mean_chunked = None
        if rank == 0:
            with h5py.File(file_path, 'r') as f:
                mean_chunked = f['mean'][:, t].copy()

        # Compare on root rank
        if rank == 0:
            # Ensure both arrays are not None
            if mean_original is None or mean_chunked is None:
                print("Error: One or both means could not be read from file.")
                return False

            # Compare using numpy's isclose for floating-point comparison
            are_equal = np.allclose(mean_original, mean_chunked, rtol=rtol, atol=atol)
            if are_equal:
                print(f"Means for time step {t} are equivalent within tolerance (rtol={rtol}, atol={atol}).")
            else:
                print(f"Means for time step {t} differ beyond tolerance (rtol={rtol}, atol={atol}).")
                # Optionally print max difference
                max_diff = np.max(np.abs(mean_original - mean_chunked))
                print(f"Maximum absolute difference: {max_diff}")
            return are_equal
        return None  # Non-root ranks return None

    def icesee_get_index(self, **kwargs):
        """
        Get the index of the variables in the vector dynamically.
        #TODO remove the loading of vec for memory optimizations

        Parameters:
            - vec: The vector to be distributed
            - vec_inputs: List of input variable names, e.g., ['h', 'u', 'v', 'smb']
            - hdim: Size of each variable in vec_inputs
            - dim_list: List of sizes of each rank; ordered through the ranks using MPI_gather on root and broadcast
            - comm: MPI communicator containing the rank and size of the processors

        Returns:
            - var_indices: Dictionary where keys are variable names and values are their respective slices from `vec`
            - index_map: Dictionary where keys are variable names and values are the indices corresponding to their slices
        """
        # -- get parameters
        vec_inputs = kwargs.get("vec_inputs", None)
        params = kwargs.get("params", None)
        nd = kwargs.get("nd", params.get("nd", None))

        if params["default_run"]:
            comm = kwargs.get("subcomm", None)
        else:
            comm = kwargs.get("comm_world", None)
        
        # get size of input vector based on user inputs
        len_vec = params["total_state_param_vars"]

        # print(f"[ICESEE] dim_list: {kwargs['dim_list']}")
        dim_list_param = np.array(kwargs.get('dim_list', None)) // len_vec  # Get the size of each variable slice
        # hdim = vec.shape[0] // len_vec  # Compute the size of each variable in vec_inputs
        hdim = nd // len_vec  # Compute the size of each variable in vec_inputs

        if comm is None:
            # Non-MPI case
            rank = 0
            dim = dim_list_param[rank]
            offsets = [0]  # No offsets needed
        else:
            # MPI case
            size_world = kwargs.get("comm_world").Get_size()  # Get the total number of processors
            # if params["even_distribution"] or (params["default_run"] and size_world <= params["Nens"]):
            if params["even_distribution"]:
                rank = 0 # Set rank to 0 for even distribution
                dim = dim_list_param[rank]
                offsets = [0]
            else:
                rank = comm.Get_rank()  # Get the rank of the current processor
                dim = dim_list_param[rank]
                offsets = np.cumsum(np.insert(dim_list_param, 0, 0))  # Compute offsets per processor

        start_idx = offsets[rank]  # Get the start index of the current processor
    
        # Dynamically determine start indices for each variable
        # var_indices = {}
        index_map = {}
        var_start = 0  # Initial start index

        for var in vec_inputs:
            start = var_start + start_idx
            end = start + dim
            # var_indices[var] = vec[start:end]
            index_map[var] = np.arange(start, end)  # Store index range for easy fetching
            var_start += hdim  # Move to the next variable slice

        local_size_per_rank = kwargs.get('dim_list', None)
        return index_map, local_size_per_rank[rank]

    def generate_observation_schedule_(self,**kwargs):
        """
        Generate observation times and indices from a given array of time points.

        Parameters:
            t (list or np.ndarray): Array of time points.
            freq_obs (int): Frequency of observations in the same unit as `t`.
            obs_max_time (int): Maximum observation time in the same unit as `t`.

        Returns:
            obs_t (list): Observation times.
            obs_idx (list): Indices corresponding to observation times in `t`.
        """
        # unpack kwargs
        t = kwargs["t"]

        # Convert input to a numpy array for easier manipulation
        t = np.array(t)
        
        # Generate observation times
        obs_t = np.arange(self.params["obs_start_time"], self.params["obs_max_time"] + self.params["freq_obs"], self.params["freq_obs"])
        print(f"Observation times (obs_t): {obs_t} for obs_start_time = {self.params['obs_start_time']}, obs_max_time = {self.params['obs_max_time']}, freq_obs = {self.params['freq_obs']}")
        # obs_t = np.linspace(obs_start_time, obs_max_time, int(obs_max_time/freq_obs)+1)
        
        # Find indices of observation times in the original array
        obs_idx = np.array([np.where(t == time)[0][0] for time in obs_t if time in t]).astype(int)
        print(f"Observation indices (obs_idx): {obs_idx} for obs_t = {obs_t} in t = {t}")

        # print(f"Number of observation instants: {len(obs_idx)} at times: {t[obs_idx]} for t = {t}")

        # number of observation instants
        num_observations = len(obs_idx)

        return obs_t, obs_idx, num_observations

    def generate_observation_schedule(self, **kwargs):
        """
        Generate observation times and indices from a given array of time points.

        Parameters:
            t (list or np.ndarray): Array of time points.
            freq_obs (int): Frequency of observations in the same unit as `t`.
            obs_start_time (int): Start time for observations in the same unit as `t`.
            obs_max_time (int): Maximum observation time in the same unit as `t`.

        Returns:
            obs_t (list): Observation times.
            obs_idx (list): Indices corresponding to observation times in `t`.
            num_observations (int): Number of observation instants.
        """
        # Unpack kwargs
        t = np.array(kwargs["t"])  # Convert to numpy array
        freq_obs = self.params["freq_obs"]
        obs_start_time = self.params["obs_start_time"]
        obs_max_time = self.params["obs_max_time"]

        # Cap obs_max_time to the maximum value in t
        max_t = np.max(t)
        obs_max_time = min(obs_max_time, max_t)

        # Generate observation times, ensuring not to exceed obs_max_time
        obs_t = np.arange(obs_start_time, obs_max_time + freq_obs, freq_obs)
        obs_t = obs_t[obs_t <= obs_max_time]  # Filter to exclude times > obs_max_time
        # print(f"Observation times (obs_t): {obs_t} for obs_start_time = {obs_start_time}, obs_max_time = {obs_max_time}, freq_obs = {freq_obs}")

        # Find indices of the closest times in t
        obs_idx = []
        for time in obs_t:
            # Find the index of the closest time in t
            idx = np.argmin(np.abs(t - time))
            # Include the index regardless of exact match, as long as it's the closest
            obs_idx.append(idx)
        obs_idx = np.array(obs_idx, dtype=int)
        # print(f"Observation indices (obs_idx): {obs_idx} for obs_t = {obs_t} in t = {t}")
        # print(f"Matched times in t: {t[obs_idx]}")

        # Number of observation instants
        num_observations = len(obs_idx)

        return obs_t, obs_idx, num_observations

    def generate_observation_schedule_(self, **kwargs):
        """
        Generate observation indices aligned to the provided time grid t.

        Behavior:
        - Assumes t is (approximately) uniformly spaced.
        - First observation is the first grid time >= obs_start_time (ceil-to-grid).
        - Subsequent observations use a fixed index stride ~ freq_obs/dt (rounded).
        - Stops at obs_max_time (inclusive, with small tolerance).

        Returns:
        obs_t   : 1D np.ndarray of observation times (taken from t)
        obs_idx : 1D np.ndarray of integer indices into t
        num_observations : int
        """
        # unpack
        t = np.asarray(kwargs["t"])
        params = self.params  # or kwargs if you pass params there

        obs_start_time = float(params["obs_start_time"])
        obs_max_time   = float(params["obs_max_time"])
        freq_obs       = float(params["freq_obs"])

        # basic validations
        if t.ndim != 1 or t.size < 2:
            raise ValueError("t must be a 1D array with at least two points.")
        if freq_obs <= 0:
            raise ValueError("freq_obs must be > 0.")
        if obs_start_time > t[-1]:
            # nothing to schedule
            return np.array([], dtype=t.dtype), np.array([], dtype=int), 0

        # estimate dt and confirm near-uniform grid
        dt_est = float(np.mean(np.diff(t)))
        if dt_est <= 0:
            raise ValueError("t must be strictly increasing.")
        # tolerance for float comparisons
        atol = max(1e-10, 1e-6 * max(1.0, abs(t[-1])))

        # find first index at or after obs_start_time (ceil-to-grid)
        i0 = int(np.searchsorted(t, obs_start_time, side="left"))
        if i0 >= t.size:
            return np.array([], dtype=t.dtype), np.array([], dtype=int), 0

        # convert freq to an index stride (round to nearest integer step)
        step_idx = max(1, int(round(freq_obs / dt_est)))

        # build candidate indices from i0 with stride step_idx
        # stop when time would exceed obs_max_time (with small tolerance)
        idx = []
        i = i0
        while i < t.size and (t[i] <= obs_max_time + atol):
            idx.append(i)
            i += step_idx

        obs_idx = np.array(idx, dtype=int)
        obs_t = t[obs_idx]
        num_observations = int(obs_idx.size)

        # Optional: log what happened for debugging
        print(f"[schedule] dt≈{dt_est:.6g}, start_i={i0}, step_idx={step_idx}, "
              f"obs_t={obs_t}, obs_idx={obs_idx}")

        return obs_t, obs_idx, num_observations

    def _create_synthetic_observations(self, **kwargs):

        synthetic_obs_zarr_path = kwargs.get('synthetic_obs_zarr_path')
        error_R_zarr_path = kwargs.get('error_R_zarr_path')
        nd = kwargs.get('nd')
        nt = kwargs.get('nt')

        #  get the observation times and indices
        obs_t, ind_m, m_obs = self.generate_observation_schedule(**kwargs)
        # print(f"Number of observation instants: {m_obs} at times: {obs_t}, in indices: {ind_m} ind_m_type: {ind_m.dtype}")

        m = m_obs
        m_R = m_obs*2 +1

        rank = self.mpi_comm.Get_rank()
        size = self.mpi_comm.Get_size()
        try:
            # Clear existing Zarr stores to avoid metadata conflicts (only rank 0)
            if rank == 0:
                if os.path.exists(synthetic_obs_zarr_path):
                    shutil.rmtree(synthetic_obs_zarr_path)
                if os.path.exists(error_R_zarr_path):
                    shutil.rmtree(error_R_zarr_path)
            self.mpi_comm.Barrier()

            if rank == 0:
                hu_obs = zarr.create_array(store=synthetic_obs_zarr_path, shape=(nd, m), chunks=(min(1000, nd), min(50, m)), dtype='f8', overwrite=True)
                error_R = zarr.create_array(store=error_R_zarr_path, shape=(nd, m_R), chunks=(min(1000, nd), min(50, m_R)), dtype='f8', overwrite=True)
            
            # else: 
            #     hu_obs = zarr.open_array(store=synthetic_obs_zarr_path, mode='r+')
            #     error_R = zarr.open_array(store=error_R_zarr_path, mode='r+')

            # self.mpi_comm.Barrier()

            self.mpi_comm.Barrier()
            hu_obs = zarr.open_array(store=synthetic_obs_zarr_path, mode='r+')
            error_R = zarr.open_array(store=error_R_zarr_path, mode='r+')
            self.mpi_comm.Barrier()

            if kwargs.get('joint_estimation', False) or self.params.get('localization_flag', False):
                hdim = nd // self.params["total_state_param_vars"]
            else:
                hdim = nd // self.params["total_state_param_vars"]

            # compute error_R for the other four serial variants
            if rank == 0:
                for i, sig in enumerate(self.params["sig_obs"]):
                    start_idx = i*hdim
                    end_idx = start_idx + hdim
                    error_R[start_idx:end_idx,:] = np.ones((hdim,1)) * sig
            self.mpi_comm.Barrier() 

            # compute synthetic observations
            # open statevec_true Zarr array
            statevec_true = zarr.open_array(store="output/statevec_true.zarr", mode='r+')

            indx_map, _ = self.icesee_get_index(**kwargs)

            # if rank==0:
            #     km = 0
            #     for step in range(nt):
            #         if (km<m_obs) and (step+1 == ind_m[km]):
            #             for key in kwargs['vec_inputs']:
            #                 # hu_obs[indx_map[key],km] = statevec_true[indx_map[key],step+1]
            #                 hu_obs[indx_map[key],km] = statevec_true[indx_map[key],step+1] + np.random.normal(0,error_R[indx_map[key],km],len(indx_map[key]))

            #             km += 1

           
            # Distribute observation instants across processes
            obs_per_process = m_obs // size
            remainder = m_obs % size
            start_obs = rank * obs_per_process + min(rank, remainder)
            num_obs = obs_per_process + 1 if rank < remainder else obs_per_process

            # Partition indices to avoid row chunk conflicts
            rows_per_process = hdim // size
            row_remainder = hdim % size
            row_start = rank * rows_per_process + min(rank, row_remainder)
            row_end = row_start + (rows_per_process + 1 if rank < row_remainder else rows_per_process)

            # Set unique random seed
            # np.random.seed(1234 + rank)

            # Compute synthetic observations and write to hu_obs
            for km in range(start_obs, start_obs + num_obs):
                if km < m_obs:
                    step = ind_m[km] - 1
                    if 0 <= step < nt:
                        for key in kwargs['vec_inputs']:
                            indices = indx_map[key]
                            # Filter indices to local row range
                            local_indices = indices[(indices >= row_start) & (indices < row_end)]
                            if len(local_indices) > 0:
                                state_data = statevec_true[local_indices, step]
                                error_data = error_R[local_indices, km]
                                result = state_data + np.random.normal(0, error_data, len(local_indices))
                                if result.shape != (len(local_indices),):
                                    raise ValueError(f"Rank {rank}: Shape mismatch at km={km}: expected {len(local_indices)}, got {result.shape}")
                                # print(f"Rank {rank}: Writing to hu_obs[{local_indices}, {km}]")
                                hu_obs[local_indices, km] = result
            # Synchronize to ensure all writes are complete
            self.mpi_comm.Barrier()

        except Exception as e:
            # if self.mpi_comm.Get_rank() == 0:
            print(f"Error in _create_synthetic_observations: {e}")
            tb_str = "".join(traceback.format_exception(*sys.exc_info()))
            print(f"Traceback details:\n{tb_str}")
            # self.mpi_comm.Abort(1)


        return obs_t, m_obs

    # observation operator
    def H_matrix(self,**kwargs):
        try:
            zarr_path = kwargs.get('H_matrix_zarr_path')
            nd = kwargs.get('nd')
            m_obs = kwargs.get('m_obs')
            m = m_obs * 2 + 1
            # Calculate distance between measurements
            print(f"Creating H matrix at {zarr_path} with shape m = {m}, nd = {nd} and m_obs = {m_obs}")
            di = int((nd - 2) / (2 *m_obs))

            # Fill the H matrix
            H_matrix_file = zarr.create_array(store=zarr_path, shape=(m,nd), chunks=(min(50,m),min(1000,nd)), dtype='f8', overwrite=True)
            for i in range(1,m_obs + 1):
                H_matrix_file[i - 1, i * di - 1] = 1
                H_matrix_file[m_obs + i - 1, int((nd - 2) / 2) + i * di - 1] = 1

            H_matrix_file[m_obs * 2, nd - 2] = 1  # Final element

            # Check if we have parameter estimation
            if self.params.get('joint_estimation', False):
                ndim = nd // self.params["total_state_param_vars"]
                state_variables_size = ndim * self.params["num_state_vars"]
                num_params_size = nd - state_variables_size
                # H_param = np.zeros(num_params_size)
                # H_matrix_file[:, state_variables_size:] = H_param
                H_matrix_file[:, state_variables_size:] = 0
        except Exception as e:
            print(f"Error in H_matrix: {e}")
            tb_str = "".join(traceback.format_exception(*sys.exc_info()))
            print(f"Traceback details:\n{tb_str}")
            self.mpi_comm.Abort(1)

    # Ensemble pertubations
    def Eta_matrix(self, t, zarr_path="output/H_matrix.zarr"):
     
        try:
            # open H from zarr file
            H_matrix_file = zarr.open_array(zarr_path, mode='r')
            H_local = H_matrix_file[:,self.nd_start_world:self.nd_end_world]

            #  read the forecast mean from h5py file
            mean_file_path = f"{self.base_path}/{self.file_prefix}_mean.h5"
            with h5py.File(mean_file_path,'r', driver='mpio', comm=self.mpi_comm) as f:
                ens_mean = f['mean'][self.nd_start_world:self.nd_end_world, t]
                
            # get local ensemble perturbations
            ens_idx=1
            state = self.read_analysis(t,ens_idx)
            ens_pertubations = state - ens_mean
            Eta_local = np.dot(H_local, ens_pertubations) # (m, nens)
        except Exception as e:
            # if self.mpi_comm.Get_rank() == 0:
            print(f"Error in Eta_matrix: {e}")
            tb_str = "".join(traceback.format_exception(*sys.exc_info()))
            print(f"Traceback details:\n{tb_str}")
            # self.mpi_comm.Abort(1)

        # print(f"\nens_pertubations shape = {ens_pertubations.shape}, Eta_local shape = {Eta_local.shape}  H_local shape = {H_local.shape}\n")
    
    # compute innivations
    def d_matrix(self, t, zarr_path="output/d_matrix.zarr"):
        # d = H @ synthetic_obs
        try:
            # open H from zarr file
            H_matrix_file = zarr.open_array("output/H_matrix.zarr", mode='r')
            H_local = H_matrix_file[:,self.nd_start_world:self.nd_end_world]

            #  read the forecast mean from h5py file
            mean_file_path = f"{self.base_path}/{self.file_prefix}_mean.h5"
            with h5py.File(mean_file_path, driver='mpio', comm=self.mpi_comm) as f:
                ens_mean = f['mean'][self.nd_start_world:self.nd_end_world, t]
                
            # get local ensemble perturbations
            ens_idx=1
            state = self.read_analysis(t,ens_idx)
            ens_pertubations = state - ens_mean
            Eta_local = np.dot(H_local, ens_pertubations) # (m, nens)

            # read synthetic observations from zarr file
            hu_obs = zarr.open_array(store="output/synthetic_obs.zarr", mode='r')
            hu_obs_local = hu_obs[self.nd_start_world:self.nd_end_world,:]

            # compute d matrix
            d_local = hu_obs_local - np.dot(H_local, ens_mean)[:, np.newaxis]  # (m, 1)
            
            # write d matrix to zarr file
            d_file = zarr.open_array(store=zarr_path, mode='a')
            d_file[self.nd_start_world:self.nd_end_world,:] = d_local

        except Exception as e:
            if self.mpi_comm.Get_rank() == 0:
                print(f"Error in d_matrix: {e}")
                tb_str = "".join(traceback.format_exception(*sys.exc_info()))
                print(f"Traceback details:\n{tb_str}")
            self.mpi_comm.Abort(1)

        # print(f"\nens_pertubations shape = {ens_pertubations.shape}, Eta_local shape = {Eta_local.shape}  H_local shape = {H_local.shape}\n")




    def close(self):
        self._close_batch()

