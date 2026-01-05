import h5py
import numpy as np
from mpi4py import MPI
import os
import glob
import gc
import zarr

class EnKFIO:
    def __init__(self, file_prefix, nd, nens, nt, tobserve, subcomm, mpi_comm, serial_file_creation=True, base_path="enkf_data", batch_size=50):
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
        self.tobserve = tobserve
        self.nt_obs = len(tobserve)
        self.base_path = base_path
        self.file_prefix = file_prefix
        self.batch_size = batch_size
        self.mpi_comm = mpi_comm
        self.comm = subcomm
        self.rank = self.comm.Get_rank()
        self.size = self.comm.Get_size()
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
        if t not in self.tobserve:
            raise ValueError(f"Time step {t} is not an observation time")
        self._ensure_batch(t)
        batch_idx = t - self.current_batch_start
        start = MPI.Wtime()
        data = self.datasets[batch_idx][self.nd_start:self.nd_end, ens_idx]
        read_time = MPI.Wtime() - start
        return data

    def write_analysis(self, t, data, ens_idx):
        if t not in self.tobserve:
            raise ValueError(f"Time step {t} is not an observation time")
        self._ensure_batch(t)
        batch_idx = t - self.current_batch_start
        start = MPI.Wtime()
        self.datasets[batch_idx][self.nd_start:self.nd_end, ens_idx] = data
        write_time = MPI.Wtime() - start

    def write_matrix(self, t, dataset_name, data, ens_idx):
        if t not in self.tobserve:
            raise ValueError(f"Time step {t} is not an observation time")
        self._ensure_batch(t)
        batch_idx = t - self.current_batch_start
        start = MPI.Wtime()
        self.files[batch_idx][dataset_name][self.nd_start:self.nd_end, ens_idx] = data
        write_time = MPI.Wtime() - start

    def read_matrix(self, t, dataset_name, ens_idx):
        if t not in self.tobserve:
            raise ValueError(f"Time step {t} is not an observation time")
        self._ensure_batch(t)
        batch_idx = t - self.current_batch_start
        start = MPI.Wtime()
        data = self.files[batch_idx][dataset_name][self.nd_start:self.nd_end, ens_idx]
        read_time = MPI.Wtime() - start
        return data

    def gather_matrix(self, t, dataset_name):
        if t not in self.tobserve:
            raise ValueError(f"Time step {t} is not an observation time")
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

    def close(self):
        self._close_batch()

