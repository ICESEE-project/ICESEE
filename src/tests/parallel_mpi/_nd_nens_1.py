import h5py
import numpy as np
from mpi4py import MPI
import os
import glob

class EnKFIO:
    def __init__(self, file_prefix, nd, nens, nt, tobserve, subcomm, mpi_comm, base_path="enkf_data", batch_size=50):
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
        # self.batch_size = batch_size if (nt > 500 and nd > 10000) else nt
        self.batch_size = batch_size
        # self.comm = mpi_comm
        self.comm = subcomm
        self.rank = self.comm.Get_rank()
        self.size = self.comm.Get_size()
        self.subcomm = subcomm

        # Partition nd across processes
        # self.nd_local = nd // self.size
        # self.nd_start = self.rank * self.nd_local
        # self.nd_end = self.nd_start + self.nd_local

        # self.nd_local = nd // self.subcomm.Get_size()
        # self.nd_start = self.subcomm.Get_rank() * self.nd_local
        # self.nd_end = self.nd_start + self.nd_local

        # Divide nd among ranks
        nd_local_base = nd // self.size
        remainder = nd % self.size

        # Assign extra row to first `remainder` ranks
        if self.rank < remainder:
            self.nd_local = nd_local_base + 1
            self.nd_start = self.rank * (nd_local_base + 1)
        else:
            self.nd_local = nd_local_base
            self.nd_start = remainder * (nd_local_base + 1) + (self.rank - remainder) * nd_local_base

        self.nd_end = self.nd_start + self.nd_local

        # Create directory and clean up old files
        # if self.rank == 0:
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
        self._create_batch(0)

    def _create_batch(self, t_start):
        start = MPI.Wtime()
        self._close_batch()
        self.files = []
        self.datasets = []
        self.current_batch_start = t_start
        nfiles = min(self.batch_size, self.nt - t_start)
        for t in range(t_start, t_start + nfiles):
            fname = f"{self.base_path}/{self.file_prefix}_{t:04d}.h5"
            f = h5py.File(fname, 'w', driver='mpio', comm=self.comm)
            # dset = f.create_dataset(
            #     'state', (self.nens, self.nd),
            #     chunks=(self.nens, self.nd_local), dtype='f8'
            # )
            # dset = f.create_dataset(
            #     'states', (self.nd, self.nens),
            #     chunks=(self.nd_local, self.nens), dtype='f8'
            # )
            row_chunck = min(1024, self.nd)  # or some small row block
            col_chunk = 1  # or some small column block
            dset = f.create_dataset(
                'states', (self.nd, self.nens),
                chunks=(row_chunck, col_chunk),
                compression="gzip", compression_opts=4,
                dtype='f8'
            )
            self.files.append(f)
            self.datasets.append(dset)
        if self.rank == 0:
            print(f"Batch creation time (t_start={t_start}): {MPI.Wtime() - start:.2f} seconds")

    def _close_batch(self):
        for f in self.files:
            f.close()
        self.files = []
        self.datasets = []

    def _ensure_batch(self, t):
        batch_start = (t // self.batch_size) * self.batch_size
        if batch_start != self.current_batch_start:
            self._create_batch(batch_start)

    def read_forecast(self, t,ens_idx):
        self._ensure_batch(t)
        batch_idx = t - self.current_batch_start
        start = MPI.Wtime()
        # data = self.datasets[batch_idx][:, self.nd_start:self.nd_end]
        # data = self.datasets[batch_idx][self.nd_start:self.nd_end, :]
        data = self.datasets[batch_idx][self.nd_start:self.nd_end, ens_idx]
        # data = self.datasets[batch_idx][:, ens_idx]
        read_time = MPI.Wtime() - start
        if self.rank == 0:
            print(f"t={t}, Read: {read_time:.2f}s")
        # return data.T
        return data
    
    def write_forecast(self, t, data,ens_idx):
        self._ensure_batch(t)
        batch_idx = t - self.current_batch_start
        start = MPI.Wtime()
        # self.datasets[batch_idx][:, self.nd_start:self.nd_end] = data.T
        # self.datasets[batch_idx][self.nd_start:self.nd_end, :] = data
        self.datasets[batch_idx][self.nd_start:self.nd_end, ens_idx] = data    
        # self.datasets[batch_idx][:, ens_idx] = data
        write_time = MPI.Wtime() - start
        if self.rank == 0:
            print(f"t={t}, Write: {write_time:.2f}s")

    def read_analysis(self, t, ens_idx):
        if t not in self.tobserve:
            raise ValueError(f"Time step {t} is not an observation time")
        self._ensure_batch(t)
        batch_idx = t - self.current_batch_start
        start = MPI.Wtime()
        # data = self.datasets[batch_idx][:, self.nd_start:self.nd_end]
        # data = self.datasets[batch_idx][self.nd_start:self.nd_end, :]
        data = self.datasets[batch_idx][self.nd_start:self.nd_end, ens_idx]
        read_time = MPI.Wtime() - start
        if self.rank == 0:
            print(f"t={t}, Analysis read: {read_time:.2f}s")
        # return data.T
        return data

    def write_analysis(self, t, data, ens_idx):
        if t not in self.tobserve:
            raise ValueError(f"Time step {t} is not an observation time")
        self._ensure_batch(t)
        batch_idx = t - self.current_batch_start
        start = MPI.Wtime()
        # self.datasets[batch_idx][:, self.nd_start:self.nd_end] = data.T
        # self.datasets[batch_idx][self.nd_start:self.nd_end, :] = data
        self.datasets[batch_idx][self.nd_start:self.nd_end, ens_idx] = data
        write_time = MPI.Wtime() - start
        if self.rank == 0:
            print(f"t={t}, Analysis write: {write_time:.2f}s")

    def close(self):
        self._close_batch()

# Dummy forecast and analyze functions for testing
def forecast(state):
    return state  # Placeholder: modify state as needed

def analyze(state):
    return state  # Placeholder: modify state as needed

if __name__ == "__main__":
    # Example parameters
    nd = 36231
    nens = 100
    nt = 500
    # tobserve = [0, 100, 200, 300, 400, 500]
    tobserve = [0, 5, 10, 15, 19]
    comm = MPI.COMM_WORLD

    start_time = MPI.Wtime()

    if nens >= comm.Get_size():
        subcomm_size = min(comm.Get_size(), nens)
        color = comm.Get_rank() % subcomm_size
        key = comm.Get_rank() // subcomm_size
        rounds = (comm.Get_size() + subcomm_size - 1) // subcomm_size
    else:
        color = comm.Get_rank() % nens if comm.Get_rank() < nens else MPI.UNDEFINED
        key = comm.Get_rank() // nens
        rounds = 1  # Only one round of processing needed

    # Create subcommunicator
    subcomm = comm.Split(color, key)

    enkf_io = EnKFIO('enkf', nd, nens, nt, tobserve, subcomm, comm, batch_size=50)

    # Simulation loop
    for t in range(nt):

        for round_idx in range(rounds):
            ens_id = color + round_idx * subcomm_size
            if ens_id < nens:
                state = enkf_io.read_forecast(t,ens_id)
                print (f"Rank {subcomm.Get_rank()}, time {t}, ens_id {ens_id}, state shape: {state.shape}")
                state = forecast(state)
                enkf_io.write_forecast(t + 1 if t < nt - 1 else t, state,ens_id)

        # state = enkf_io.read_forecast(t)
        # state = forecast(state)
        # enkf_io.write_forecast(t + 1 if t < nt - 1 else t, state)
        if t in tobserve:
            state = enkf_io.read_analysis(t,ens_id)
            print(f"Rank {subcomm.Get_rank()}, Analysis read at time {t}, state shape: {state.shape}")
            state = analyze(state)
            enkf_io.write_analysis(t, state, ens_id)

    enkf_io.close()

    end_time = MPI.Wtime()
    walltime = comm.reduce(end_time - start_time, op=MPI.MAX, root=0)
    exec_time = comm.reduce(end_time - start_time, op=MPI.SUM, root=0)
    if comm.Get_rank() == 0:
        print(f"Total execution time: {exec_time:.2f} seconds, Wall time: {walltime:.2f} seconds")
        print("Test completed successfully.")