import h5py
import os
import glob
import numpy as np
from mpi4py import MPI

class EnKFIO:
    def __init__(self, file_prefix, nd, nens, nt, tobserve, mpi_comm,base_path="enkf_data"):
        """
        Initialize EnKF I/O manager for (nd, nens) data in nt files.
        
        Args:
            file_prefix (str): Prefix for HDF5 files (e.g., 'enkf' -> 'enkf_0000.h5').
            nd (int): State dimension.
            nens (int): Number of ensemble members.
            nt (int): Total number of time steps.
            tobserve (list): List of observation time steps (0-based indices).
            mpi_comm: MPI communicator (e.g., MPI.COMM_WORLD).
        """
        self.nd = nd
        self.nens = nens
        self.nt = nt
        self.tobserve = tobserve
        self.nt_obs = len(tobserve)
        self.base_path = base_path
        self.comm = mpi_comm
        self.rank = mpi_comm.Get_rank()
        self.size = mpi_comm.Get_size()
        self.file_prefix = file_prefix

        # Partition nd across processes
        self.nd_local = nd // self.size
        self.nd_start = self.rank * self.nd_local
        self.nd_end = self.nd_start + self.nd_local

        # Create directory
        if self.rank == 0:
            os.makedirs(base_path, exist_ok=True)
        self.comm.Barrier()

        if self.rank == 0:
            # remove all file with enkf_*
             # List all files starting with 'ensemble_' and 'temp_analysis' in base_path
            patterns = [f"{self.base_path}/enkf_*.h5"]
            for pattern in patterns:
                matching_files = glob.glob(pattern)
                for file_path in matching_files:
                    try:
                        os.remove(file_path)
                        # print(f"Deleted: {file_path}")
                    except OSError as e:
                        print(f"Error deleting {file_path}: {e}")

        self.comm.Barrier()

        # Create and keep open nt HDF5 files
        self.files = []
        self.datasets = []
        self.file_create_track = 0
        
    def create_files(self, t):
        # check if t>500 and nd > 100000 then create files in groups for every 50 time steps
        if self.nt > 500 and self.nd > 100000:
            start_group_size = 50
            nt = start_group_size
            if t> nt:
                # Create files in groups of start_group_size
                for t in range(t+nt):
                    fname = f"{self.base_path}/{self.file_prefix}_{t:04d}.h5"
                    f = h5py.File(fname, 'w', driver='mpio', comm=self.comm)
                    dset = f.create_dataset(
                        'state', (self.nd, self.nens), chunks=(self.nd_local, self.nens), dtype='f8'
                    )
                    self.files.append(f)
                    self.datasets.append(dset)
            else:
                # Create files in groups of start_group_size
                for t in range(nt):
                    fname = f"{self.base_path}/{self.file_prefix}_{t:04d}.h5"
                    f = h5py.File(fname, 'w', driver='mpio', comm=self.comm)
                    dset = f.create_dataset(
                        'state', (self.nd, self.nens), chunks=(self.nd_local, self.nens), dtype='f8'
                    )
                    self.files.append(f)
                    self.datasets.append(dset)
        else:
            # Create files for each time step
            self.files = []
            self.datasets = []
            for t in range(self.nt):
                fname = f"{self.base_path}/{self.file_prefix}_{t:04d}.h5"
                f = h5py.File(fname, 'w', driver='mpio', comm=self.comm)
                dset = f.create_dataset(
                    'state', (self.nd, self.nens), chunks=(self.nd_local, self.nens), dtype='f8'
                )
                self.files.append(f)
                self.datasets.append(dset)


    def read_forecast(self, t):
        """
        Read (nd_local, nens) for forecast at time t.
        
        Args:
            t (int): Time step index (0-based).
        
        Returns:
            numpy.ndarray: Local forecast state [nd_local × nens].
        """
        return self.datasets[t][self.nd_start:self.nd_end, :]  # [nd_local × nens]

    def write_forecast(self, t, data):
        """
        Write (nd_local, nens) for forecast at time t.
        
        Args:
            t (int): Time step index (0-based).
            data (numpy.ndarray): Local forecast state [nd_local × nens].
        """
        self.datasets[t][self.nd_start:self.nd_end, :] = data

    def read_analysis(self, t):
        """
        Read (nd_local, nens) for analysis at time t (must be in tobserve).
        
        Args:
            t (int): Time step index (0-based).
        
        Returns:
            numpy.ndarray: Local analysis state [nd_local × nens].
        """
        if t not in self.tobserve:
            raise ValueError(f"Time step {t} is not an observation time")
        return self.datasets[t][self.nd_start:self.nd_end, :]

    def write_analysis(self, t, data):
        """
        Write (nd_local, nens) for analysis at time t (must be in tobserve).
        
        Args:
            t (int): Time step index (0-based).
            data (numpy.ndarray): Local analysis state [nd_local × nens].
        """
        if t not in self.tobserve:
            raise ValueError(f"Time step {t} is not an observation time")
        self.datasets[t][self.nd_start:self.nd_end, :] = data

    def close(self):
        """
        Close all HDF5 files.
        """
        for f in self.files:
            f.close()

if __name__ == "__main__":
    import numpy as np
    from mpi4py import MPI
    # from enkf_io import EnKFIO  # Import the class above


    start_time = MPI.Wtime()    

    # Example parameters
    nd = 100000  # State dimension
    nens = 100   # Ensemble members
    nt = 1000     # Total time steps
    tobserve = [0, 100, 200, 300, 400, 500, 600, 700, 800, 900]  # Observation times
    comm = MPI.COMM_WORLD

    # Initialize I/O manager (creates nt files)
    enkf_io = EnKFIO('enkf', nd, nens, nt, tobserve, comm)

    # Simulation loop
    for t in range(nt):
        # Forecast step
        state = enkf_io.read_forecast(t)  # [nd_local × nens]
        # state = forecast(state)  # Your forecast function
        enkf_io.write_forecast(t + 1 if t < nt - 1 else t, state)

        # Analysis step (at observation times)
        if t in tobserve:
            state = enkf_io.read_analysis(t)  # Read forecast for analysis
            # state = analyze(state)  # Your analysis function
            enkf_io.write_analysis(t, state)

    # Clean up (close all files)
    enkf_io.close()

    end_time = MPI.Wtime()
    comm = MPI.COMM_WORLD
    walltime = comm.reduce(end_time - start_time, op=MPI.MAX, root=0)
    exec_time = comm.reduce(end_time - start_time, op=MPI.SUM, root=0)
    if enkf_io.rank == 0:
        print(f"Total execution time: {exec_time:.2f} seconds, Wall time: {walltime:.2f} seconds")
        print("Test completed successfully.")