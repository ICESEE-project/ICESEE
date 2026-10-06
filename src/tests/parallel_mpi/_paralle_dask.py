import h5py
import numpy as np
from mpi4py import MPI
import glob
import os
import time

class EnKFIO:
    """Manages parallel I/O for EnKF with subcommunicators, rounds, and optimized file access."""
    
    def __init__(self, nd, nt, nens, base_path="enkf_data", chunk_size=1000, num_creator_ranks=2, serial_creation=False, pre_open_files=False):
        """
        Initialize EnKF I/O with subcommunicators, dynamically setting nd if files exist.
        
        Args:
            nd (int): Number of state dimensions (overridden if files exist).
            nt (int): Number of timesteps.
            nens (int): Number of ensemble members.
            base_path (str): Directory for HDF5 files.
            chunk_size (int): Chunk size for state dimension.
            num_creator_ranks (int): Number of ranks to use for file creation.
            serial_creation (bool): If True, create files serially on rank 0.
            pre_open_files (bool): If True, open all local files in __init__.
        """
        self.nt = nt
        self.nens = nens
        self.base_path = base_path
        self.num_creator_ranks = min(num_creator_ranks, MPI.COMM_WORLD.Get_size(), nens) if not serial_creation else 1
        self.serial_creation = serial_creation
        self.pre_open_files = pre_open_files
        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.size = self.comm.Get_size()
        self.files = {}
        self.datasets = {}
        
        if self.rank == 0:
            # List all files starting with 'ensemble_' and 'temp_analysis' in base_path
            patterns = [f"{self.base_path}/ensemble_*.h5", f"{self.base_path}/temp_analysis*.h5"]
            for pattern in patterns:
                matching_files = glob.glob(pattern)
                for file_path in matching_files:
                    try:
                        os.remove(file_path)
                        # print(f"Deleted: {file_path}")
                    except OSError as e:
                        print(f"Error deleting {file_path}: {e}")

        self.comm.Barrier()
        
        self.nd = nd
        self.chunk_size = min(chunk_size, self.nd)
        
        # Create subcommunicator
        self.subcomm, self.local_ens = self._create_subcomm_and_distribute()
        print(f"Rank {self.rank}: size={self.size}, local_ens={self.local_ens}, subcomm_size={self.subcomm.Get_size() if self.subcomm != MPI.COMM_NULL else 'NULL'}")
        
        # Create directory
        if self.rank == 0:
            os.makedirs(base_path, exist_ok=True)
        self.comm.Barrier()
        
        # Create files
        self._create_files()
        
        # Pre-open files for local ensembles
        if self.pre_open_files and self.subcomm != MPI.COMM_NULL:
            for ens_idx in self.local_ens:
                file_path = f"{self.base_path}/ensemble_{ens_idx}.h5"
                try:
                    self.files[ens_idx] = h5py.File(file_path, "r+", driver=None)
                    self.datasets[ens_idx] = self.files[ens_idx]["state"]
                    print(f"Rank {self.rank} pre-opened {file_path}, shape={self.datasets[ens_idx].shape}")
                except Exception as e:
                    print(f"Rank {self.rank} error pre-opening {file_path}: {e}")
                    self.comm.Abort(1)

    def _create_subcomm_and_distribute(self):
        """Create subcommunicator and distribute ensemble members."""
        if self.nens >= self.size:
            # Divide ranks into subcommunicators
            subcomm_size = min(self.size, self.nens)  # Use at most `nens` groups
            color = self.rank % subcomm_size  # Group ranks into subcommunicators
            key = self.rank // subcomm_size  # Ordering within each subcommunicator
            rounds = (self.nens + subcomm_size - 1) // subcomm_size  # Ceiling division
        else:
            # More processes than ensembles, map processes to ensembles efficiently
            color = self.rank % self.nens if self.rank < self.nens else MPI.UNDEFINED
            key = self.rank // self.nens
            rounds = 1  # Only one round of processing needed
        
        # Create subcommunicator
        subcomm = self.comm.Split(color, key)
        
        local_ens = []
        if subcomm != MPI.COMM_NULL:
            sub_rank = subcomm.Get_rank()
            sub_size = subcomm.Get_size()
            if self.nens >= self.size:
                # Distribute ensembles across rounds
                for round_id in range(rounds):
                    current_ens = color + round_id * subcomm_size
                    if current_ens < self.nens:
                        local_ens.append(current_ens)
            else:
                local_ens = [color] if color < self.nens else []
        
        return subcomm, local_ens

    def _create_files(self):
        """Create HDF5 files, optionally serially on rank 0."""
        skip_creation = False
        if self.rank == 0:
            skip_creation = all(os.path.exists(f"{self.base_path}/ensemble_{i}.h5") for i in range(self.nens))
            if skip_creation:
                print(f"All {self.nens} files exist, skipping creation")
        
        skip_creation = self.comm.bcast(skip_creation, root=0)
        
        if skip_creation:
            self.comm.Barrier()
            return
        
        if self.rank == 0:
            print(f"Creating {self.nens} files {'serially on rank 0' if self.serial_creation else f'using {self.num_creator_ranks} ranks'}, nd={self.nd}, nt={self.nt}")
        start_time = time.time()
        
        if self.serial_creation:
            if self.rank == 0:
                for ens_idx in range(self.nens):
                    file_start = time.time()
                    file_path = f"{self.base_path}/ensemble_{ens_idx}.h5"
                    try:
                        open_start = time.time()
                        with h5py.File(file_path, "w", libver="latest") as f:
                            open_time = time.time() - open_start
                            create_start = time.time()
                            dset = f.create_dataset(
                                "state",
                                (self.nd, self.nt),
                                dtype=np.float64,
                                chunks=(self.chunk_size, 1),
                                fillvalue=None
                            )
                            create_time = time.time() - create_start
                        file_time = time.time() - file_start
                        print(f"Rank {self.rank} created {file_path} in {file_time:.3f} seconds (open: {open_time:.3f}s, create: {create_time:.3f}s), shape=({self.nd}, {self.nt})")
                    except Exception as e:
                        print(f"Rank {self.rank} error creating {file_path}: {e}")
                        self.comm.Abort(1)
        else:
            if self.rank < self.num_creator_ranks:
                files_per_rank = self.nens // self.num_creator_ranks
                remainder = self.nens % self.num_creator_ranks
                start_idx = self.rank * files_per_rank + min(self.rank, remainder)
                count = files_per_rank + 1 if self.rank < remainder else files_per_rank
                end_idx = start_idx + count
                
                for ens_idx in range(start_idx, end_idx):
                    file_start = time.time()
                    file_path = f"{self.base_path}/ensemble_{ens_idx}.h5"
                    try:
                        open_start = time.time()
                        with h5py.File(file_path, "w", libver="latest") as f:
                            open_time = time.time() - open_start
                            create_start = time.time()
                            dset = f.create_dataset(
                                "state",
                                (self.nd, self.nt),
                                dtype=np.float64,
                                chunks=(self.chunk_size, 1),
                                fillvalue=None
                            )
                            create_time = time.time() - create_start
                        file_time = time.time() - file_start
                        print(f"Rank {self.rank} created {file_path} in {file_time:.3f} seconds (open: {open_time:.3f}s, create: {create_time:.3f}s), shape=({self.nd}, {self.nt})")
                    except Exception as e:
                        print(f"Rank {self.rank} error creating {file_path}: {e}")
                        self.comm.Abort(1)
        
        self.comm.Barrier()
        total_time = time.time() - start_time
        if self.rank == 0:
            print(f"Total file creation time: {total_time:.2f} seconds")

    def write_forecast(self, ens_idx, t, state):
        """Write state for an ensemble member at a timestep."""
        if self.subcomm == MPI.COMM_NULL or ens_idx not in self.local_ens:
            return
        file_path = f"{self.base_path}/ensemble_{ens_idx}.h5"
        if ens_idx not in self.files:
            try:
                self.files[ens_idx] = h5py.File(file_path, "r+", driver=None)
                self.datasets[ens_idx] = self.files[ens_idx]["state"]
                print(f"Rank {self.rank} opened {file_path} for writing, shape={self.datasets[ens_idx].shape}")
            except Exception as e:
                print(f"Rank {self.rank} error opening {file_path} for writing: {e}")
                self.comm.Abort(1)
        dataset_shape = self.datasets[ens_idx].shape
        if state.shape != (dataset_shape[0],):
            raise ValueError(f"Rank {self.rank}: State shape {state.shape} does not match dataset shape ({dataset_shape[0]},) for ens_idx={ens_idx}, t={t}")
        self.datasets[ens_idx][:, t] = state

    def read_forecast(self, ens_idx, t):
        """Read state for an ensemble member at a timestep."""
        if self.subcomm == MPI.COMM_NULL or ens_idx not in self.local_ens:
            return None
        file_path = f"{self.base_path}/ensemble_{ens_idx}.h5"
        if ens_idx not in self.files:
            try:
                self.files[ens_idx] = h5py.File(file_path, "r+", driver=None)
                self.datasets[ens_idx] = self.files[ens_idx]["state"]
                print(f"Rank {self.rank} opened {file_path} for reading, shape={self.datasets[ens_idx].shape}")
            except Exception as e:
                print(f"Rank {self.rank} error opening {file_path} for reading: {e}")
                self.comm.Abort(1)
        return self.datasets[ens_idx][:, t]

    def gather_analysis(self, t):
        """Aggregate states to a temporary shared file for a timestep using subcomm."""
        temp_file_path = f"{self.base_path}/temp_analysis_t{t}.h5"
        
        # Let rank 0 create the temp file
        if self.rank == 0:
            try:
                with h5py.File(temp_file_path, "w", libver="latest") as f:
                    f.create_dataset(
                        "analysis",
                        (self.nd, self.nens),
                        dtype=np.float64,
                        chunks=(self.chunk_size, 1)
                    )
            except Exception as e:
                print(f"Rank {self.rank} error creating {temp_file_path}: {e}")
                self.comm.Abort(1)
        self.comm.Barrier()

        # Read local ensemble states
        local_nens = len(self.local_ens)
        local_data = np.zeros((self.nd, local_nens), dtype=np.float64) if local_nens > 0 else np.zeros((self.nd, 0))
        for i, ens_idx in enumerate(self.local_ens):
            file_path = f"{self.base_path}/ensemble_{ens_idx}.h5"
            try:
                with h5py.File(file_path, "r", libver="latest", swmr=True) as f:
                    # Read chunk-wise to manage memory
                    for start in range(0, self.nd, self.chunk_size):
                        end = min(start + self.chunk_size, self.nd)
                        local_data[start:end, i] = f["state"][start:end, t]
            except Exception as e:
                print(f"Rank {self.rank} error reading {file_path} for ens_idx={ens_idx}, t={t}: {e}")
                self.comm.Abort(1)

        # Determine column indices for writing
        if self.subcomm != MPI.COMM_NULL:
            sub_rank = self.comm.Get_rank()
            sub_size = self.comm.Get_size()
            members_per_rank = self.nens // sub_size
            remainder = self.nens % sub_size
            start_idx = sub_rank * members_per_rank + min(sub_rank, remainder)
            count = members_per_rank + 1 if sub_rank < remainder else members_per_rank
            end_idx = start_idx + count
        else:
            start_idx, end_idx = 0, 0

        # Write to temporary HDF5 file using parallel I/O
        # if local_data.size > 0:
        #     try:
        #         with h5py.File(temp_file_path, "a", driver="mpio", comm=self.comm) as f:
        #             dset = f["analysis"]
        #             # Write chunk-wise to manage memory
        #             for start in range(0, self.nd, self.chunk_size):
        #                 end = min(start + self.chunk_size, self.nd)
        #                 dset[start:end, start_idx:end_idx] = local_data[start:end, :]
        #     except Exception as e:
        #         print(f"Rank {self.rank} error writing to {temp_file_path}: {e}")
        #         self.comm.Abort(1)

        # self.subcomm.Barrier()
        
        # Read the full matrix for return (optional, can skip if only needed on disk)
        try:
            with h5py.File(temp_file_path, "r", libver="latest", swmr=True) as f:
                analysis_data = np.array(f["analysis"])
        except Exception as e:
            print(f"Rank {self.rank} error reading {temp_file_path}: {e}")
            self.comm.Abort(1)
        
        # Clean up temporary file
        self.comm.Barrier()
        if self.rank == 0:
            try:
                os.remove(temp_file_path)
                print(f"Rank {self.rank} deleted {temp_file_path}")
            except OSError as e:
                print(f"Rank {self.rank} error deleting {temp_file_path}: {e}")

        print(f"Rank {self.rank} completed gather_analysis for t={t}, shape: {analysis_data.shape}")
        return analysis_data

    def close(self):
        """Close all open HDF5 files."""
        for ens_idx in list(self.files.keys()):
            try:
                self.files[ens_idx].close()
            except Exception as e:
                print(f"Rank {self.rank} error closing file for ens_idx={ens_idx}: {e}")
            del self.files[ens_idx]
            del self.datasets[ens_idx]

# Test case
if __name__ == "__main__":
    import glob
    import os
    import time
    import numpy as np
    from mpi4py import MPI

    # Parameters
    nd = 425000  # Number of state dimensions
    nt = 250
    nens = 100
    chunk_size = 1000
    num_creator_ranks = 2
    serial_creation = False
    pre_open_files = False  # Pre-open files to reduce concurrent opens
    
    start_time = MPI.Wtime()
    
    print(f"\nTesting with nens={nens}")
    try:
        enkf_io = EnKFIO(nd, nt, nens, chunk_size=chunk_size, num_creator_ranks=num_creator_ranks, 
                        serial_creation=serial_creation, pre_open_files=pre_open_files)
    except Exception as e:
        print(f"Rank {MPI.COMM_WORLD.Get_rank()} error initializing EnKFIO: {e}")
        MPI.COMM_WORLD.Abort(1)

    # Initialize ensemble members
    for ens_idx in enkf_io.local_ens:
        state = np.random.rand(enkf_io.nd)
        try:
            enkf_io.write_forecast(ens_idx, 0, state)
        except Exception as e:
            print(f"Rank {enkf_io.rank} error writing forecast for ens_idx={ens_idx}, t={0}: {e}")
            enkf_io.comm.Abort(1)
    
    for t in range(nt):
        for ens_idx in enkf_io.local_ens:
            # Read the state from the file
            state = enkf_io.read_forecast(ens_idx, t)
            # Write state back to the file
            enkf_io.write_forecast(ens_idx, t, state)

        if t % 2 == 0:
            try:
                all_states = enkf_io.gather_analysis(t)
            except Exception as e:
                print(f"Rank {enkf_io.rank} error in gather_analysis: {e}")
                enkf_io.comm.Abort(1)
    
    enkf_io.close()

    end_time = MPI.Wtime()
    comm = MPI.COMM_WORLD
    walltime = comm.reduce(end_time - start_time, op=MPI.MAX, root=0)
    exec_time = comm.reduce(end_time - start_time, op=MPI.SUM, root=0)
    if enkf_io.rank == 0:
        print(f"Total execution time: {exec_time:.2f} seconds, Wall time: {walltime:.2f} seconds")
        print("Test completed successfully.")