
import sys, os
import warnings
import numpy as np
warnings.filterwarnings('ignore') 

os.environ["OMP_NUM_THREADS"] = "1"

# firedrake imports
import firedrake
from firedrake import *

import dask.array as da
import dask
from dask.distributed import Client


from mpi4py import MPI

Lx, Ly = 50e2, 12e2
nx, ny = 12,8

b_in, b_out = 20, -40
s_in, s_out = 850, 50


# mpi communicator
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

# split communicator
# === each rank gets a color and a key (best for cases when Nens > size_world)
# color = rank % size
# key = rank
# comm = comm.Split(color, key)
# ========================================

import h5py
import numpy as np
from mpi4py import MPI
import os
import time
import glob

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
        # if self.rank == 0:
        #     print(f"Rank {self.rank}: size={self.size}, local_ens={self.local_ens}, subcomm_size={self.subcomm.Get_size() if self.subcomm != MPI.COMM_NULL else 'NULL'}")
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

    def _create_subcomm_and_distribute_(self):
        """Create subcommunicator and distribute ensemble members."""
    
        # if nens > self.size:
        #     # Divide ranks into `size` subcommunicators
        #     subcomm_size = min(self.size, self.nens)  # Use at most `Nens` groups
        #     color = self.rank % subcomm_size  # Group ranks into `subcomm_size` subcommunicators
        #     key = self.rank // subcomm_size  # Ordering within each subcommunicator
        #     subcomm = comm.Split(color, key)

        #     sub_rank = subcomm.Get_rank()  # Rank within subcommunicator
        #     sub_size = subcomm.Get_size()  # Size of subcommunicator

        #     # Determine the number of processing rounds
        #     rounds = (self.nens + subcomm_size - 1) // subcomm_size  # Ceiling divisio
        #     local_ens = []
        #     for round_id in range(rounds):
        #         ensemble_id = color + round_id * subcomm_size
        #         local_ens.append(ensemble_id) 

        # else:
        #     # Standard case where each rank maps 1-to-1 with an ensemble
        #     color = self.rank % nens  
        #     key = self.rank // nens   
        #     subcomm = comm.Split(color, key)

        #     sub_rank = subcomm.Get_rank()
        #     sub_size = subcomm.Get_size()
        #     local_ens = [color]

        if self.nens <= self.size:
            subcomm_size = 1
            color = self.rank if self.rank < self.nens else MPI.UNDEFINED
        else:
            subcomm_size = self.size
            color = 0
        subcomm = self.comm.Split(color, self.rank)
        
        local_ens = []
        if subcomm != MPI.COMM_NULL:
            sub_rank = subcomm.Get_rank()
            sub_size = subcomm.Get_size()
            
            if self.nens <= self.size:
                if self.rank < self.nens:
                    local_ens = [self.rank]
            else:
                ens_per_rank = self.nens // sub_size
                remainder = self.nens % sub_size
                start = sub_rank * ens_per_rank + min(sub_rank, remainder)
                count = ens_per_rank + 1 if sub_rank < remainder else ens_per_rank
                local_ens = list(range(start, start + count))
        
        return subcomm, local_ens
    
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
            subcomm_size = None
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
                ens_id = color + sub_rank * subcomm_size
                for round_id in range(rounds):
                    current_ens = color + round_id * subcomm_size
                    if current_ens < self.nens:
                        local_ens.append(current_ens)
            else:
                # Each rank gets at most one ensemble
                # if self.rank < self.nens:
                #     local_ens = [self.rank]
                # Each subcommunicator corresponds to one ensemble
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
        # if self.subcomm == MPI.COMM_NULL or ens_idx not in self.local_ens:
        #     return
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
        
        local_nens = len(self.local_ens) if self.subcomm != MPI.COMM_NULL else 0
        local_counts = np.zeros(self.size, dtype=np.int32)
        if self.subcomm != MPI.COMM_NULL:
            local_counts[self.rank] = local_nens
        self.comm.Allreduce(MPI.IN_PLACE, local_counts, op=MPI.SUM)
        displacements = np.cumsum(local_counts) - local_counts
        start_col = displacements[self.rank]
        
        total_nens = np.sum(local_counts)
        if total_nens != self.nens:
            if self.rank == 0:
                print(f"Error: Total gathered ensembles ({total_nens}) != nens ({self.nens})")
            self.comm.Abort(1)
        
        # let rank 0 create the temp file
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
        
        # All ranks write their local ensemble states to the temp file
        if self.subcomm != MPI.COMM_NULL and local_nens > 0:
            try:
                with h5py.File(temp_file_path, "r+", driver="mpio", comm=self.subcomm) as f:
                    dset = f["analysis"]
                    local_states = np.array([self.read_forecast(ens_idx, t) for ens_idx in self.local_ens]).T
                    local_states = np.ascontiguousarray(local_states)
                    with dset.collective:
                        dset[:, start_col:start_col + local_nens] = local_states
            except Exception as e:
                print(f"Rank {self.rank} error writing to {temp_file_path}: {e}")
                self.comm.Abort(1)
        
        self.comm.Barrier()
        
        # Here we use dask to read data from the temp file collectively while performing computations
        result = None
        if self.rank == 0:
            try:
                with h5py.File(temp_file_path, "r") as f:
                    result = f["analysis"][:, :]
                    # result = f["analysis"][:,1]
                print(f"Rank {self.rank} read aggregated data from {temp_file_path}, shape: {result.shape}")
            except Exception as e:
                print(f"Rank {self.rank} error reading {temp_file_path}: {e}")
                self.comm.Abort(1)
        
        # if self.rank == 0:
        #     try:
        #         os.remove(temp_file_path)
        #     except Exception as e:
        #         print(f"Rank {self.rank} error deleting {temp_file_path}: {e}")
        
        return result
    
    def gather_analysis_(self, t):
        """Aggregate states to a temporary shared file for a timestep using subcomm."""
        temp_file_path = f"{self.base_path}/temp_analysis_t{t}.h5"
        
        # let rank 0 create the temp file
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

        # All ranks write their local ensemble states to the temp file
        # create a list of file paths for each ensemble member
        file_paths = [f"{self.base_path}/ensemble_{ens_idx}.h5" for ens_idx in range(self.nens)]
        
        def form_matrix(file_path, t):
            """Read state vector from file and return as a dask array."""
          
            # create dask arrays for collective reading
            dask_arrays = []
            for file_path in file_paths:
                # with h5py.File(file_path, "r", driver="mpio", comm=self.subcomm) as f:
                with h5py.File(file_path, "r") as f:
                    dset = f["state"]
                    # create dask array from the slice[:,timestep], chuncked along nd
                    dask_array = da.from_array(dset[:, t], chunks=(self.chunk_size,))
                    dask_arrays.append(dask_array)

            # Stack the dask arrays along the second axis (ensemble members)
            return da.stack(dask_arrays, axis=1).compute()
        
        # Use dask to read data from the files collectively while performing computations
        client = Client(n_workers=self.size, threads_per_worker=2, processes=True, memory_limit='4GB')

        analysis_data = form_matrix(file_paths, t)

        print(f"Rank {self.rank} read aggregated data from files, shape: {analysis_data.shape}")

        client.close()  # Close the Dask client
        
        
        



    def close(self):
        """Close all open HDF5 files."""
        for ens_idx in list(self.files.keys()):
            try:
                self.files[ens_idx].close()
            except Exception as e:
                print(f"Rank {self.rank} error closing file for ens_idx={ens_idx}: {e}")
            del self.files[ens_idx]
            del self.datasets[ens_idx]

    def collect_ensemble_data_parallel(self, t, max_retries=5, retry_delay=1):
        """
        Collectively reads from nens ensemble_*.h5 files and writes to a single
        temp file of shape (nd, nens) for EnKF analysis using MPI.

        Parameters:
        - t (int): Time index to extract from each file.
        - max_retries (int): Maximum number of retries for file opening.
        - retry_delay (float): Delay between retries in seconds.
        """
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()

        nens = self.nens  # Total number of ensemble members

        local_ids = list(range(rank, nens, size))  # Local ensemble member indices for this rank
        nd = self.nd  # Number of state dimensions

        local_data = {}

        for i in local_ids:
            file_name = f"{self.base_path}/ensemble_{i}.h5"
            with h5py.File(file_name, 'r', driver='mpio', comm=comm) as f:
                # Read the state vector at time t
                data = f['state'][:, t]
                local_data[i] = data

        comm.Barrier()  # Ensure all ranks have read their data

        # write collectively to a temporary file
        temp_file_path = f"{self.base_path}/temp_analysis_t{t}.h5"
        with h5py.File(temp_file_path, 'w', driver='mpio', comm=comm) as f:
            dset = f.create_dataset('analysis', (nd, nens), dtype=np.float64)

            # Each rank writes its local data to the corresponding columns
            for i in local_ids:
                if i in local_data:
                    dset[:, i] = local_data[i]



        # # Disable HDF5 file locking for input files
        # os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

        # # Determine dimensions from the first file (only rank 0 reads for consistency)
        # nd, nt = None, None
        # if rank == 0:
        #     file_path = f"{self.base_path}/ensemble_1.h5"
        #     for attempt in range(max_retries):
        #         try:
        #             with h5py.File(file_path, 'r') as f:
        #                 data = f['data']
        #                 nd, nt = data.shape
        #                 if t >= nt or t < 0:
        #                     raise ValueError(f"Time index t={t} out of bounds for nt={nt}")
        #             break
        #         except OSError as e:
        #             if attempt < max_retries - 1:
        #                 time.sleep(retry_delay)
        #                 continue
        #             raise OSError(f"Failed to open {file_path} after {max_retries} attempts: {e}")
        
        # # Broadcast dimensions to all processes
        # comm.Barrier()  # Ensure rank 0 has finished reading
        # nd, nt = comm.bcast((nd, nt), root=0)

        # # Divide ensemble members among processes
        # local_nens = nens // size
        # remainder = nens % size
        # start_idx = rank * local_nens + min(rank, remainder)
        # local_nens += 1 if rank < remainder else 0
        # end_idx = start_idx + local_nens

        # # Create output file with collective access
        # temp_file_path = f"{self.base_path}/temp_analysis_t{t}.h5"
        # with h5py.File(temp_file_path, 'w', driver='mpio', comm=comm) as f:
        #     dset = f.create_dataset('analysis', (nd, nens), dtype=np.float64)

        #     # Each process reads its assigned ensemble files and writes to the dataset
        #     for i in range(start_idx, end_idx):
        #         file_idx = i + 1  # Files are 1-indexed
        #         file_name = f"{self.base_path}/ensemble_{file_idx}.h5"
        #         for attempt in range(max_retries):
        #             try:
        #                 with h5py.File(file_name, 'r') as f_in:
        #                     data = f_in['data'][:, t]  # Read state vector at time t (shape: nd)
        #                     # Collectively write to the corresponding column
        #                     with dset.collective:
        #                         dset[:, i] = data
        #                 break
        #             except OSError as e:
        #                 if attempt < max_retries - 1:
        #                     time.sleep(retry_delay)
        #                     continue
        #                 raise OSError(f"Failed to open {file_name} after {max_retries} attempts: {e}")
        
        # comm.Barrier()  # Ensure all writes are complete
# Test case
if __name__ == "__main__":
    # Parameters
    nd = 425000  # Number of state dimensions
    nt = 250
    nens = 100
    chunk_size = 1000
    num_creator_ranks = 2
    serial_creation = False
    pre_open_files = False  # Pre-open files to reduce concurrent opens
    
    start_time = MPI.Wtime()
    
    # for nens in nens_values:
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
            # read the state from the file
            state = enkf_io.read_forecast(ens_idx, t)

            # write state back to the file
            enkf_io.write_forecast(ens_idx, t, state)


        # if t % 2 == 0:
        #     try:
        #         all_states = enkf_io.gather_analysis(t)
        #         # all_states = enkf_io.collect_ensemble_data_parallel(t)
        #         # if enkf_io.rank == 0:
        #         #     print(f"Analysis step at t={t}, shape: {all_states.shape if all_states is not None else 'None'}")
        #     except Exception as e:
        #         # print(f"Rank {enkf_io.rank} error in gather_analysis: {e}")
        #         print(f"Rank {enkf_io.rank} error in collect_ensemble_data_parallel: {e}")
        #         enkf_io.comm.Abort(1)
    
    enkf_io.close()
    
    # if enkf_io.rank == 0:
    #     for ens_idx in range(nens):
    #         try:
    #             with h5py.File(f"enkf_data/ensemble_{ens_idx}.h5", "r") as f:
    #                 data = f["state"][:, nt-1]
    #                 print(f"Ensemble {ens_idx} at t={nt-1}, mean: {np.mean(data):.4f}, shape: {data.shape}")
    #         except Exception as e:
    #             print(f"Error reading {ens_idx}: {e}")
        
        # if enkf_io.rank == 0:
        #     for ens_idx in range(nens):
        #         try:
        #             os.remove(f"enkf_data/ensemble_{ens_idx}.h5")
        #         except Exception:
        #             pass
        #     try:
        #         os.rmdir("enkf_data")
        #     except Exception:
        #         pass

    end_time = MPI.Wtime()
    comm = MPI.COMM_WORLD
    walltime = comm.reduce(end_time - start_time, op=MPI.MAX, root=0)
    exec_time = comm.reduce(end_time - start_time, op=MPI.SUM, root=0)
    if enkf_io.rank == 0:
        print(f"Total execution time: {exec_time:.2f} seconds, Wall time: {walltime:.2f} seconds")
        print("Test completed successfully.")

# # --- method 1: run model squentially with all ranks for Nens times
# def func(unsplitted_comm,h):
#     # mesh
#     mesh = firedrake.RectangleMesh(nx, ny, Lx, Ly, quadrilateral=True,comm=unsplitted_comm)


#     Q = firedrake.FunctionSpace(mesh, "CG", 2)
#     V = firedrake.VectorFunctionSpace(mesh, "CG", 2)
#     x,y = firedrake.SpatialCoordinate(mesh)

#     b = firedrake.interpolate(b_in - (b_in - b_out) * x / Lx, Q)
#     s0 = firedrake.interpolate(s_in - (s_in - s_out) * x / Lx, Q)
#     h0 = firedrake.interpolate(s0 - b, Q)

#     # get the size of the function space
#     h = Function(Q)
#     # print(f"Size of the function space: {h.dat.data.size} on rank {rank}")
#     return h0.dat.data_ro

# def forcast_squential(comm,Nens, ensemble):
#     # ensemble = np.empty((425, Nens))
#     for ens in range(Nens):
#         comm.barrier() # make sure all ranks are synchronized before running the model
#         h = func(comm)
#         # comm.barrier() # make sure all ranks have run the model before gathering the results
#         h_all = comm.gather(h, root=0)
#         if rank == 0:
#             h_all = np.hstack(h_all)
#             # determine the shape of ensemble (this will be done one at t==0 not here, since the model willbe updating the ensemble at each time step)
#             # if ens == 0:
#             #     ensemble = np.empty((len(h_all), Nens))
#             ensemble[:,ens] = h_all
#             print(f"shape of the gathered array: {ensemble.shape}")
#     return ensemble if rank == 0 else None

# def forecast_parallel_doc(comm, Nens, h):
#     """
#     Runs multiple ensemble members in parallel using `Nens` subcommunicators.
    
#     If `Nens > size`, ensembles are handled in batches using the same subcommunicators.

#     Parameters:
#         comm: MPI communicator (MPI.COMM_WORLD typically)
#         Nens: Number of ensemble members
#         h: Initial state for processing

#     Returns:
#         h_all_global: Final gathered results (only on rank 0)
#     """
#     rank = comm.Get_rank()
#     size = comm.Get_size()

#     import h5py

#     if Nens > size:
#         # Divide ranks into `size` subcommunicators
#         subcomm_size = min(size, Nens)  # Use at most `Nens` groups
#         color = rank % subcomm_size  # Group ranks into `subcomm_size` subcommunicators
#         key = rank // subcomm_size  # Ordering within each subcommunicator
#         subcomm = comm.Split(color, key)

#         sub_rank = subcomm.Get_rank()  # Rank within subcommunicator
#         sub_size = subcomm.Get_size()  # Size of subcommunicator

#         # Determine the number of processing rounds
#         rounds = (Nens + subcomm_size - 1) // subcomm_size  # Ceiling divisio
#     else:
#         # Standard case where each rank maps 1-to-1 with an ensemble
#         color = rank % Nens  
#         key = rank // Nens   
#         subcomm = comm.Split(color, key)

#         sub_rank = subcomm.Get_rank()
#         sub_size = subcomm.Get_size()

#     if sub_rank == 0:
#         file_name = f"ensemble_data_{color}.h5"
#         # all sub ranks open file uisng MPI I/O
#         with h5py.File(file_name, "w", driver='mpio', comm=subcomm) as f:
#             # Create dataset for ensemble members
#             f.create_dataset('ensemble', (425, Nens), dtype='f8')

#     if Nens > size:
       
#         h_list = []
#         for round_id in range(rounds):
#             ensemble_id = color + round_id * subcomm_size  # Global ensemble index

#             if ensemble_id < Nens:  # Only process valid ensembles
#                 print(f"Rank {rank} processing ensemble {ensemble_id} in round {round_id + 1}/{rounds}")

#                 subcomm.Barrier()  # Synchronize before execution

#                 # Run the function in parallel within each subcommunicator
#                 ensemble_vec = func(subcomm, h)  # Each subcommunicator runs the function independently

#                 subcomm.Barrier()  # Ensure all ranks finish

#                 # Gather results within each subcommunicator
#                 # h_all_sub = subcomm.gather(h, root=0)
#                 local_shape = ensemble_vec.size
#                 local_shapes = subcomm.gather(local_shape, root=0)

#                 if sub_rank == 0:
#                     total_size = sum(local_shapes)
#                     counts = local_shapes
#                     displs = [sum(counts[:i]) for i in range(subcomm.Get_size())]
#                     gathered_ensemble = np.empty(total_size, dtype=ensemble_vec.dtype)
#                 else:
#                     gathered_ensemble = None
#                     counts = None
#                     displs = None

#                 subcomm.Gatherv(ensemble_vec[:], [gathered_ensemble, counts, displs, MPI.DOUBLE], root=0)

#                 if sub_rank == 0:
#                     h_all_sub = np.column_stack(gathered_ensemble)  # Stack gathered data
#                     print(f"Subcomm {color} (Rank {sub_rank}) gathered array shape: {h_all_sub.shape}")

#                 h_list.append(h_all_sub if sub_rank == 0 else None)  # Collect only from sub_rank 0

#         # Gather results from subcommunicators to the global rank 0
#         h_all_global = comm.gather(h_list, root=0)  

#     else:
        
#         subcomm.Barrier()  # Synchronize

#         h = func(subcomm, h)  

#         subcomm.Barrier()

#         h_all_sub = subcomm.gather(h, root=0)

#         if sub_rank == 0:
#             h_all_sub = np.hstack(h_all_sub)  
#             print(f"Subcomm {color} (Rank {sub_rank}) gathered array shape: {h_all_sub.shape}")

#         h_all_global = comm.gather(h_all_sub, root=0)

#     if rank == 0:
#         # h_all_global = [arr for sublist in h_all_global for arr in sublist if arr is not None]
#         if Nens < size:
#             h_all_global = [arr for arr in h_all_global if arr is not None]
#         else:
#             h_all_global = [arr for sublist in h_all_global for arr in sublist if arr is not None]
#         # h = np.column_stack(h_all_global)
#         h = np.column_stack(h_all_global)

#         print(f"Final ensemble shape: {h.shape}")

#         return h_all_global if rank == 0 else None

# # run the main function
# if __name__ == "__main__":
#     import time
#     comm = MPI.COMM_WORLD
#     rank = comm.Get_rank()
#     size = comm.Get_size()

#     Nens = 4
#     ensemble= np.empty((425, Nens))
#     start = time.time()
#     # forcast_squential(comm,Nens, ensemble)
#     # forcast_parallel(comm,Nens)
#     forecast_parallel_doc(comm, Nens, ensemble)
#     # forecast_parallel_dynamic(comm, Nens, ensemble)

#     stop = time.time()
#     # get total time from all ranks
#     total_time = comm.reduce(stop-start, op=MPI.MAX, root=0)
#     if rank == 0:
#         print(f" Total time taken using {size} porcs: {total_time/60} minutes")

