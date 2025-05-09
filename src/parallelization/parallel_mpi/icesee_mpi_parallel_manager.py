# =============================================================================
# @author: Brian Kyanjo
# @date: 2025-02-10
# @description: Initializes model communicators: size, rank, and variables for 
#               MPI parallelization to be shared by all modules and the ICESS 
#               package.
# =============================================================================

from mpi4py import MPI
import numpy as np
import math
import copy
import gc

class ParallelManager:
    """
    This class provides variables for MPI parallelization to be shared 
    between model-related routines.

    Implements Singleton pattern to ensure a single instance is used.
    """

    _instance = None  # Singleton instance

    def __new__(cls):
        """Singleton pattern to ensure only one instance of ParallelManager exists."""
        if cls._instance is None:
            cls._instance = super(ParallelManager, cls).__new__(cls)
            cls._instance._initialized = False  # Track initialization state
        return cls._instance

    def __init__(self):
        """Initializes MPI communicator, size, rank, and variables."""
        if self._initialized:
            return  # Avoid re-initialization

        self._initialized = True  # Mark as initialized
        self._initialize_variables()

    def _initialize_variables(self):
        """Initializes default variables for MPI communication."""
        # Global communicator
        self.COMM_WORLD = None  # MPI communicator for all PEs
        self.rank_world = None  # Rank in MPI_COMM_WORLD
        self.size_world = None  # Size of MPI_COMM_WORLD

        # ICESEE sub-communicators
        self.comm_sub = None  # Sub-communicator for ensemble tasks
        self.rank_sub = None  # Rank in the sub-communicator
        self.size_sub = None  # Size of the sub-communicator
        self.rounds   = None  # Number of processing rounds
        self.color    = None  # Sub-communicator group identifier
        self.key      = None  # Ordering within each sub-communicator

        # Model communicator (forecast step)
        self.COMM_model = None  # MPI communicator for model tasks
        self.rank_model = None  # Rank in the COMM_model communicator
        self.size_model = None  # Size of the COMM_model communicator


        # ICESS communicator for the analysis step
        self.COMM_filter = None  # MPI communicator for filter PEs
        self.rank_filter = None  # Rank in COMM_filter
        self.size_filter = None  # Number of PEs in COMM_filter
        self.n_filterpes = 1  # Number of parallel filter analysis tasks

        # ICESS communicator for coupling filter and model and also for initializations
        self.COMM_couple  = None  # MPI communicator for coupling filter and model
        self.rank_couple  = None  # Rank in COMM_couple
        self.size_couple  = None  # Number of PEs in COMM_couple
        self.ens_id       = None  # Index of ensemble member (1,...,Nens)
        self.model_nprocs = None  # Number of PEs in the model communicator

        # ICESS variables
        self.n_modeltasks = None  # Number of parallel model tasks
        self.modelpe = False  # Whether PE is in COMM_model
        self.filterpe = False  # Whether PE is in COMM_filter
        self.task_id = None  # Index of model task (1,...,n_modeltasks)
        self.MPIerr = 0  # Error flag for MPI
        self.local_size_model = None  # Number of PEs per ensemble

        # Shared memory variables
        self. node_comm = None  # MPI communicator for shared memory
        self.node_rank = None  # Rank in the node communicator
        self.node_size = None  # Size of the node communicator
        self.shared_array = None  # Shared memory array
        self.win = None  # MPI window for shared memory
        self.mem_buf = None  # Memory buffer for shared memory
        self.mem_size = None  # Size of the shared memory buffer (Nens*state_block_size)
        self.local_start = None  # Start index of local ensemble
        self.local_stop = None  # Stop index of local ensemble


    def init_parallel_non_mpi_model(self):
        """
        Initializes MPI in a non-MPI model.

        Determines the number of PEs (`size_world`) and the rank of a PE (`rank_world`).
        The model is executed within the scope of `COMM_model`.
        """
        if not MPI.Is_initialized():
            MPI.Init()

        self.COMM_model = MPI.COMM_WORLD
        self.size_world = self.COMM_model.Get_size()
        self.rank_world = self.COMM_model.Get_rank()

        # Use MPI_COMM_WORLD as the model communicator
        self.size_model = self.size_world
        self.rank_model = self.rank_world

        #  return the initialized parallel manager
        return self
    
    def initialize_seed(self, comm_world, base_seed=None):
        """
        Initialize the random seed across all MPI processes.
        
        Parameters:
        - comm_world (MPI communicator): Global MPI communicator.
        - base_seed (int or None): Base seed value. If None, generate from rank 0.
        
        Returns:
        - seed (int): The synchronized seed for this process.
        """

        rank_world = comm_world.Get_rank()

        if rank_world == 0:
            # Rank 0 generates or uses the base seed
            seed = base_seed if base_seed is not None else np.random.randint(0, 1000000)
        else:
            seed = None
        
        # Broadcast the seed to all ranks
        seed = comm_world.bcast(seed, root=0)
        
        # Set the seed for NumPy's RNG
        np.random.seed(seed + rank_world)  # Offset by rank for unique sequences on every rank
        
        return seed



    def icesee_mpi_init(self, params):
        """
        Initializes MPI in an ICESEE application.
        Ensures MPI is initialized safely and retrieves essential communication parameters.
        """
        import os, shutil
        if not MPI.Is_initialized():
            try:
                MPI.Init()
            except Exception as e:
                raise RuntimeError(f"[ICESEE] MPI failed to initialize: {e}")

        self.COMM_WORLD = MPI.COMM_WORLD
        self.size_world = self.COMM_WORLD.Get_size()
        self.rank_world = self.COMM_WORLD.Get_rank()

        # if self.model_nprocs is None: use model_nprocs = size_world or subcomm_size
        self.model_nprocs = params.get("model_nprocs")

        # remove data file
        if self.rank_world == 0:
            if os.path.exists(params.get("data_path")):
                if os.path.isdir(params.get("data_path")):
                    shutil.rmtree(params.get("data_path"))
                else:
                    os.remove(params.get("data_path"))
        # Synchronize all ranks
        self.COMM_WORLD.Barrier()

        Nens = params.get("Nens", 1)  # Number of ensemble members

        if params.get("sequential_run", False):
            if self.rank_world == 0: print("[ICESEE] Running sequential mode")
            self.ens_id = None
            return self.rank_world, self.size_world, self.COMM_WORLD, self.ens_id
        
        if params.get("default_run", False):
            if self.rank_world == 0: print("[ICESEE] Running default parallel mode")
            if Nens >= self.size_world: 
                # Divide ranks into `size` subcommunicators
                subcomm_size = min(self.size_world, Nens)  # Use at most `Nens` groups
                self.color = self.rank_world % subcomm_size  # Group ranks into `subcomm_size` subcommunicators
                self.key = self.rank_world // subcomm_size  # Ordering within each subcommunicator
                
                #  here ens_id = number
                self.ens_id = self.color # only needed for initializations as we will only have color ranks available either way
            else:
                # More processes than ensembles, map processes to ensembles efficiently
                self.color = self.rank_world % Nens  
                self.key = self.rank_world // Nens
                self.ens_id = self.color
            
            self.comm_sub = self.COMM_WORLD.Split(self.color, self.key)
            self.rank_sub = self.comm_sub.Get_rank()
            self.size_sub = self.comm_sub.Get_size()

            return self.rank_sub, self.size_sub, self.comm_sub, self.ens_id
        
        if params.get("even_distribution", False):
            if self.rank_world == 0: print("[ICESEE] Running even distribution mode")
            # split the global communicator into size subcommunicators
            self.color = self.rank_world % self.size_world
            self.comm_sub = self.COMM_WORLD.Split(self.color)
            self.rank_sub = self.comm_sub.Get_rank()
            self.size_sub = self.comm_sub.Get_size()
            self.ens_id = self.color

            return self.rank_sub, self.size_sub, self.comm_sub, self.ens_id
        
    def icesee_mpi_ens_distribution(self, params):
        """
        Runs multiple ensemble members in parallel using `Nens` subcommunicators.
        
        Each subcommunicator runs in parallel, and its ranks work together.
        The results are gathered within each subcommunicator and returned.

        Parameters:
            Nens: Number of ensemble members (subcommunicators)
            params: Dictionary of parameters for the ICESEE application

        Returns:
            A tuple containing:
                - rounds: Number of processing rounds required
                - color: Subcommunicator group identifier
                - subcomm: The MPI subcommunicator object
                - sub_rank: Rank within the subcommunicator
                - sub_size: Size of the subcommunicator
                - rank_world: Global MPI rank
                - size_world: Total MPI size
                - comm_world: Global MPI communicator
                - start: Start index for assigned ensembles (only if even_distribution is enabled)
                - stop: Stop index for assigned ensembles (only if even_distribution is enabled)
        """
        Nens = params.get("Nens", 1)  # Number of ensemble members

        comm_world = MPI.COMM_WORLD        # Global communicator
        size_world = comm_world.Get_size() # Number of MPI processes
        rank_world = comm_world.Get_rank() # Rank of this MPI process

        if params.get("sequential_run", False):
            return None, None, None, None, None, None, rank_world, size_world, comm_world, None, None

        if params.get("default_run", False):

            if Nens >= size_world: 
                # Divide ranks into `size` subcommunicators
                subcomm_size = min(size_world, Nens)  # Use at most `Nens` groups
                color = rank_world % subcomm_size  # Group ranks into `subcomm_size` subcommunicators
                key = rank_world // subcomm_size  # Ordering within each subcommunicator
                
                # Determine how many rounds of processing are needed
                rounds = (Nens + subcomm_size - 1) // subcomm_size  # Ceiling division
            
            else:
                # More processes than ensembles, map processes to ensembles efficiently
                color = rank_world % Nens  
                key = rank_world // Nens   
                
                rounds = 1  # Only one round of processing needed
                subcomm_size = None
            
            subcomm = comm_world.Split(color, key)
            # get rank and size for each subcommunicator
            sub_rank = subcomm.Get_rank() # Rank within the subcommunicator
            sub_size = subcomm.Get_size() # Size of the subcommunicator

            return rounds, color, sub_rank, sub_size, subcomm, subcomm_size, rank_world, size_world, comm_world, None, None
        
        if params.get("even_distribution", False):
            # --- Properly Distribute Tasks for All Cases ---
            if Nens >= self.size_world:
                if Nens > size_world:
                    # Case 1: More ensembles than processes → Distribute as evenly as possible
                    mem_per_task = Nens // size_world  # Base number of tasks per process
                    remainder = Nens % size_world       # Extra tasks to distribute

                    if rank_world < remainder:
                        # the first remainder gets mem_per_task+1 tasks each
                        start = rank_world * (mem_per_task + 1)
                        stop = start + (mem_per_task + 1)
                        # stop = start + mem_per_task
                    else:
                        #  the remaining (size - remainder) get mem_per_task tasks each
                        start = rank_world * mem_per_task + remainder
                        stop = start + mem_per_task
                        # stop = start + mem_per_task-1

                elif Nens == size_world:
                    #  Assign at most one task per rank
                    if rank_world < Nens:
                        start, stop = rank_world, rank_world + 1
                    else:
                        # Extra ranks do nothing
                        start, stop = 0, 0

                # split the global communicator into size subcommunicators
                subcomm = comm_world.Split(rank_world % size_world)
                sub_rank = subcomm.Get_rank()
                sub_size = subcomm.Get_size()

                # print(f"[Rank {rank_world}] Processing ensembles {start} to {stop}")

                return None, None, sub_rank, sub_size, subcomm, None, rank_world, size_world, comm_world, start, stop
            else:
                # raise an error if the number of ensembles is less than the number of processes
                raise ValueError("Number of ensembles must be greater than the number of processes or use the default_run mode")
        return None


    # --- Parallel load distribution ---
    def ensembles_load_distribution(self, Nens = None, ensemble =np.empty((100,10))  ,comm=None):
        """
        Distributes ensemble members among MPI processes based on rank and size."""

        global_shape,Nens = ensemble.shape
        rank = comm.Get_rank()
        size = comm.Get_size()

        # --- Properly Distribute Tasks for All Cases ---
        if Nens >= self.size_world:
            if Nens > size:
                # Case 1: More ensembles than processes → Distribute as evenly as possible
                mem_per_task = Nens // size  # Base number of tasks per process
                remainder = Nens % size       # Extra tasks to distribute

                if rank < remainder:
                    # the first remainder gets mem_per_task+1 tasks each
                    start = rank * (mem_per_task + 1)
                    stop = start + (mem_per_task + 1)
                    # stop = start + mem_per_task
                else:
                    #  the remaining (size - remainder) get mem_per_task tasks each
                    start = rank * mem_per_task + remainder
                    stop = start + mem_per_task
                    # stop = start + mem_per_task-1

            elif Nens == size:
                #  Assign at most one task per rank
                if rank < Nens:
                    start, stop = rank, rank + 1
                else:
                    # Extra ranks do nothing
                    start, stop = 0, 0

            # split the global communicator into size subcommunicators
            comm = comm.Split(rank % size)

            subcomm = comm

            print(f"[Rank {rank}] Processing ensembles {start} to {stop}")
        else:
            # Case 2: More processes than ensembles → Assign at most one task per rank
            if True:
                if rank < Nens:
                    start, stop = rank, rank + 1
                else:
                    # Extra ranks do nothing
                    start, stop = 0, 0
                
                # workflow design
                # 1. Split the global communicator (comm) into Nens subcommunicators (subcomm), so each ensemble member has a dedicated group of processes.
                # 2. Each subcommunicator collectively updates its assigned ensemble member (one column of ensemble_state).
                # 3. One process per subcommunicator (sub_rank == 0) gathers the updated results and sends them back to the global root (rank 0).
                # 4. The global root collects the updated ensemble_state and proceeds with the analysis step.
                # 5. The global root broadcasts the updated ensemble_state to all processes.
                # 6. All processes proceed with the analysis step.

        # form  local ensembles
        # ensemble_local = np.zeros((global_shape, stop-start))
        ensemble_local = ensemble[:global_shape,start:stop]
        # for memory issues return a deepcopy of the ensemble_local
        return copy.deepcopy(ensemble_local), start, stop, subcomm
    
    # --- state vector load distribution ---
    def state_vector_load_distribution(self, state_vector,comm):
        """
        Distributes state vector among MPI processes based on rank and size."""

        global_shape, Nens = state_vector.shape
        rank = comm.Get_rank()
        size = comm.Get_size()

        # --- Properly Distribute Tasks for All Cases ---
        if global_shape > self.size_world:
           workloads = [global_shape // size for i in range(size)]
           for i in range(global_shape % size):
               workloads[i] += 1
           start = 0
           for i in range(rank):
               start += workloads[i]
           stop = start + workloads[rank]
        else:
            # Case 2: More processes than state variables → Assign at most one task per rank
            if rank < global_shape:
                start, stop = rank, rank + 1
            else:
                # Extra ranks do nothing
                start, stop = 0, 0
        
        state_vector_local = copy.deepcopy(state_vector[start:stop,:])
        
        return state_vector_local, start, stop
    
    # row vector load distribution
    def icesee_mpi_row_distribution(self, comm, num_rows):
        """ Distribute row vector among MPI processes based on rank and size.
            Only valid if num_rows > size_world.
        """

        rank = comm.Get_rank()
        size = comm.Get_size()

        if num_rows > size:
            # --- Properly Distribute Tasks for All Cases ---
            # Determine rows assigned to each rank
            rows_per_rank = num_rows // size
            extra = num_rows % size  # Handle uneven splits

            if rank < extra:
                local_rows = rows_per_rank + 1
                start_row = rank * (rows_per_rank + 1)
            else:
                local_rows = rows_per_rank
                start_row = rank * rows_per_rank + extra

            end_row = start_row + local_rows
        else:
            # Case 2: More processes than rows → Assign at most one row per rank
            if rank < num_rows:
                local_rows = 1
                start_row = rank
            else:
                # Extra ranks do nothing
                local_rows = 0
                start_row = 0

            end_row = start_row + local_rows

        return rows_per_rank, extra, local_rows, start_row, end_row

    
    # --- memory formulation ---
    def memory_usage(self, global_shape, Nens, bytes_per_element=8):
        """
        Computes the memory usage of an ensemble in bytes."""
        return global_shape * Nens * bytes_per_element/1e9  # Convert to GB

    # ---- Collective Communication Operations ----
    # -- method to gather data from all ranks (many to many)
    def all_gather_data(self, comm, data):
        """
        Gathers data from all ranks using collective communication."""

        size = comm.Get_size()  # Number of MPI processes
        data = np.asarray(data)  # Ensure data is a NumPy array

        # Get the shape of the incoming data
        local_shape = data.shape  # Should be (18915, 16) per rank
        global_shape = (size,) + local_shape  # Expected (size, 18915, 16)

        # Allocate the buffer for gathering
        gathered_data = np.zeros(global_shape, dtype=np.float64)

        # Use Allgather to collect data from all ranks
        comm.Allgather([data, MPI.DOUBLE], [gathered_data, MPI.DOUBLE])

        return gathered_data
    
    # -- method to gather data from all ranks (many to one)
    def gather_data(self, comm, data, root=0):
        """
        Gathers data from all ranks using collective communication."""
        data = np.asarray(data)
        size = comm.Get_size()
        if comm.Get_rank() == root:
            gathered_data = np.empty((size,) + data.shape, dtype=np.float64)
        else:
            gathered_data = None
        comm.Gather([data, MPI.DOUBLE], [gathered_data, MPI.DOUBLE], root=root)
        return gathered_data

    def all_reduce_sum(self, comm, data):
        """
        Reduces data from all ranks using collective communication."""

        size = comm.Get_size()
        data = np.asarray(data)  # Ensure it's an array

        # Allocate buffer for reduced data
        reduced_data = np.zeros_like(data)
        comm.Allreduce([data, MPI.DOUBLE], [reduced_data, MPI.DOUBLE], op=MPI.SUM)

        return reduced_data
    
    # -- method to scatter data to all ranks
    def scatter_data(self,comm, data):
        """
        Scatters data from one rank to all other ranks using collective communication."""

        rank = comm.Get_rank()
        size = comm.Get_size()

        # Ensure data is correctly divided
        local_rows = data.shape[0] // size
        recv_data = np.zeros((local_rows, data.shape[1]), dtype=np.float64)

        # Scatter from Rank 0
        comm.Scatter([data, MPI.DOUBLE], [recv_data, MPI.DOUBLE], root=0)

        return recv_data

    
    # -- method to Bcast data to all ranks
    def broadcast_data(self, comm, data, root=0):
        """
        Broadcasts data from one rank to all other ranks using collective communication."""
        data = np.asarray(data)  # Ensure it's an array
        comm.Bcast([data, MPI.DOUBLE], root=root)
        return data
    
    # -- method to exchange data between all ranks
    def alltoall_exchange(self, comm, data):
        """
        Exchanges data between all ranks using collective communication."""

        size = comm.Get_size()
        local_rows = data.shape[0] // size

        # Each rank prepares send buffer with `size` chunks
        sendbuf = np.split(data, size, axis=0)
        sendbuf = np.concatenate(sendbuf, axis=0)  # Flatten for Alltoall

        # Allocate receive buffer
        recvbuf = np.empty_like(sendbuf)

        # Perform Alltoall communication
        comm.Alltoall([sendbuf, MPI.DOUBLE], [recvbuf, MPI.DOUBLE])

        return recvbuf.reshape(size, local_rows, -1)  # Reshape into proper format

    # --- Point-to-Point Communication Operations ---
    def send_receive_data(self, comm, local_data, source=0, dest=1):
        """
        Sends data from one rank to another using point-to-point communication."""
        rank = comm.Get_rank()

        if rank == source:
            comm.Send([local_data, MPI.DOUBLE], dest=dest)
            print(f"[Rank {rank}] Sent data to Rank {dest}")

        elif rank == dest:
            recv_data = np.empty_like(local_data)
            comm.Recv([recv_data, MPI.DOUBLE], source=source)
            print(f"[Rank {rank}] Received data from Rank {source}")
            return recv_data
        
    # ==== Analysis parallel computations ====
    # -- method to compute the ensemble mean
    def compute_covariance(self,ensemble, mean, node_comm):
        """
        Computes the ensemble covariance matrix in parallel using shared memory.
        """
        N, state_dim = ensemble.shape
        local_ensemble_centered = ensemble - mean  # Center the ensemble
        local_cov = local_ensemble_centered.T @ local_ensemble_centered / (N - 1)

        # Sum covariance across ranks
        global_cov = np.zeros_like(local_cov)
        node_comm.Allreduce([local_cov, MPI.DOUBLE], [global_cov, MPI.DOUBLE], op=MPI.SUM)
        return global_cov
    
    # -- method to compute the ensemble mean
    def compute_mean_from_local_matrix(self,ensemble, node_comm):
        """
        Computes the ensemble mean in parallel using shared memory.
        """
        local_mean = np.mean(ensemble, axis=1)  # Each rank computes mean over its part
        global_mean = np.zeros_like(local_mean)
        node_comm.Allreduce([local_mean, MPI.DOUBLE], [global_mean, MPI.DOUBLE], op=MPI.SUM)
        global_mean /= node_comm.size  # Compute the final mean
        return global_mean
    
    # --- method for matrix-vector multiplication
    def matvec_product(self, matrix, vector, node_comm):
        """
        Computes the matrix-vector product in parallel using shared memory.
        """
        local_product = matrix @ vector
        global_product = np.zeros_like(local_product)
        node_comm.Allreduce([local_product, MPI.DOUBLE], [global_product, MPI.DOUBLE], op=MPI.SUM)
        return global_product
    
    # --- method for matrix-matrix multiplication
    def matmat_product(self, matrix1, matrix2, node_comm):
        """
        Computes the matrix-matrix product in parallel using shared memory.
        """
        local_product = matrix1 @ matrix2
        global_product = np.zeros_like(local_product)
        node_comm.Allreduce([local_product, MPI.DOUBLE], [global_product, MPI.DOUBLE], op=MPI.SUM)
        return global_product
    
    # --- method for Kalamn gain computation
    def compute_kalman_gain(self, ensemble, obs_cov, obs_op, state_dim, node_comm):
        """
        Computes the Kalman gain in parallel using shared memory.
        """
        N, _ = ensemble.shape
        local_cov = ensemble.T @ ensemble / (N - 1)  # Compute ensemble covariance
        local_gain = local_cov @ obs_op.T @ np.linalg.inv(obs_op @ local_cov @ obs_op.T + obs_cov)
        global_gain = np.zeros((state_dim, obs_cov.shape[0]))
        node_comm.Allreduce([local_gain, MPI.DOUBLE], [global_gain, MPI.DOUBLE], op=MPI.SUM)
        return global_gain


    def compute_mean_matrix_from_root(self, X, nd, Nens, comm,root=0):
        """ Compute the mean of each row of a matrix X distributed across processes.
        X has to be available on rank 0 and will be distributed across processes."""
        size = comm.Get_size()
        rank = comm.Get_rank()

        # --- Step 1: Distribute rows of X across processes ---
        local_nd = nd // size  # Number of rows per process
        remainder = nd % size  # Extra rows to distribute

        # Define local row counts
        if rank < remainder:
            local_nd += 1
        start = rank * (nd // size) + min(rank, remainder)
        end = start + local_nd

        # Allocate local array
        local_X = np.empty((local_nd, Nens), dtype=np.float64)

        # Scatter the data
        counts = [((nd // size) + (1 if r < remainder else 0)) * Nens for r in range(size)]
        displs = [sum(counts[:r]) for r in range(size)]

        comm.Scatterv([X, (counts, displs), MPI.DOUBLE], local_X, root=root)

        # --- Step 2: Compute local mean along Nens ---
        local_mean = np.mean(local_X, axis=1)  # Shape: (local_nd,)

        # --- Step 3: Gather local means ---
        # Adjust counts and displs for number of rows (not elements)
        counts_rows = [((nd // size) + (1 if r < remainder else 0)) for r in range(size)]  # Rows per rank
        displs_rows = [sum(counts_rows[:r]) for r in range(size)]  # Offsets in rows

        # Allocate array for all local means on rank 0
        if rank == root:
            all_means = np.empty(nd, dtype=np.float64)
        else:
            all_means = None

        # Gather local means
        comm.Gatherv(local_mean, [all_means, (counts_rows, displs_rows), MPI.DOUBLE], root=root)

        # free memory
        del local_X, local_mean; gc.collect()

        return all_means
    
    def gatherV(self, comm, data, root=0):
        """
        Gathers data from all ranks using collective communication."""
        data = np.asarray(data)
        size = comm.Get_size()
        if comm.Get_rank() == root:
            gathered_data = np.empty((size,) + data.shape, dtype=np.float64)
        else:
            gathered_data = None
        comm.Gather([data, MPI.DOUBLE], [gathered_data, MPI.DOUBLE], root=root)
        return gathered_data





def icesee_mpi_parallelization(Nens, global_shape=1024, n_modeltasks=None, screen_output=True):
    """
    Initializes MPI communicators for parallel model tasks and determines `n_modeltasks` dynamically if not provided.

    Parameters:
        - Nens (int): Number of ensemble members.
        - n_modeltasks (int, optional): Number of parallel model tasks. If `None`, it is auto-determined.
        - screen_output (bool): Whether to print MPI configuration.

    Returns:
        - parallel_manager (ParallelManager): Initialized parallel manager instance.
    """

    parallel_manager = ParallelManager()
    COMM_WORLD = MPI.COMM_WORLD

    # Initialize MPI processing element (PE) information
    parallel_manager.size_world = COMM_WORLD.Get_size()
    parallel_manager.rank_world = COMM_WORLD.Get_rank()

    # --- check if size_world is divisible by Nens for Nens > size_world ---    
    # if Nens > parallel_manager.size_world:
    #     if parallel_manager.size_world < 6:
    #         n_modeltasks = 1

        # size_world should be divisible by Nens
        # if Nens % parallel_manager.size_world != 0:
            # effective_nprocs = find_largest_divisor(Nens, parallel_manager.size_world)
            # if parallel_manager.rank_world == 0:
            #     print(f"\n [ICESEE] Adjusting number of MPI processes from {parallel_manager.size_world} to {effective_nprocs} for even distribution of Ensemble ({Nens})\n")

            # # split MPI processes: only the first effective_nprocs will be used
            # color = 0 if parallel_manager.rank_world < effective_nprocs else 1
            # COMM_WORLD = COMM_WORLD.Split(color, key=parallel_manager.rank_world)
            # parallel_manager.size_world = COMM_WORLD.Get_size()
            # parallel_manager.rank_world = COMM_WORLD.Get_rank()
            # # redefine the n_modeltasks
            # n_modeltasks = max(2, n_modeltasks) if n_modeltasks is not None else None
            # # check if comm_world is redefined
            # if parallel_manager.rank_world == 0:
            #     print(f"\n [ICESEE] Reinitialized communicators with { parallel_manager.size_world} MPI processes\n")
            # raise ValueError(f"Number of MPI processes ({parallel_manager.size_world}) must be divisible by the number of ensemble members ({Nens})")

    # elif parallel_manager.size_world > Nens:
        # size_world should be divisible by Nens
        # if parallel_manager.size_world % Nens != 0:
            # raise ValueError(f"Number of MPI processes ({parallel_manager.size_world}) must be divisible by the number of ensemble members ({Nens})")

    # Display initialization message
    # if parallel_manager.rank_world == 0:
    #     print("\n [ICESEE] Initializing communicators...\n")

    # --- Determine `n_modeltasks` dynamically if not provided ---
    if n_modeltasks is None:
        # if Nens > parallel_manager.size_world:
        #     # Case: More ensembles than processes
        #     n_modeltasks = min(Nens // (parallel_manager.size_world // 2), parallel_manager.size_world // 2)
        # elif Nens <= parallel_manager.size_world:
        #     # Case: Fewer ensembles than processes
        #     n_modeltasks = min(Nens, parallel_manager.size_world // 4)
        # else:
        #     # Case: Roughly equal processes and ensembles
        #     n_modeltasks = min(Nens, parallel_manager.size_world // 2)

        # n_modeltasks = parallel_manager.size_world/(np.log2(parallel_manager.size_world+1))
        # n_modeltasks = int(parallel_manager.size_world / max(1, int(np.ceil(Nens / parallel_manager.size_world))))
        n_modeltasks = math.gcd(Nens, parallel_manager.size_world)
        
        # Ensure `n_modeltasks` is at least 1
        n_modeltasks = max(1, n_modeltasks)

    
    # else:
        # Ensure `n_modeltasks` does not exceed available resources
        # parallel_manager.n_modeltasks = min(n_modeltasks, parallel_manager.size_world, Nens)

    # update the parallel_manager with the number of model tasks
    parallel_manager.n_modeltasks = n_modeltasks

    # --- Check number of parallel ensemble tasks ---
    if parallel_manager.n_modeltasks > parallel_manager.size_world:
        parallel_manager.n_modeltasks = parallel_manager.size_world
    
    # Adjust `n_modeltasks` to ensemble size
    if Nens > 0 and parallel_manager.n_modeltasks > Nens:
        parallel_manager.n_modeltasks = Nens

    # --- Print Optimization Choice ---
    if parallel_manager.rank_world == 0:
        print(f"[ICESEE] Optimized Model Tasks: {parallel_manager.n_modeltasks} "
              f"(for {parallel_manager.size_world} MPI ranks, {Nens} ensembles)")

    # Generate communicator for ensemble tasks
    COMM_ensemble = COMM_WORLD
    size_ens = parallel_manager.size_world
    rank_ens = parallel_manager.rank_world

    # Allocate and distribute PEs per model task
    parallel_manager.local_size_model = np.full(
        parallel_manager.n_modeltasks,
        parallel_manager.size_world // parallel_manager.n_modeltasks,
        dtype=int
    )
    remainder = parallel_manager.size_world % parallel_manager.n_modeltasks
    parallel_manager.local_size_model[:remainder] += 1

    # Assign each PE to a model task
    pe_index = 0
    for i in range(parallel_manager.n_modeltasks):
        for j in range(parallel_manager.local_size_model[i]):
            if rank_ens == pe_index:
                parallel_manager.task_id = i + 1  # Convert to 1-based index
                break
            pe_index += 1
        if parallel_manager.task_id is not None:
            break  # Exit outer loop

    # Create COMM_MODEL communicator
    parallel_manager.COMM_model = COMM_ensemble.Split(color=parallel_manager.task_id, key=rank_ens)
    parallel_manager.size_model = parallel_manager.COMM_model.Get_size()
    parallel_manager.rank_model = parallel_manager.COMM_model.Get_rank()

    # Assign filter PEs
    parallel_manager.filterpe = (parallel_manager.task_id == 1)

    # Create COMM_FILTER communicator
    my_color = parallel_manager.task_id if parallel_manager.filterpe else MPI.UNDEFINED
    parallel_manager.COMM_filter = COMM_WORLD.Split(color=my_color, key=parallel_manager.rank_world)
    
    if parallel_manager.filterpe:
        parallel_manager.size_filter = parallel_manager.COMM_filter.Get_size()
        parallel_manager.rank_filter = parallel_manager.COMM_filter.Get_rank()

    # Create COMM_COUPLE communicator
    color_couple = parallel_manager.rank_model + 1
    parallel_manager.COMM_couple = COMM_ensemble.Split(color=color_couple, key=parallel_manager.rank_world)
    parallel_manager.rank_couple = parallel_manager.COMM_couple.Get_rank()
    parallel_manager.size_couple = parallel_manager.COMM_couple.Get_size()

    # --- split COMM_WORLD into sub-communicators able to use shared memory. 
    # - [node_comm = MPI.COMM_WORLD.Split_type(MPI.COMM_TYPE_SHARED)] 
    # - and then use MPI.Win.Allocate_shared() within each node-local sub-communicator to easily access shared memory. (ranks can efficiently communicate with each othercwithout excessive communication overhead)
    
    # define the size of the shared memory (number of elements)
    # get the ensemble size
    parallel_manager.mem_size = Nens * global_shape

    # get the size(in bytes) of a double precision floating point number
    disp_unit = MPI.DOUBLE.Get_size()

    # split COMM_WORLD into node-local sub-communicators
    parallel_manager.node_comm = COMM_WORLD.Split_type(MPI.COMM_TYPE_SHARED)

    # allocate shared memory
    if False:
        parallel_manager.win = MPI.Win.Allocate_shared(
            parallel_manager.mem_size * disp_unit if parallel_manager.node_comm.rank == 0 else 0,
            disp_unit, 
            comm = parallel_manager.node_comm)
        
        # querying shared memory buffer
        parallel_manager.mem_buf, itemsize = parallel_manager.win.Shared_query(0)
        assert itemsize == MPI.DOUBLE.Get_size()

        # convert the memory buffer to a numpy array
        parallel_manager.mem_buf = np.array(parallel_manager.mem_buf, dtype='B', copy=False)

        # create a shared array
        # parallel_manager.shared_array = np.ndarray(buffer=parallel_manager.mem_buf, dtype='d', shape=( parallel_manager.mem_size,))
        parallel_manager.shared_array = np.ndarray(buffer=parallel_manager.mem_buf, dtype='d', shape=(global_shape, Nens))

        parallel_manager.node_rank = parallel_manager.node_comm.Get_rank()
        parallel_manager.node_size = parallel_manager.node_comm.Get_size()

        # perform ensemble analysis step using shared memory
        parallel_manager.local_start = parallel_manager.node_comm.rank * (Nens // parallel_manager.node_comm.size)
        parallel_manager.local_stop = (parallel_manager.node_comm.rank + 1) * (Nens // parallel_manager.node_comm.size)

    # Display MPI ICESEE Configuration
    if screen_output:
        display_pe_configuration(parallel_manager)

    return parallel_manager


def display_pe_configuration(parallel_manager):
    """
    Displays MPI parallel configuration in a structured format.
    Uses tabulate for better readability.
    """
    from tabulate import tabulate

    COMM_WORLD = MPI.COMM_WORLD
    rank = parallel_manager.rank_world

    rank_info = [
        rank,
        parallel_manager.rank_filter if parallel_manager.filterpe else "-",
        parallel_manager.task_id,
        parallel_manager.rank_model,
        parallel_manager.COMM_couple.Get_rank(),
        parallel_manager.node_rank,
        "✔" if parallel_manager.filterpe else "✘"
    ]

    # Gather all rank data
    pe_data = COMM_WORLD.gather(rank_info, root=0)

    COMM_WORLD.Barrier()

    if rank == 0:
        print("\n[ICESEE] Parallel Execution Configuration:\n")
        headers = ["World Rank", "Filter Rank", "Model Task", "Model Rank", "Couple Rank", "Shared Rank", "Filter PE"]
        print(tabulate(pe_data, headers=headers, tablefmt="double_grid"))

    COMM_WORLD.Barrier()


def find_largest_divisor(Nens, size_world):
    """
    Finds the largest divisor of `Nens` that is less than or equal to `size_world`.
    
    Parameters:
        - Nens (int): Number of ensemble members.
        - size_world (int): Total number of MPI processes.

    Returns:
        - best_size (int): Largest divisor of `Nens` ≤ `size_world`.
    """
    best_size = 1  # Start with the minimum valid size
    for i in range(1, size_world + 1):
        if Nens % i == 0:
            best_size = i  # Update with the largest divisor found
    return best_size


# --- Main Execution (test parallel manager) ---
if __name__ == "__main__":
    Nens = 8  # Number of ensemble members
    n_modeltasks = 2  # Set to `None` for automatic determination
    global_shape = 1024  # Number of state variables

    parallel_manager = icesee_mpi_parallelization(Nens=Nens, global_shape=global_shape, n_modeltasks=n_modeltasks, screen_output=True)

    # Finalize MPI
    if MPI.Is_initialized():
        MPI.Finalize()
