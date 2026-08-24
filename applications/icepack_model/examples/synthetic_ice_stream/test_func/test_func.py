def  compute_forecast_mean_chunked_v2(self, k, flag=None):
    """
    Simple & hang-free:
    - running sum in RAM (length = local_rows)
    - collective dataset creation
    - collective write with empty selection for zero-row ranks
    """
    from mpi4py import MPI
    import numpy as np, h5py,shutil,os
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


def compute_forecast_mean_v3(self, k, flag=None):
    """
    Parallel ensemble mean for column k.
    Reads local (rows x nens) block in one shot and writes mean[:, k].
    Bitwise-equivalent to np.mean(block, axis=1, dtype=np.float64).
    """
    import numpy as np, h5py, shutil, sys, traceback
    from mpi4py import MPI
    import h5py.h5p as h5p, h5py.h5s as h5s, h5py.h5fd as h5fd

    comm = self.mpi_comm
    rank = comm.Get_rank()

    nd0, nd1 = self.nd_start_world, self.nd_end_world
    local_rows = max(nd1 - nd0, 0)
    nens = int(self.nens)
    nt   = int(self.nt)

    # Map global timestep to dataset index (if you batch timesteps)
    batch_idx = int(k - self.current_batch_start)

    try:
        # ---------- 1) Read (rows x nens) in one I/O ----------
        # Expect shape (nd, nens) for this timestep/batch.
        # If your storage is transposed, slice accordingly once here.
        # Ensure float64 to match numpy reference.
        if local_rows > 0:
            block = np.array(self.datasets[batch_idx][nd0:nd1, :nens], dtype=np.float64, copy=False)
            if block.shape != (local_rows, nens):
                raise ValueError(f"Expected (local_rows,nens)=({local_rows},{nens}), got {block.shape}")
            local_mean = block.mean(axis=1, dtype=np.float64)
        else:
            local_mean = np.empty((0,), dtype=np.float64)

        # ---------- 2) Prepare output file/dataset ----------
        file_path = f"{self.base_path}/{self.file_prefix}_mean.h5"

        if flag == "initial" and rank == 0 and k == 0:
            try:
                shutil.rmtree(file_path)
            except OSError:
                pass
        comm.Barrier()

        dxpl = h5p.create(h5p.DATASET_XFER)
        dxpl.set_dxpl_mpio(h5fd.MPIO_COLLECTIVE)

        with h5py.File(file_path, "a", driver="mpio", comm=comm) as f:
            # Single, synchronized branch for creation
            exists = ("mean" in f)
            exists_any = comm.bcast(exists if rank == 0 else None, root=0)
            if not exists_any:
                if rank == 0:
                    chunk_rows = min(self.nd, 4096)
                    f.create_dataset(
                        "mean",
                        shape=(self.nd, nt),
                        chunks=(chunk_rows, 1),
                        dtype="f8",
                    )
                comm.Barrier()

            dset = f["mean"]

            # ---------- 3) Collective write of column k ----------
            file_space = dset.id.get_space()
            if local_rows > 0:
                file_space.select_hyperslab((nd0, k), (local_rows, 1))
                mem_space = h5s.create_simple((local_rows,))
                buf = np.ascontiguousarray(local_mean)
            else:
                # Empty selection for ranks owning no rows
                mem_space = h5s.create_simple((0,))
                file_space.select_none()
                buf = np.empty((0,), dtype=np.float64)

            dset.id.write(mem_space, file_space, buf, dxpl=dxpl)

        comm.Barrier()

    except Exception as e:
        print(f"[rank {rank}] compute_forecast_mean_v3 failed: {e}")
        print("".join(traceback.format_exception(*sys.exc_info())))
        comm.Abort(1)


def compute_forecast_mean_v3_(self, k, flag=None):
    """
    Simple & hang-free:
    - running sum in RAM (length = local_rows)
    - collective dataset creation
    - collective write with empty selection for zero-row ranks
    """
    from mpi4py import MPI
    import numpy as np, h5py, shutil,os
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


def compute_forecast_mean_v3__(self, k, flag=None):
    """
    Compute ensemble mean at timestep k across all ranks and write collectively.
    Produces identical results to np.mean(X[:, :, k], axis=1).
    """
    import numpy as np, h5py, shutil, sys, traceback
    from mpi4py import MPI
    import h5py.h5p as h5p, h5py.h5s as h5s, h5py.h5fd as h5fd

    comm = self.mpi_comm
    rank = comm.Get_rank()
    size = comm.Get_size()

    nt = self.nt
    nd0, nd1 = self.nd_start_world, self.nd_end_world
    local_rows = max(nd1 - nd0, 0)
    nens = self.nens

    try:
        # ---- Local running sum ----
        local_sum = np.zeros(local_rows, dtype=np.float64)
        batch_idx = k - self.current_batch_start

        for ens_idx in range(nens):
            v = np.array(self.datasets[batch_idx][nd0:nd1, ens_idx], dtype=np.float64)
            if v.shape != (local_rows,):
                raise ValueError(f"Dataset shape mismatch: {v.shape} expected {(local_rows,)}")
            local_sum += v

        local_mean = local_sum / float(nens) if local_rows > 0 else np.empty((0,), dtype=np.float64)

        # ---- Parallel HDF5 write ----
        file_path = f"{self.base_path}/{self.file_prefix}_mean.h5"
        if flag == 'initial' and rank == 0 and k == 0:
            try: shutil.rmtree(file_path)
            except OSError: pass
        comm.Barrier()

        dxpl = h5p.create(h5p.DATASET_XFER)
        dxpl.set_dxpl_mpio(h5fd.MPIO_COLLECTIVE)

        with h5py.File(file_path, 'a', driver='mpio', comm=comm) as f:
            if 'mean' not in f:
                chunk_rows = min(self.nd, 4096)
                f.create_dataset('mean', (self.nd, self.nt),
                                chunks=(chunk_rows, 1), dtype='f8')
                comm.Barrier()
            dset = f['mean']

            file_space = dset.id.get_space()
            if local_rows > 0:
                file_space.select_hyperslab((nd0, k), (local_rows, 1))
                mem_space = h5s.create_simple((local_rows,))
                buf = np.ascontiguousarray(local_mean)
            else:
                mem_space = h5s.create_simple((1,))
                file_space.select_none()
                buf = np.empty((0,), dtype=np.float64)

            dset.id.write(mem_space, file_space, buf, dxpl=dxpl)

        comm.Barrier()

    except Exception as e:
        print(f"[Rank {rank}] Error in compute_forecast_mean_chunked_v2: {e}")
        print("".join(traceback.format_exception(*sys.exc_info())))
        comm.Abort(1)