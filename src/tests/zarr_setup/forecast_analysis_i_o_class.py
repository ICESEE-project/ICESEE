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
            self.datasets[batch_idx][self.nd_start:self.nd_end, ens_idx] = data
            write_time = MPI.Wtime() - start
        except Exception as e:
            print(f"Error occurred in write_analysis: {e}")
            tb_str = "".join(traceback.format_exception(*sys.exc_info()))
            print(f"Traceback details:\n{tb_str}")
            self.mpi_comm.Abort(1)

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

    def compute_forecast_mean_chunked_(self, t, ens_chunk_size=None, use_collective_io=True, max_ranks=20):
        self._ensure_batch(t)
        comm = self.mpi_comm
        rank = comm.Get_rank()
        size = comm.Get_size()

        import numpy as np
        import h5py.h5p as h5p
        import h5py.h5s as h5s
        import h5py.h5fd as h5fd

        try:
            # Use provided max_ranks, with fallback to serial for very small sizes
            if size <= 2:
                active_ranks = 1
            else:
                active_ranks = min(max_ranks, size)
                # Optional: Add min_rows_per_rank logic if nd is known
                # e.g., active_ranks = min(max_ranks, max(1, self.nd // min_rows_per_rank))

            # Exclude ranks with no work
            color = 0 if (rank < active_ranks) else MPI.UNDEFINED
            sub_comm = comm.Split(color, rank)
            active = sub_comm != MPI.COMM_NULL

            # Timing instrumentation
            t_io = 0.0
            t_comp = 0.0
            t_comm = 0.0
            start_total = MPI.Wtime()

            if active:
                sub_rank = sub_comm.Get_rank()
                sub_size = sub_comm.Get_size()

                # Distribute rows
                local_rows = (self.nd_end_world - self.nd_start_world) if rank == 0 else 0
                local_rows = sub_comm.bcast(local_rows, root=0)
                rows_per_rank = local_rows // sub_size
                extra_rows = local_rows % sub_size
                nd_start = sub_rank * rows_per_rank + min(sub_rank, extra_rows)
                nd_end = nd_start + rows_per_rank + (1 if sub_rank < extra_rows else 0)

                batch_idx = t - self.current_batch_start
                ds = self.datasets[batch_idx]

                # Dynamic chunk size
                if ens_chunk_size is None:
                    bytes_per_element = 8
                    target_memory = 1e9  # 1GB target
                    ens_chunk_size = max(1, int(target_memory / (nd_end - nd_start) / bytes_per_element))
                    ens_chunk_size = min(ens_chunk_size, self.nens)
                    if sub_rank == 0:
                        print(f"Dynamic ens_chunk_size: {ens_chunk_size}")

                local_sum = np.zeros(nd_end - nd_start, dtype='f8')

                # Set up I/O
                dxpl = h5p.create(h5p.DATASET_XFER)
                if size <= 2:
                    dxpl.set_dxpl_mpio(h5fd.MPIO_INDEPENDENT)  # Serial I/O for small sizes
                elif use_collective_io:
                    dxpl.set_dxpl_mpio(h5fd.MPIO_COLLECTIVE)
                else:
                    dxpl.set_dxpl_mpio(h5fd.MPIO_INDEPENDENT)


                # Double buffering
                buffers = [np.empty((nd_end - nd_start, ens_chunk_size), dtype='f8') for _ in range(2)]
                current_buffer = 0

                for start_ens in range(0, self.nens, ens_chunk_size):
                    end_ens = min(start_ens + ens_chunk_size, self.nens)
                    chunk_cols = end_ens - start_ens

                    if chunk_cols < ens_chunk_size:
                        buffers[current_buffer] = np.empty((nd_end - nd_start, chunk_cols), dtype='f8')

                    # I/O
                    t_start_io = MPI.Wtime()
                    file_space = ds.id.get_space()
                    file_space.select_hyperslab((self.nd_start_world + nd_start, start_ens), (nd_end - nd_start, chunk_cols))
                    mem_space = h5s.create_simple((nd_end - nd_start, chunk_cols))
                    ds.id.read(mem_space, file_space, buffers[current_buffer], dxpl=dxpl)
                    t_io += MPI.Wtime() - t_start_io

                    # Compute on previous buffer
                    if start_ens > 0:
                        t_start_comp = MPI.Wtime()
                        local_sum += np.sum(buffers[1 - current_buffer], axis=1)
                        t_comp += MPI.Wtime() - t_start_comp

                    current_buffer = 1 - current_buffer

                # Final compute
                t_start_comp = MPI.Wtime()
                local_sum += np.sum(buffers[current_buffer], axis=1)
                t_comp += MPI.Wtime() - t_start_comp

                local_mean = local_sum / self.nens

                # Output
                file_path = f"{self.base_path}/{self.file_prefix}_mean.h5"
                t_start_comm = MPI.Wtime()
                sub_comm.Barrier()
                t_comm += MPI.Wtime() - t_start_comm

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

                t_start_comm = MPI.Wtime()
                sub_comm.Barrier()
                t_comm += MPI.Wtime() - t_start_comm

                if sub_rank == 0:
                    print(f"Total time: {MPI.Wtime() - start_total:.2f}s, "
                        f"I/O: {t_io:.2f}s, Compute: {t_comp:.2f}s, Comm: {t_comm:.2f}s")

                # Note: For async I/O, explore h5py.AsyncRequest with mpi4py futures in Python 3.12+
                # Not implemented here due to experimental status

            t_start_comm = MPI.Wtime()
            comm.Barrier()
            t_comm += MPI.Wtime() - t_start_comm

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

    def _create_synthetic_observations(self, **kwargs):
        try:
            synthetic_obs_zarr_path = kwargs.get('synthetic_obs_zarr_path')
            error_R_zarr_path = kwargs.get('error_R_zarr_path')
            nd = kwargs.get('nd')
            nt = kwargs.get('nt')

            obs_t, ind_m, m_obs = self.generate_observation_schedule(**kwargs)
            m = m_obs
            m_R = m_obs*2 +1

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

    def Eta_matrix(self, k, **kwargs):
        try:
            H_matrix_zarr_path = kwargs.get('H_matrix_zarr_path', "output/H_matrix.zarr")
            Eta_matrix_zarr_path = kwargs.get('Eta_matrix_zarr_path', "output/Eta_matrix.zarr")

            H_matrix_file = zarr.open_array(H_matrix_zarr_path, mode='r')
            H_local = H_matrix_file[:, self.nd_start_world:self.nd_end_world]

            mean_file_path = f"{self.base_path}/{self.file_prefix}_mean.h5"
            with h5py.File(mean_file_path, 'r', driver='mpio', comm=self.mpi_comm) as f:
                ens_mean = f['mean'][self.nd_start_world:self.nd_end_world, k]
                
            ens_idx = 1
            state = self.read_analysis(k, ens_idx)
            ens_pertubations = state - ens_mean
            Eta_local = np.dot(H_local, ens_pertubations)
        except Exception as e:
            print(f"Error in Eta_matrix: {e}")
            tb_str = "".join(traceback.format_exception(*sys.exc_info()))
            print(f"Traceback details:\n{tb_str}")
            self.mpi_comm.Abort(1)

    def d_matrix(self, km, **kwargs):
        try:
            H_matrix_zarr_path = kwargs.get('H_matrix_zarr_path', "output/H_matrix.zarr")
            synthetic_obs_zarr_path = kwargs.get('synthetic_obs_zarr_path', "output/synthetic_obs.zarr")

            H_matrix = zarr.open_array(H_matrix_zarr_path, mode='r')
            H_local = H_matrix[:, self.nd_start_world:self.nd_end_world]

            #  get the synthetic observations
            synthetic_obs = zarr.open_array(synthetic_obs_zarr_path, mode='r')
            synthetic_obs_local = synthetic_obs[self.nd_start_world:self.nd_end_world, :]
            print(f"\nShape of synthetic observations: {synthetic_obs_local.shape} H_local shape: {H_local.shape}\n")

            #TODO -> will start from here after lunch


            #  compute the d matrix
            # d_local = np.dot(H_local, synthetic_obs_local)
        except Exception as e:
            print(f"Error in d_matrix: {e}")
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