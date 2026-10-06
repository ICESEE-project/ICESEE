from mpi4py import MPI
import h5py
import zarr
import numpy as np
import click
import os
import shutil
import time
import uuid

# from diagtimer import DiagnosticTimer
# diagtimer.py
from contextlib import contextmanager
from time import time
import pandas as pd

class DiagnosticTimer:
    def __init__(self):
        self.diagnostics = []

    @contextmanager
    def time(self, **kwargs):
        tic = time()
        yield
        toc = time()
        kwargs['runtime'] = toc - tic
        self.diagnostics.append(kwargs)

    def dataframe(self):
        return pd.DataFrame(self.diagnostics)

# Zarr API updates
# try:
#     from zarr.storage import DirectoryStore, NestedDirectoryStore
# except Exception:  # very old zarr
#     DirectoryStore = zarr.DirectoryStore
#     NestedDirectoryStore = zarr.NestedDirectoryStore

# Compressors now come from numcodecs
try:
    from numcodecs import GZip  # or Zlib; GZip matches your original intent
except Exception:
    GZip = None  # fall back to no compression if unavailable


@click.command()
@click.option('--nsteps', default=1, help="Number of iterations to perform")
@click.option('--size', default=1_000_000, help="Length of each row of array")
@click.option('--output_dir', default='./', help="Where to write the data")
@click.option('--compression', default='none', type=click.Choice(['none', 'gzip']))
@click.option('--nested', is_flag=True, help="Use Zarr NestedDirectoryStore")
def main(nsteps, size, output_dir, compression, nested):
    timer = DiagnosticTimer()

    comm = MPI.COMM_WORLD
    nprocs = comm.Get_size()
    rank = comm.Get_rank()

    x = np.arange(size)
    dtype = 'f8'
    shape = (nsteps, nprocs, size)

    # IMPORTANT: each (n, rank, :) is one chunk → no overlap across ranks
    chunks = (1, 1, size)
    chunk_size_bytes = np.dtype(dtype).itemsize * size

    # HDF5: filters still not supported for parallel writes; keep none
    hdf_compression_kw = {}

    # Zarr compression setup (current Zarr uses numcodecs compressor objects)
    if compression == 'none' or GZip is None:
        zarr_compressor = None
        zarr_compression_type, zarr_compression_level = None, None
    else:
        zarr_compressor = GZip(level=4)
        zarr_compression_type, zarr_compression_level = 'gzip', 4

    zarr_log_options = dict(
        nprocs=nprocs,
        size_in_bytes=chunk_size_bytes,
        format='zarr',
        compression=zarr_compression_type,
        compression_level=zarr_compression_level,
        nested=nested,
    )
    hdf_log_options = dict(
        nprocs=nprocs,
        size_in_bytes=chunk_size_bytes,
        format='hdf',
        compression=None,
        compression_level=None,
        nested=None,
    )

    # Unique filename base
    uid = comm.bcast(str(uuid.uuid1())[:8] if rank == 0 else None, root=0)
    fname_base = f'parallel_test_{uid}'

    # -------------------------------
    # HDF5 init (parallel) — unchanged
    # -------------------------------
    hfname = os.path.join(output_dir, f'{fname_base}.hdf5')
    hfile = h5py.File(hfname, 'w', driver='mpio', comm=MPI.COMM_WORLD)
    hdset = hfile.create_dataset('test', shape, dtype=dtype, chunks=chunks, **hdf_compression_kw)

    # -------------------------------
    # Zarr init (modern API)
    # -------------------------------
    import fsspec

    def get_store(path):
        # Local directory store via fsspec mapper (works for v2/v3)
        return fsspec.get_mapper(path)

    zname = os.path.join(output_dir, f'{fname_base}.zarr')
    store = get_store(zname)

    # Only rank 0 creates metadata/array, then others reopen in r+.
    if rank == 0:
        # zarr.open returns a zarr.Array (v2) or compatible (v3) when given shape/chunks
        z_array = zarr.open(
            store,
            mode='w',
            shape=shape,
            chunks=chunks,
            dtype=dtype,
            compressor=zarr_compressor,
        )
    comm.Barrier()
    # All ranks reopen for writing
    z_array = zarr.open(store, mode='r+')

    # --------------------------------
    # Write loop (each rank writes its chunk)
    # --------------------------------
    rng = np.random.default_rng(seed=rank + 12345)
    for n in range(nsteps):
        # Make deterministic-ish data per rank
        phase = 2 * np.pi * rng.random()
        data = (np.cos(20 * np.pi * x / size + phase) + 0.1 * rng.random(size)).astype(dtype)

        # Zarr write (no overlapping chunks across ranks)
        with timer.time(operation='write', step=n, **zarr_log_options):
            z_array[n, rank, :] = data
            comm.Barrier()  # keep steps aligned for cleaner timing

        # HDF5 write
        with timer.time(operation='write', step=n, **hdf_log_options):
            hdset[n, rank, :] = data
            comm.Barrier()

    hfile.close()

    # --------------------------------
    # Readback & verify
    # --------------------------------
    hfile = h5py.File(hfname, 'r', driver='mpio', comm=MPI.COMM_WORLD)
    hdset = hfile['test']
    z_array = zarr.open(store, mode='r')

    if rank == 0:
        print('format,nprocs,operation,runtime,rank')

    for n in range(nsteps):
        with timer.time(operation='read', step=n, **zarr_log_options):
            tic = time.time()
            zdata = np.array(z_array[n, rank, :], copy=False)
            readtime = time.time() - tic
            comm.Barrier()
        print(f'zarr,{nprocs},read,{readtime},{rank}')

        with timer.time(operation='read', step=n, **hdf_log_options):
            tic = time.time()
            hdata = np.array(hdset[n, rank, :], copy=False)
            readtime = time.time() - tic
            comm.Barrier()
        print(f'hdf,{nprocs},read,{readtime},{rank}')

        np.testing.assert_allclose(zdata, hdata)

    hfile.close()

    # Cleanup & export timings
    if rank == 0:
        shutil.rmtree(zname, ignore_errors=True)
        try:
            os.remove(hfname)
        except FileNotFoundError:
            pass
        df = timer.dataframe()
        df.to_csv(f'parallel_read_write_{time.strftime("%Y-%m-%d_%H%M.%S")}.csv', index=False)


if __name__ == "__main__":
    main()