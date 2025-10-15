# H5py Spack Environment

Spack environment for HDF5 and h5py with parallel I/O support.

## Purpose

This environment provides HDF5 with MPI support and the h5py Python interface, which are critical for ICESEE's parallel I/O operations.

## Contents

The `spack.yaml` file specifies:
- HDF5 with MPI support
- h5py Python package
- Compatible MPI implementation
- Required dependencies

## Features

### Parallel HDF5
- MPI-enabled for parallel I/O
- Optimized for HPC filesystems
- Thread-safe configuration

### h5py
- Python interface to HDF5
- Parallel I/O capabilities
- NumPy integration

## Installation

### Activate and Install
```bash
cd config/spack_env/h5py-env
spack env activate .
spack install
```

### Verify Installation
```bash
python -c "import h5py; print(h5py.version.info)"
python -c "import h5py; print('Parallel:', h5py.get_config().mpi)"
```

## Usage

After installation, the environment provides:
- HDF5 libraries with MPI support
- h5py Python module
- Proper library paths and environment variables

### In Python
```python
import h5py

# Create parallel HDF5 file
with h5py.File('data.h5', 'w', driver='mpio', comm=MPI.COMM_WORLD) as f:
    dset = f.create_dataset('data', shape=(1000, 1000), dtype='f')
    dset[rank*chunk:(rank+1)*chunk, :] = local_data
```

## Configuration Options

The Spack spec may include:
- `+mpi`: Enable MPI support (required)
- `+fortran`: Fortran interface (optional)
- `+hl`: High-level API (recommended)

## Troubleshooting

### MPI Not Found
Ensure MPI is available:
```bash
spack find mpi
```

If not installed:
```bash
spack install mpi
```

### h5py Import Error
Check Python can find h5py:
```bash
python -c "import h5py"
```

If it fails, ensure the environment is activated:
```bash
spack env activate .
```

## Performance Tips

For optimal parallel I/O:
- Use collective operations when possible
- Align chunk sizes with filesystem blocks
- Tune MPI-IO hints
- Use parallel HDF5 filters

## System-Specific Notes

Some HPC systems have optimized HDF5 builds. Consider using system modules:
```bash
module avail hdf5
```

If suitable modules exist, you may not need this Spack environment.

## Requirements

- Spack
- MPI implementation
- C compiler
- Python development headers

## References

- [HDF5 Documentation](https://www.hdfgroup.org/solutions/hdf5/)
- [h5py Documentation](https://docs.h5py.org/)
- [Parallel HDF5 Guide](https://docs.h5py.org/en/stable/mpi.html)
