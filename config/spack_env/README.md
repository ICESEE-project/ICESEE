# Spack Environment Configurations

This directory contains Spack environment specifications for managing ICESEE dependencies.

## Purpose

Spack is a package manager designed for HPC systems that helps manage complex software dependencies. These environment files define reproducible software environments for ICESEE.

## Structure

Each subdirectory contains a Spack environment definition for specific dependency sets:
- `h5py-env/` - Environment for HDF5 and h5py with parallel support

## Using Spack Environments

### Install Spack
```bash
git clone https://github.com/spack/spack.git
source spack/share/spack/setup-env.sh
```

### Activate an Environment
```bash
cd config/spack_env/h5py-env
spack env activate .
```

### Install Dependencies
```bash
spack install
```

### Deactivate Environment
```bash
spack env deactivate
```

## Benefits

Using Spack environments provides:
- **Reproducibility**: Exact versions and configurations
- **Portability**: Works across different HPC systems
- **Customization**: Optimize builds for specific hardware
- **Dependency Management**: Automatic handling of complex dependencies

## Modifying Environments

To modify an environment:
1. Edit `spack.yaml` in the environment directory
2. Update package versions or add new packages
3. Reinstall: `spack install`

## Creating New Environments

For new dependency sets:
```bash
cd config/spack_env
spack env create new-env
cd new-env
# Edit spack.yaml
spack install
```

## System-Specific Notes

Some HPC systems have Spack already installed. Check with:
```bash
which spack
```

If available, use the system Spack instead of installing your own.

## Requirements

- Spack package manager
- C/C++/Fortran compilers
- Sufficient disk space for builds

## Resources

- [Spack Documentation](https://spack.readthedocs.io/)
- [Spack Environments Guide](https://spack.readthedocs.io/en/latest/environments.html)
- [ICESEE Wiki](https://github.com/ICESEE-project/ICESEE/wiki)
