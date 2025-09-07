from spack import *
import os
import subprocess

class PyVenvPackages(PythonPackage):
    """Custom package to create a Python virtual environment and install pip packages."""

    homepage = "https://example.com"
    url = "file:///storage/home/hcoda1/8/bkyanjo3/r-arobel3-0/h5py-env/dummy.tar.gz"
    version('0.1', sha256='YOUR_CHECKSUM_HERE', expand=False)  # Replace with actual checksum

    # Dependencies
    depends_on('python@3.11.8', type=('build', 'run'))
    depends_on('py-pip', type=('build', 'run'))
    depends_on('py-setuptools', type=('build', 'run'))
    depends_on('py-wheel', type=('build', 'run'))
    depends_on('openmpi@5.0.7', type=('build', 'run'))  # For mpi4py and bigmpi4py
    depends_on('hdf5@1.12.1+mpi', type=('build', 'run'))  # For py-h5py

    # Override do_stage to handle the dummy tarball
    def do_stage(self, mirror_only=False):
        source_path = self.url.replace('file://', '')
        stage_path = join_path(self.stage.path, 'spack-src', 'dummy.tar.gz')

        # Ensure the stage directory exists
        mkdirp(join_path(self.stage.path, 'spack-src'))

        # Ensure the dummy tarball exists
        if not os.path.exists(source_path):
            # Create a minimal dummy tarball if it doesn't exist
            dummy_dir = join_path(self.stage.path, 'dummy')
            mkdirp(dummy_dir)
            with open(join_path(dummy_dir, 'dummy.txt'), 'w') as f:
                f.write('dummy')
            tar = which('tar')
            tar('-czf', stage_path, '-C', dummy_dir, 'dummy.txt')
        else:
            copy(source_path, stage_path)

    def install(self, spec, prefix):
        # Create virtual environment in the project directory
        venv_dir = '/storage/home/hcoda1/8/bkyanjo3/r-arobel3-0/h5py-env/.venv'
        python = spec['python'].command

        if not os.path.exists(venv_dir):
            python('-m', 'venv', venv_dir)

        # Paths to pip and python in the virtual environment
        pip = join_path(venv_dir, 'bin', 'pip')
        python_bin = join_path(venv_dir, 'bin', 'python')

        # Set environment variables for MPI and HDF5
        env['LD_LIBRARY_PATH'] = ':'.join([
            spec['openmpi'].prefix.lib,
            spec['hdf5'].prefix.lib,
            os.environ.get('LD_LIBRARY_PATH', '')
        ])
        env['CPATH'] = ':'.join([
            spec['openmpi'].prefix.include,
            spec['hdf5'].prefix.include,
            os.environ.get('CPATH', '')
        ])

        # Upgrade pip, setuptools, and wheel
        subprocess.run([pip, 'install', '--upgrade', 'pip', 'setuptools', 'wheel'], check=True)

        # Install pip packages
        pip_packages = [
            'numpy',
            'scipy>=1.15.1',
            'jupyter>=1.1.1',
            'matplotlib>=3.10.0',
            'pandas>=2.2.3',
            'bigmpi4py==1.2.4',
            'jax==0.5.3',
            'jaxlib==0.5.3',
            'jaxopt==0.8.3',
            'jaxtyping==0.3.0',
            'ray',
            'cython',
            'dask',
            'psutil==6.1.1',
            'tqdm==4.67.1',
            'pyyaml',
            'nose',
            'gstools',
            'mpi4py',
            'zarr'
        ]

        subprocess.run([pip, 'install'] + pip_packages, check=True)

        # Mark installation as complete
        touch(join_path(venv_dir, '.installed'))

        # Create a minimal prefix to satisfy Spack
        mkdirp(prefix.bin)
        touch(join_path(prefix.bin, 'dummy'))
