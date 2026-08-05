
import h5py
import glob
import os

def inspect_hdf5_file(folder_path, file_pattern='icesee_enkf_ens_*.h5'):
    """
    Inspect the first HDF5 file matching the pattern and print its datasets.
    
    Parameters:
    - folder_path: Directory containing the HDF5 files.
    - file_pattern: Pattern to match HDF5 files (default: 'icesee_enkf_ens_*.h5').
    """
    file_list = sorted(glob.glob(os.path.join(folder_path, file_pattern)))
    
    if not file_list:
        print(f"No files found matching pattern {file_pattern} in {folder_path}")
        return
    
    first_file = file_list[0]
    print(f"Inspecting file: {first_file}")
    
    with h5py.File(first_file, 'r') as f:
        print("Available datasets in the file:")
        for key in f.keys():
            print(f" - {key} (shape: {f[key].shape}, dtype: {f[key].dtype})")

if __name__ == "__main__":
    folder_path = '_modelrun_datasets'
    inspect_hdf5_file(folder_path)