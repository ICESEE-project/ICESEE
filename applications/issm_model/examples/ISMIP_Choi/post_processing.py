
import numpy as np
import h5py
import glob
import os
import re
import zarr
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed

def load_single_file(file_path, dataset_name):
    """
    Load data from a single HDF5 file.
    
    Parameters:
    - file_path: Path to the HDF5 file.
    - dataset_name: Name of the dataset within the HDF5 file.
    
    Returns:
    - data: NumPy array from the dataset.
    """
    with h5py.File(file_path, 'r') as f:
        if dataset_name not in f:
            raise KeyError(f"Dataset '{dataset_name}' not found in {file_path}. Available datasets: {list(f.keys())}")
        return f[dataset_name][:]

def load_ensemble_data(folder_path, file_pattern='icesee_enkf_ens_*.h5', dataset_name='states', max_workers=None):
    """
    Load ensemble data from HDF5 files in parallel and create a (nd, nens, nt) dataset.
    
    Parameters:
    - folder_path: Directory containing the HDF5 files.
    - file_pattern: Pattern to match HDF5 files (default: 'icesee_enkf_0*.h5').
    - dataset_name: Name of the dataset within each HDF5 file (default: 'states').
    - max_workers: Number of parallel workers (default: None, uses all available).
    
    Returns:
    - data: NumPy array of shape (nd, nens, nt).
    - time_indices: List of time indices extracted from filenames.
    """
    # Get list of files matching the pattern
    file_list = sorted(glob.glob(os.path.join(folder_path, file_pattern)),
                      key=lambda x: int(re.search(r'\d+', os.path.basename(x)).group()))
    
    if not file_list:
        raise ValueError(f"No files found matching pattern {file_pattern} in {folder_path}")
    
    # Extract time indices from filenames
    time_indices = [int(re.search(r'\d+', os.path.basename(f)).group()) for f in file_list]
    
    # Load the first file to determine dimensions
    first_file = file_list[0]
    try:
        first_data = load_single_file(first_file, dataset_name)
    except KeyError as e:
        print(str(e))
        return None, None
    
    nd, nens = first_data.shape
    
    # Initialize output array
    nt = len(file_list)
    # data = np.zeros((nd, nens, nt))
    data = zarr.zeros((nd, nens, nt), chunks=(nd, nens, 1), dtype=first_data.dtype, store=folder_path + '/ensemble_data.zarr', overwrite=True)

    # Load data[0] from first file
    data[:, :, 0] = first_data
    
    # Load remaining files in parallel
    if nt > 1:
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            # Submit tasks for remaining files (index 1 to nt-1)
            future_to_index = {
                executor.submit(load_single_file, file_path, dataset_name): i
                for i, file_path in enumerate(file_list[1:], 1)
            }
            
            # Collect results as they complete
            for future in as_completed(future_to_index):
                t = future_to_index[future]
                try:
                    loaded_data = future.result()
                    if loaded_data.shape != (nd, nens):
                        raise ValueError(f"Shape mismatch in {file_list[t]}: expected {(nd, nens)}, got {loaded_data.shape}")
                    data[:, :, t] = loaded_data
                except Exception as exc:
                    print(f'File {file_list[t]} generated an exception: {exc}')
                    return None, None
    
    return data, time_indices

def main():
    parser = argparse.ArgumentParser(description="Process ensemble data from HDF5 files into a (nd, nens, nt) dataset.")
    parser.add_argument('--dataset-name', type=str, default='states', help="Name of the dataset in HDF5 files (default: 'states')")
    parser.add_argument('--folder-path', type=str, default='_modelrun_datasets', help="Path to folder with HDF5 files")
    args = parser.parse_args()
    
    try:
        # Load data in parallel
        dataset, time_indices = load_ensemble_data(args.folder_path, dataset_name=args.dataset_name)
        
        if dataset is None or time_indices is None:
            print("Failed to load dataset. Please check the dataset name and file contents.")
            return
        
        # Print summary
        print(f"Loaded dataset with shape: {dataset.shape}")
        print(f"Time indices: {time_indices}")
        
        # Save the dataset to a new HDF5 file
        output_file = os.path.join(args.folder_path, 'icesee_ensemble_dataset.h5')
        with h5py.File(output_file, 'w') as f:
            f.create_dataset('ensemble_data', data=dataset)
            f.create_dataset('time_indices', data=time_indices)
        print(f"Dataset saved to {output_file}")

        # cleanup zarr store
        zarr_store_path = args.folder_path + '/ensemble_data.zarr'
        if os.path.exists(zarr_store_path):
            import shutil
            shutil.rmtree(zarr_store_path)
            print(f"Cleaned up temporary zarr store at {zarr_store_path}")
        
    except Exception as e:
        print(f"Error: {str(e)}")

if __name__ == "__main__":
    main()