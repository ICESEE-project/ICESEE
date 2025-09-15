# Synthetic Ice Stream Example

This example showcases synthetic ice stream modeling with various ensemble Kalman filters (EnKF, DEnKF, EnTKF, EnRSKF). You can run the example using one of two methods:

1. **Interactive Notebook**: Use the [synthetic_ice_stream_da.ipynb](./synthetic_ice_stream_da.ipynb) notebook for an exploratory, hands-on experience.
2. **Terminal Execution** (Recommended for HPC): Use the [run_da_icepack.py](./run_da_icepack.py) script for optimized, high-performance runs.

---

## Running via `run_da_icepack.py`

Follow these steps to execute the [run_da_icepack.py](./run_da_icepack.py) script. All parameters are managed in the [params.yaml](./params.yaml) file for easy configuration.

### Steps

1. **Set Up Parameters**:
   - Modify the [params.yaml](./params.yaml) file to specify your desired inputs and parameters.
   - The script retrieves these parameters using helper functions from the [config/_utility_imports](https://github.com/KYANJO/ICESEE/blob/main/config/_utility_imports.py) module.

2. **Run the Script**:
   - **Serial Execution**:
     ```bash
     python run_da_icepack.py
     ```
   - **Parallel Execution**:
     ```bash
     mpiexec -n 8 python run_da_icepack.py
     ```

3. **Select a Filter**:
   - **Note**: Filter selection is only available in serial mode. Parallel mode currently supports only `EnKF`.
   - Update the `filter_type` parameter in [params.yaml](./params.yaml) to choose a filter:
     - `EnKF`: Ensemble Kalman Filter
     - `DEnKF`: Deterministic Ensemble Kalman Filter
     - `EnTKF`: Ensemble Transform Kalman Filter
     - `EnRSKF`: Ensemble Square Root Kalman Filter

4. **View Outputs**:
   - Results are stored as `.h5` files in the `results` and `_modelrun_datasets` directories, named as:
     ```
     filter_type-model.h5
     ```

5. **Analyze Results**:
   - Use the [read_results.ipynb](./read_results.ipynb) notebook to load and visualize the results.

---

## Running with Containers

For containerized environments (e.g., HPC clusters), use Apptainer/Singularity to run the script.

### Steps

1. **Build the Container**:
   - Follow the instructions in the [/src/container/apptainer/](https://github.com/KYANJO/ICESEE/tree/main/src/container/apptainer) directory to build the `icepack.sif` container image.

2. **Execute the Script**:
   - Run the script within the container:
     ```bash
     apptainer exec icepack.sif python run_da_icepack.py
     ```