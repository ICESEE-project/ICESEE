# Synthetic Ice Stream Example

This example showcases synthetic ice stream modeling with various ensemble Kalman filters (EnKF, DEnKF, EnTKF, EnRSKF). You can run the example using one of two methods:

1. **Interactive Notebook**: Use the [synthetic_ice_stream_da.ipynb](./synthetic_ice_stream_da.ipynb) notebook for an exploratory, hands-on experience.
2. **Terminal Execution** (Recommended for HPC): Use the [run_da_icepack.py](./run_da_icepack.py) script for optimized, high-performance runs.

---

## Initialization data

Data assimilation starts from the final state (history index 20000) of the
1000-year spin-up stored in `data/extended_beta1000yrs.h5` (~39.5 GB). By
default (`compact_initialization: true` in `params.yaml`) the model reads that
state from a compact checkpoint, `data/extended_beta1000yrs_compact_idx20000.h5`
(~29 MB): the mesh, `velocity`/`thickness`/`surface` at index 20000, and
`bed`, `grounded`, `floating`, `fluidity`, `extended_beta`, copied unchanged.

Build the compact checkpoint once from the spin-up history, from this
directory:

```bash
python tools/build_compact_initialization.py \
    --source data/extended_beta1000yrs.h5 \
    --dest data/extended_beta1000yrs_compact_idx20000.h5 \
    --idx 20000
```

The full spin-up history is needed only to build the compact checkpoint, or to
initialize directly from it with
`--compact_initialization=False --initFile=data/extended_beta1000yrs.h5`.

---

## Execution mode 3: checkpoints and resuming

Mode 3 writes restart checkpoints of the whole ensemble under
`<data_path>/_mode3_state_history/`. By default only the newest two step
checkpoints and the initial ensemble are kept, so disk use does not grow with
the number of timesteps. These are restart state, not a results archive.

- Keep every step's ensemble (old behaviour, ~`Nens` × state size per step):
  `--checkpoint_keep_last=0`
- Also keep every analysis step: `--checkpoint_keep_analysis=True`
- Write a checkpoint only every N steps: `--checkpoint_every=N`

To continue an interrupted run (crash, wall-time limit, full disk), repeat the
**original command unchanged** and add `--resume_from_checkpoint=True`. The run
continues after the newest valid checkpoint and keeps the original
observations. Without that flag a run starts by deleting `data_path`. Details:
`docs/execution-mode-3-design.md` ("Restart, retention, and crash consistency").

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