## **Lorenz96 model**

This example demonstrates Lorenz96 modeling using various types of ensemble Kalman filters (EnKF, DEnKF, EnTKF, and EnRSKF). You can run the example using one of two methods:

1. **Notebook Interface**: Use the interactive notebook [EnKF_all_types.ipynb](./EnKF_all_types.ipynb) for a hands-on, exploratory approach.
2. **Terminal Execution** (Recommended for HPC platforms): Use the [run_da_lorenz96.py](./run_da_lorenz96.py) script for high-performance and feature-rich execution.

---

### **Running with [run_da_lorenz96.py](./run_da_lorenz96.py)**

To execute the [run_da_lorenz96.py](./run_da_lorenz96.py) script, follow these steps:

**Note**: All parameters are centrally managed in the [params.yaml](./params.yaml) file for ease of modification and experimentation.

1. **Configure Inputs**:
   - Edit the [params.yaml](./params.yaml) file to define your desired inputs and parameters.
   - The script dynamically fetches these parameters using helper functions from the [config/icepack_config](./../../../config/icepack_config) directory.

2. **Run the Script**:
   ```bash
   python run_da_lorenz96.py
   ```

3. **Experiment with Filters**:
   - Modify the `filter_type` parameter in the [params.yaml](./params.yaml) file to switch between different filters:
     - **EnKF**: Ensemble Kalman Filter
     - **DEnKF**: Deterministic Ensemble Kalman Filter
     - **EnTKF**: Ensemble Transform Kalman Filter
     - **EnRSKF**: Ensemble Square Root Kalman Filter

4. **Outputs**:
   - Results are saved in the `results/` directory as `.h5` files with the naming format:
     ```
     filter_type-model.h5
     ```

5. **Visualize Results**:
   - Use the notebook [read_results.ipynb](./read_results.ipynb) to load and plot the results for detailed analysis.

---
