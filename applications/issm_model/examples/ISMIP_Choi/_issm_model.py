# ==============================================================================
# @des: This file contains run functions for ISSM model python wrapper.
#       - contains different options of the EnKF data assimilation schemes.
# @date: 2025-03-26
# @author: Brian Kyanjo
# ==============================================================================

# --- python imports ---
import sys
import os
import shutil
import numpy as np
from scipy.stats import multivariate_normal,norm

# --- Utility imports ---
from ICESEE.config._utility_imports import icesee_get_index
from ICESEE.applications.issm_model.issm_utils.matlab2python.mat2py_utils import setup_reference_data, setup_ensemble_data
from ICESEE.applications.issm_model.issm_utils.matlab2python.server_utils import run_icesee_with_server

# --- model initialization ---
def initialize_model(**kwargs):
    """ des: intialize the issm model
        - calls the issm initalize_model.m matlab function to initialize the model
    """
    import h5py
    import scipy.io as sio

    # --- copy intialize_model.m to the current directory
    shutil.copyfile(os.path.join(os.path.dirname(__file__), 'initialize_model.m'), 'initialize_model.m')

    # -- get parameters from kwargs
    comm = kwargs.get('icesee_comm')
    icesee_rank = comm.Get_rank()
    icesee_size = kwargs.get('model_nprocs') 
    ens_id      = kwargs.get('ens_id')
    server      = kwargs.get('server')
    icesee_path = kwargs.get('icesee_path')
    data_path   = kwargs.get('data_path')
    vec_inputs  = kwargs.get('vec_inputs')
    use_reference_data = kwargs.get('use_reference_data', False)
    reference_data_dir = kwargs.get('reference_data_dir')
    reference_data     = kwargs.get('reference_data')

    # --- prepare the reference data if use_reference_data is True ---
    rank_data_dir, rank_data_file = setup_reference_data(reference_data_dir, reference_data, use_reference_data)
    
    #  call the issm initalize_model.m matlab function to initialize the model
    issm_cmd = f"run(\'issm_env\'); initialize_model({icesee_rank}, {icesee_size}, {ens_id})"
    # result = run_icesee_with_server(lambda: server.send_command(issm_cmd),server,False,comm)
    if not server.send_command(issm_cmd):
        print("[ICESEE DEBUG] Error sending command: {issm_cmd}")
        server.kill_matlab_processes()
        sys.exit(1)       
    

    # -- we would have broadcasted data to the remaining  ranks but now if nprocs > Nens, we need to duplicate data by copying data from ens_id_0000 to ens_id_0001, ens_id_0002, ... ens_id_000Nens
    Nens = kwargs.get('Nens')
    ensemble_dir, ensemble_kwargs = setup_ensemble_data(Nens)

    # fetch model size from output file
    output_filename = f'{icesee_path}/{data_path}/ensemble_init_{ens_id}.h5'
    # print("[ICESEE DEBUG] Attempting to open file: {output_filename}")
    if not os.path.exists(output_filename):
        print("[ICESEE ERROR] File does not exist: {output_filename}")
        return None
    # --get the size of the state vector from the output file
    with h5py.File(output_filename, 'r', driver='mpio', comm=comm) as f:
        for key in vec_inputs:
            nd = f[key][0].shape[0]
        return nd

    
# ---- ISSM model ----
def ISSM_model(**kwargs):
    """ des: run the issm model
        - calls the issm run_model.m matlab function to run the model
    """

    # --- get the number of processors ---
    # nprocs = kwargs.get('nprocs')
    k  = kwargs.get('k')
    dt = kwargs.get('dt')
    tinitial = kwargs.get('tinitial')
    tfinal = kwargs.get('tfinal')
    ens_id = kwargs.get('ens_id')
    comm = kwargs.get('comm')

    # get rank
    rank   = comm.Get_rank()
    nprocs = kwargs.get('model_nprocs')

    # --- copy run_model.m to the current directory
    shutil.copyfile(os.path.join(os.path.dirname(__file__), 'run_model.m'), 'run_model.m')

    
    # --- call the run_model.m function ---
    server   = kwargs.get('server')
    filename = kwargs.get('fname') 

    try:
        cmd = (
            f"run('issm_env'); run_model('{filename}', {ens_id}, {rank}, {nprocs}, {k}, {dt}, {tinitial}, {tfinal}); "
        )
        if not server.send_command(cmd):
            print("[ICESEE DEBUG] Error sending command: {cmd}")
    except Exception as e:
        print("[ICESEE DEBUG] Error sending command: {e}")
        server.shutdown()
        server.reset_terminal()
        sys.exit(1)


# ---- Run model for ISSM ----
def run_model(ensemble, **kwargs):
    """
    Run the ISSM model with an ensemble matrix from ICESEE.

    Args:
        ensemble: Ensemble matrix for the model.
        **kwargs: Additional parameters including:
            - nprocs: Number of processors.
            - server: Server information.
            - issm_examples_dir: Directory for ISSM examples.
            - icesee_path: Path to ICESEE directory.
            - comm: MPI communicator object.
            - vec_inputs: List of input vector keys.
            - ens_id: Ensemble ID.
            - data_path: Path to data directory.
            - k: timestep index.

    Returns:
        dict: Dictionary containing the output from the ISSM model, or None if an error occurs.
    """
    import h5py
    import numpy as np
    import os

    # Extract keyword arguments
    nprocs = kwargs.get('nprocs')
    server = kwargs.get('server')
    issm_examples_dir = kwargs.get('issm_examples_dir')
    icesee_path = kwargs.get('icesee_path')
    comm = kwargs.get('comm')
    vec_inputs = kwargs.get('vec_inputs')
    ens_id = kwargs.get('ens_id')
    data_path = kwargs.get('data_path')

    # Change to ISSM examples directory
    os.chdir(issm_examples_dir)

    # Define filename for data saving
    fname = 'enkf_state.mat'
    kwargs.update({'fname': fname})

    # Generate output filename based on ensemble ID
    input_filename = f'{icesee_path}/{data_path}/ensemble_output_{ens_id}.h5'

    # Get ensemble indices
    vecs, indx_map, _ = icesee_get_index(ensemble, **kwargs)

    #  --- Joint Estimations ---
    if kwargs["joint_estimation"]:
        bed = ensemble[indx_map['bed']]
        coefficient = ensemble[indx_map['coefficient']]
    else:
        k = kwargs.get('k', 0) 
        if k == 0:
            bed_int = ensemble[indx_map['bed']]
            coefficient_int = ensemble[indx_map['coefficient']]
        bed = bed_int
        coefficient = coefficient_int

    # Write ensemble data to HDF5 file to be accessed by ISSM on the Matlab side
    with h5py.File(input_filename, 'w', driver='mpio', comm=comm) as f:
        # for key in vec_inputs:
        #     f.create_dataset(key, data=ensemble[indx_map[key]])
        f.create_dataset('Thickness', data=ensemble[indx_map['Thickness']])
        f.create_dataset('bed', data=bed)
        f.create_dataset('coefficient', data=coefficient)

    # Run ISSM model to update state and parameters
    try:
        ISSM_model(**kwargs)
    except Exception as e:
        print(f"[ICESEE run_model Error] Error running the ISSM model: {e}")
        server.kill_matlab_processes()
        return None

    # Read output from HDF5 file to be accessed by ICESEE on the Python side
    output_filename = f'{icesee_path}/{data_path}/ensemble_output_{ens_id}.h5'
    if not os.path.exists(output_filename):
        print("[ICESEE run_model Error] File does not exist: {output_filename}")
        return None
    
    updated_state = {}
    with h5py.File(output_filename, 'r', driver='mpio', comm=comm) as f:
        updated_state['Thickness'] = f['Thickness'][0]
        # --Joint Estimations--
        if kwargs["joint_estimation"]:
            updated_state['bed'] = f['bed'][0]
            updated_state['coefficient'] = f['coefficient'][0]
        else:
            bed_int = bed
            coefficient_int = coefficient
            
    os.chdir(icesee_path)

    return updated_state
