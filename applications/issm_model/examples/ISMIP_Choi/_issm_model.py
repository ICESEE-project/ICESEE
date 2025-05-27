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
from ICESEE.applications.issm_model.issm_utils.matlab2python.mat2py_utils import subprocess_cmd_run
from ICESEE.applications.issm_model.issm_utils.matlab2python.server_utils import run_icesee_with_server

# --- model initialization ---
def initialize_model(physical_params, modeling_params, comm):
    """ des: intialize the issm model
        - calls the issm initalize_model.m matlab function to initialize the model
    """
    import h5py
    import scipy.io as sio

    # --- copy intialize_model.m to the current directory
    shutil.copyfile(os.path.join(os.path.dirname(__file__), 'initialize_model.m'), 'initialize_model.m')

    icesee_rank = comm.Get_rank()
    # icesee_size = comm.Get_size()
    icesee_size = modeling_params.get('model_nprocs')
    ens_id      = modeling_params.get('ens_id')
    

    # read the model kwargs from the file
    server      = modeling_params.get('server')
    icesee_path = modeling_params.get('icesee_path')
    data_path   = modeling_params.get('data_path')
    issm_cmd = f"run(\'issm_env\'); initialize_model({icesee_rank}, {icesee_size}, {ens_id})"
    # result = run_icesee_with_server(lambda: server.send_command(issm_cmd),server,False,comm)
    if not server.send_command(issm_cmd):
        print(f"[DEBUG] Error sending command: {issm_cmd}")
        server.kill_matlab_processes()
        sys.exit(1)
    
    # if not result:
    #     sys.exit(1)
    # -- we would have broadcasted data to the remaining  ranks but now if nprocs > Nens, we need to duplicate data by copying data from ens_id_0000 to ens_id_0001, ens_id_0002, ... ens_id_000Nens
    if icesee_rank == 0:
        data_dir = './Models/ens_id_0'
        kwargs_data = 'model_kwargs_0.mat'
        Nens = modeling_params.get('Nens')
        for ens in range(1, Nens):
            new_data_dir = f'./Models/ens_id_{ens}'
            new_kwargs_data = f'model_kwargs_{ens}.mat'
            # if os.path.exists(new_data_dir):
            #     shutil.rmtree(new_data_dir)
            shutil.copytree(data_dir, new_data_dir, dirs_exist_ok=True)
            shutil.copyfile(kwargs_data, new_kwargs_data)
            # print(f"[DEBUG] Copied {data_dir} to {new_data_dir}")
            # print(f"[DEBUG] Copied {kwargs_data} to {new_kwargs_data}")
    comm.Barrier()

    # fetch model size from output file
    try: 
        output_filename = f'{icesee_path}/{data_path}/ensemble_init_{ens_id}.h5'
        # print(f"[DEBUG] Attempting to open file: {output_filename}")
        if not os.path.exists(output_filename):
            print(f"[ERROR] File does not exist: {output_filename}")
            return None
        with h5py.File(output_filename, 'r', driver='mpio', comm=comm) as f:
            Vx = f['Vx'][0]
            # get the size of the Vx variable
            nd = Vx.shape[0]
            # print(f"[DEBUG] Vx shape: {Vx.shape}")
            return nd
    except Exception as e:
        print(f"[Initialize-model] Error reading the file: {e}")

        
    
# ---- ISSM model ----
def ISSM_model(**kwargs):
    """ des: run the issm model
        - calls the issm run_model.m matlab function to run the model
    """

    # --- get the number of processors ---
    # nprocs = kwargs.get('nprocs')
    k = kwargs.get('k')
    dt = kwargs.get('dt')
    tinitial = kwargs.get('tinitial')
    tfinal = kwargs.get('tfinal')
    # rank = kwargs.get('rank')
    ens_id = kwargs.get('ens_id')
    comm = kwargs.get('comm')
    params = kwargs.get('params')

    # get rank
    rank   = comm.Get_rank()
    # nprocs = comm.Get_size()
    nprocs = params.get('model_nprocs')

    # print(f"[DEBUG] ISSM model {rank} of {nprocs}")

    # --- copy run_model.m to the current directory
    shutil.copyfile(os.path.join(os.path.dirname(__file__), 'run_model.m'), 'run_model.m')

    
    # --- call the run_model.m function ---
    server   = kwargs.get('server')
    filename = kwargs.get('fname') 
    # print(f"file name: {filename}")

    try:
        cmd = (
            f"run('issm_env'); run_model('{filename}', {ens_id}, {rank}, {nprocs}, {k}, {dt}, {tinitial}, {tfinal}); "
        )
        if not server.send_command(cmd):
            print(f"[DEBUG] Error sending command: {cmd}")
    except Exception as e:
        print(f"[DEBUG] Error sending command: {e}")
        server.shutdown()
        server.reset_terminal()
        sys.exit(1)


# ---- Run model for ISSM ----
def run_model(ensemble, **kwargs):
    """ des: run the issm model with ensemble matrix from ICESEE
        returns: dictionary of the output from the issm model
    """
    import h5py
    import numpy as np
    import os

    # --- get the number of processors ---
    nprocs = kwargs.get('nprocs')
    server              = kwargs.get('server')
    # rank                = kwargs.get('rank')
    issm_examples_dir   = kwargs.get('issm_examples_dir')
    icesee_path         = kwargs.get('icesee_path')
    comm                = kwargs.get('comm')
    vec_inputs          = kwargs.get('vec_inputs')

    rank                = comm.Get_rank()

    #  --- change directory to the issm directory ---
    os.chdir(issm_examples_dir)

    # --- filename for data saving
    fname = 'enkf_state.mat'
    kwargs.update({'fname': fname})
    
    ens_id = kwargs.get('ens_id')

    # try: 
    if True:
        # Generate output filename based on rank
        icesee_path    = kwargs.get('icesee_path')
        data_path      = kwargs.get('data_path')
        input_filename = f'{icesee_path}/{data_path}/ensemble_output_{ens_id}.h5'

        #  write the ensemble to the h5 file
        k = kwargs.get('k')

        # -- call the icess_get_index function to get the index of the ensemble
        # print(f"[DEBUG] run_model {rank} of {comm.Get_size()}")
        # print(f"[DEBUG] Ensemble: {ensemble[:5]}")
        vecs, indx_map, _ = icesee_get_index(ensemble, **kwargs)

        # if k > 0:
        if True:
            # -- create our ensemble for test purposes
             with h5py.File(input_filename, 'w', driver='mpio', comm=comm) as f:
                for key in vec_inputs:
                    f.create_dataset(key, data=ensemble[indx_map[key]])

        # --- call the issm model  to update the state and parameters variables ---
        ISSM_model(**kwargs)

        # -- change directory back to the original directory
        os.chdir(icesee_path)
        
        #  --- read the output from the h5 file ISSM model ---
        try:
        # if True:
            output_filename = f'{icesee_path}/{data_path}/ensemble_output_{ens_id}.h5'
            # check if the file exists
            if not os.path.exists(output_filename):
                print(f"[ERROR] File does not exist: {output_filename}")
                return None
            updated_state = {}
            with h5py.File(output_filename, 'r',driver='mpio',comm=comm) as f:
                for key in vec_inputs:
                    updated_state[key] = f[key][0]

        except Exception as e:
            print(f"[Run model] Error reading the file: {e}")
        #     return None
    
    # -- change directory back to the original directory
    os.chdir(icesee_path)
    return updated_state
