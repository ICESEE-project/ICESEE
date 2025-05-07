# ==============================================================================
# @des: This file contains run functions for ISSM data assimilation.
#       - contains different options of the EnKF data assimilation schemes.
# @date: 2025-03-25
# @author: Brian Kyanjo
# ==============================================================================

import os
import numpy as np
import h5py


# --- import utility functions ---
from ICESEE.applications.issm_model.examples.ISMIP._issm_model import *
from ICESEE.config._utility_imports import icesee_get_index

# --- Forecast step ---
def forecast_step_single(ensemble=None, **kwargs):
    """ensemble: packs the state variables and parameters of a single ensemble member
    Returns: ensemble: updated ensemble member
    """
    #  -- control time stepping   
    # k = kwargs.get('k')
    time = kwargs.get('t')
    k    = kwargs.get('k')
    
    kwargs.update({'tinitial': time[k], 'tfinal': time[k+1]})

    realization = run_model(ensemble, **kwargs)

    # remove the copied directories
    # if k == 0:
    #     comm = kwargs.get('comm')
    #     comm.Barrier()
    #     rank = comm.Get_rank()
    #     if rank == 0:
    #         issm_path = kwargs.get('issm_examples_dir')
    #         data_dir = './Models/ens_id_0000'
    #         Nens = kwargs.get('Nens')
    #         for ens in range(1, Nens):
    #             new_data_dir = f'{issm_path}/Models/ens_id_000{ens}'
    #             if os.path.exists(new_data_dir):
    #                 shutil.rmtree(new_data_dir)
    #                 print(f"[DEBUG] Removed {new_data_dir}") 
    #     else:
    #         pass

    #  call the run_model fun to push the state forward in time
    return realization

# --- generate true state ---
def generate_true_state(**kwargs):
    """des: generate the true state of the model
    Returns: true_state: the true state of the model
    """
    params = kwargs.get('params')
    time   = kwargs.get('t')
    server = kwargs.get('server')
    
    issm_examples_dir   = kwargs.get('issm_examples_dir')
    icesee_path         = kwargs.get('icesee_path')
    data_path           = kwargs.get('data_path')
    comm                = kwargs.get('comm')

    # get the rank of the current process
    rank = comm.Get_rank()

    #  --- change directory to the issm directory ---
    os.chdir(issm_examples_dir)

    # --- filename for data saving
    fname = 'true_state.mat'
    kwargs.update({'fname': fname})
    ens_id = kwargs.get('ens_id')

    try:
    # if True:
        # --- fetch treu state vector
        statevec_true = kwargs.get('statevec_true')

        # -- call the icesee_get_index function to get the index of the state vector
        vecs, indx_map, dim_per_proc = icesee_get_index(statevec_true, **kwargs)

        # -- fetch data from inital state
        try: 
            output_filename = f'{icesee_path}/{data_path}/ensemble_init_{ens_id}.h5'
            # print(f"[DEBUG] Attempting to open file: {output_filename}")
            if not os.path.exists(output_filename):
                print(f"[ERROR] File does not exist: {output_filename}")
                return None
            with h5py.File(output_filename, 'r', driver='mpio', comm=comm) as f:
                # -- fetch state variables
                statevec_true[indx_map["Vx"],0] = f['Vx'][0]
                statevec_true[indx_map["Vy"],0] = f['Vy'][0]
                statevec_true[indx_map["Vz"],0] = f['Vz'][0]
                statevec_true[indx_map["Pressure"],0] = f['Pressure'][0]
        except Exception as e:
            print(f"[Generate-True-State: read init file] Error reading the file: {e}")
                            
        # -- mimic the time integration loop to save vec on every time step
        for k in range(kwargs.get('nt')):
            kwargs.update({'k': k})
            time = kwargs.get('t')
            kwargs.update({'tinitial': time[k], 'tfinal': time[k+1]})
            # --- write the state back to h5 file for ISSM model
            input_filename = f'{icesee_path}/{data_path}/ensemble_output_{rank}.h5'
            with h5py.File(input_filename, 'w', driver='mpio', comm=comm) as f:
                f.create_dataset('Vx', data=statevec_true[indx_map["Vx"],k])
                f.create_dataset('Vy', data=statevec_true[indx_map["Vy"],k])
                f.create_dataset('Vz', data=statevec_true[indx_map["Vz"],k])
                f.create_dataset('Pressure', data=statevec_true[indx_map["Pressure"],k])

            # -- call the run_model function to push the state forward in time
            ISSM_model(**kwargs)
           
            try:
                output_filename = f'{icesee_path}/{data_path}/ensemble_output_{ens_id}.h5'
                with h5py.File(output_filename, 'r', driver='mpio', comm=comm) as f:
                    statevec_true[indx_map["Vx"],k+1] = f['Vx'][0]
                    statevec_true[indx_map["Vy"],k+1] = f['Vy'][0]
                    statevec_true[indx_map["Vz"],k+1] = f['Vz'][0]
                    statevec_true[indx_map["Pressure"],k+1] = f['Pressure'][0]
                    vx = f['Vx'][0]
            except Exception as e:
                print(f"[Generate-True-State: read output file] Error reading the file: {e}")
                # return None
        
        updated_state = {'Vx': statevec_true[indx_map["Vx"],:],
                        'Vy': statevec_true[indx_map["Vy"],:],
                        'Vz': statevec_true[indx_map["Vz"],:],
                        'Pressure': statevec_true[indx_map["Pressure"],:]}

        #  --- change directory back to the original directory ---
        os.chdir(icesee_path)
        
        return updated_state
    
    except Exception as e:
        print(f"[Generate true state] Error sending command: {e}")
        # Ensure directory is changed back even on error
        os.chdir(icesee_path)
        return None
    

def generate_nurged_state(**kwargs):
    """generate the nurged state of the model"""
    params = kwargs.get('params')
    time   = kwargs.get('t')
    server = kwargs.get('server')
    issm_examples_dir   = kwargs.get('issm_examples_dir')
    icesee_path         = kwargs.get('icesee_path')
    data_path           = kwargs.get('data_path')
    comm                = kwargs.get('comm')

    #  --- change directory to the issm directory ---
    os.chdir(issm_examples_dir)

    # get the rank of the current process
    rank = comm.Get_rank()

    # --- filename for data saving
    fname = 'nurged_state.mat'
    kwargs.update({'fname': fname})
    ens_id = kwargs.get('ens_id')

    try:
    # if True:
        # --- fetch treu state vector
        statevec_nurged = kwargs.get('statevec_nurged')

        # -- call the icesee_get_index function to get the index of the state vector
        vecs, indx_map, dim_per_proc = icesee_get_index(statevec_nurged, **kwargs)

        # -- fetch data from inital state
        try: 
            output_filename = f'{icesee_path}/{data_path}/ensemble_init_{rank}.h5'
            # print(f"[DEBUG] Attempting to open file: {output_filename}")
            if not os.path.exists(output_filename):
                print(f"[ERROR] File does not exist: {output_filename}")
                return None
            with h5py.File(output_filename, 'r', driver='mpio', comm=comm) as f:
                # -- fetch state variables
                statevec_nurged[indx_map["Vx"],0] = f['Vx'][0]
                statevec_nurged[indx_map["Vy"],0] = f['Vy'][0]
                statevec_nurged[indx_map["Vz"],0] = f['Vz'][0]
                statevec_nurged[indx_map["Pressure"],0] = f['Pressure'][0]
        except Exception as e:
            print(f"[DEBUG] Error reading the file: {e}")
                            
        # -- mimic the time integration loop to save vec on every time step
        for k in range(kwargs.get('nt')):
            kwargs.update({'k': k})
            time = kwargs.get('t')
            kwargs.update({'tinitial': time[k], 'tfinal': time[k+1]})

            # --- write the state back to h5 file for ISSM model
            input_filename = f'{icesee_path}/{data_path}/ensemble_output_{ens_id}.h5'
            with h5py.File(input_filename, 'w', driver='mpio', comm=comm) as f:
                f.create_dataset('Vx', data=statevec_nurged[indx_map["Vx"],k])
                f.create_dataset('Vy', data=statevec_nurged[indx_map["Vy"],k])
                f.create_dataset('Vz', data=statevec_nurged[indx_map["Vz"],k])
                f.create_dataset('Pressure', data=statevec_nurged[indx_map["Pressure"],k])

            # -- call the run_model function to push the state forward in time
            ISSM_model(**kwargs)

            try:
                output_filename = f'{icesee_path}/{data_path}/ensemble_output_{ens_id}.h5'
                with h5py.File(output_filename, 'r', driver='mpio', comm=comm) as f:
                    statevec_nurged[indx_map["Vx"],k+1] = f['Vx'][0]
                    statevec_nurged[indx_map["Vy"],k+1] = f['Vy'][0]
                    statevec_nurged[indx_map["Vz"],k+1] = f['Vz'][0]
                    statevec_nurged[indx_map["Pressure"],k+1] = f['Pressure'][0]
            except Exception as e:
                print(f"[DEBUG] Error reading the file: {e}")
                # return None
        
        updated_state = {'Vx': statevec_nurged[indx_map["Vx"],:],
                        'Vy': statevec_nurged[indx_map["Vy"],:],
                        'Vz': statevec_nurged[indx_map["Vz"],:],
                        'Pressure': statevec_nurged[indx_map["Pressure"],:]}

        #  --- change directory back to the original directory ---
        os.chdir(icesee_path)
        
        # return updated_state
        return statevec_nurged
    
    except Exception as e:
        print(f"[DEBUG] Error sending command: {e}")
        # Ensure directory is changed back even on error
        os.chdir(icesee_path)
        return None
        
#  --- initialize ensemble members ---
def initialize_ensemble(ens, **kwargs):
    """des: initialize the ensemble members
    Returns: ensemble: the ensemble members
    """
    import h5py
    import os, sys

    server              = kwargs.get('server')
    issm_examples_dir   = kwargs.get('issm_examples_dir')
    icesee_path         = kwargs.get('icesee_path')
    data_path           = kwargs.get('data_path')
    comm                = kwargs.get('comm')

    #  --- change directory to the issm directory ---
    os.chdir(issm_examples_dir)
    ens_id = kwargs.get('ens_id')

    # get the rank of the current process
    rank = comm.Get_rank()
   
    try:
        output_filename = f'{icesee_path}/{data_path}/ensemble_init_{ens_id}.h5'
        # with h5py.File(output_filename, 'r', driver='mpio', comm=comm) as f:
        with h5py.File(output_filename, 'r') as f:
            # --- fetch state variables
            Vx = f['Vx'][0]
            Vy = f['Vy'][0]
            Vz = f['Vz'][0]
            Pressure = f['Pressure'][0]

        updated_state = {'Vx': Vx,
                        'Vy': Vy,
                        'Vz': Vz,
                        'Pressure': Pressure}
        os.chdir(icesee_path)

        return updated_state
        
    except Exception as e:
        print(f"[Initialze_ensemble] Error sending command: {e}")
        server.shutdown()
        server.reset_terminal()
        sys.exit(1)
    
    # # --- change directory back to the original directory ---
    # os.chdir(icesee_path)