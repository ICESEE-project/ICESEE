# =============================================================================
# @author: Brian Kyanjo
# @date: 2024-11-06
# @description: Conatins the class for analysis steps all Ensemble Kalman Filter (EnKF)
#              EnKF include both Stochastic and Deterministic EnKF options: EnKF,
#              EnSRF, DEnKF, and EnTKF.
#             - The class further includes parallelized forecast and analysis steps with
#               MPI, Dask, Ray, and Multiprocessing.
# =============================================================================

import os, sys
import re
import numpy as np
import scipy
from scipy.stats import multivariate_normal

# from _utility_imports import *
# src_dir = os.path.join(project_root, 'src', 'parallelization')
# sys.path.insert(0, src_dir)
from ICESEE.src.parallelization.parallel_mpi.icesee_mpi_parallel_manager import ParallelManager
from ICESEE.config._utility_imports import icesee_get_index
from ICESEE.src.run_model_da._error_generation import generate_enkf_field


# Move `worker` to global scope
def worker(args):
    """Wrapper function for multiprocessing."""
    ens_idx, ensemble_member, nd, Q_err, params, model_kwargs= args
    return EnsembleKalmanFilter.forecast_step_single(ensemble_member[:, ens_idx], **model_kwargs)

class EnsembleKalmanFilter:
    def __init__(self, Observation_vec=None, Cov_obs=None, Cov_model=None, \
                 Observation_function=None, Obs_Jacobian=None, parameters=None, taper_matrix=None, parallel_manager = None, parallel_flag="serial"):
        """
        Initializes the Analysis class for the Ensemble Kalman Filter (EnKF).
        
        Parameters:
        Observation_function: Callable - Observation function mapping state space to observation space.
        observation_vec: ndarray - Observation vector (m x 1).
        Obs_Jacobian: Callable - Jacobian of the observation function.
        Cov_obs: ndarray - Observation covariance matrix (m x m).
        Cov_model: ndarray - Model covariance matrix (n x n).
        taper: ndarray - Covariance taper matrix (n x n).
        params: dict - Dictionary containing parameters like "m_obs" and others.\
        parallel_flag: str - Flag for parallelization (serial,MPI, Dask, Ray, Multiprocessing).
        """
        self.Observation_vec        = Observation_vec
        self.Cov_obs                = Cov_obs
        self.Cov_model              = Cov_model
        self.Obs_Jacobian           = Obs_Jacobian
        self.parameters             = parameters
        self.taper_matrix           = taper_matrix
        self.Observation_function   = Observation_function
        self.parallel_flag          = parallel_flag
        self.parallel_manager       = parallel_manager
    
    # Forecast step
    def forecast_step(self, ensemble=None, forecast_step_single=None, **model_kwargs):
        """
        Forecast step for the Ensemble Kalman Filter (EnKF).
        
        Parameters:
            forecast_step_single: Callable - Function for the forecast step of each ensemble member.
            Q_err: ndarray - Process noise matrix.
            parallel_flag: str - Flag for parallelization (serial,MPI, Dask, Ray, Multiprocessing).
            **model_kwargs: dict - Keyword arguments for the model.
        
        Returns:
            ensemble: ndarray - Updated ensemble matrix.
        """
        Q_err = model_kwargs.get("Q_err")
        L_C           = model_kwargs.get("L_C", None)
        Lx             = model_kwargs.get("Lx", 1.0)
        Ly             = model_kwargs.get("Ly", 1.0)
        len_scale      = model_kwargs.get("len_scale", 1.0)
        Q_rho          = model_kwargs.get("Q_rho", 1.0)
        params         = model_kwargs.get("params", {})
        alpha         = model_kwargs.get("alpha", 0.0)
        dt             = model_kwargs.get("dt", 1.0)
        rho           = model_kwargs.get("rho", 1.0)
        noise         = model_kwargs.get("noise", None)

        
        if re.match(r"\Aserial\Z", self.parallel_flag, re.IGNORECASE):
            # Serial forecast step
            nd, Nens = ensemble.shape # Get the number of ensemble members
            if model_kwargs["joint_estimation"] or params["localization_flag"]:
                hdim = ensemble.shape[0] // params["total_state_param_vars"]
            else:
                hdim = ensemble.shape[0] // params["num_state_vars"]
            state_block_size = hdim * params["num_state_vars"]
            vecs, indx_map, dim_per_proc = icesee_get_index(**model_kwargs)
            # Loop over the ensemble members
            for ens in range(Nens):
                updated_state = forecast_step_single(ensemble=ensemble[:, ens], **model_kwargs)
                for key,value in updated_state.items():
                    ensemble[indx_map[key],ens] = value

                # add process noise
                noise_all = []
                q0 = []
                for ii, sig in enumerate(params["sig_Q"]):
                    if ii <=params["num_state_vars"]:
                        for jj, key in enumerate(model_kwargs['vec_inputs']):
                            if ii == jj:
                                model_kwargs.update({"ii_sig": ii, "Lx_dim": np.sqrt(Lx*Ly), "noise_dim": len(indx_map[key]), "hdim":hdim, "num_vars":params["total_state_param_vars"]})
                                W = generate_enkf_field(**model_kwargs)
                                noise_ = alpha*noise[indx_map[key]] + np.sqrt(1 - alpha**2)*W
                                q0.append(noise_)

                                Z = np.sqrt(dt)*sig*rho*noise_
                                noise_all.append(Z)
                            
                noise_ = np.concatenate(noise_all, axis=0)

                # only update the state variables with noise
                for ii, key in enumerate(model_kwargs['vec_inputs']):
                    if ii < params["num_state_vars"]:
                        ensemble[indx_map[key],ens] += noise_[indx_map[key]]

                # ensemble[:state_block_size,ens] += noise_[:state_block_size]
                noise = np.concatenate(q0, axis=0)
                model_kwargs.update({"noise": noise})
                del noise_all, q0, noise_, W

            return ensemble

        # Using divide and conquer parallelization with MPI for non-MPI application in forecast_step_single
        elif re.match(r"\ANon-mpi-model\Z", self.parallel_flag, re.IGNORECASE):
            """
            Called only when the numerical model in the forecast step is not MPI parallelized
            """
            from mpi4py import MPI

            comm = MPI.COMM_WORLD
            rank = comm.Get_rank()

            # Get the number of ensemble members
            nd, Nens = ensemble.shape

            size = min(comm.Get_size(), Nens)  # Limit processes to ensemble size

            # Balanced workload
            chunk_sizes = [(Nens // size) + (1 if i < (Nens % size) else 0) for i in range(size)]
            displacements = [sum(chunk_sizes[:i]) for i in range(size)]
            start_idx = displacements[rank]
            end_idx = start_idx + chunk_sizes[rank]

            # Local chunk
            local_ensemble = ensemble[:, start_idx:end_idx]

            # Perform forecast step
            for ens in range(local_ensemble.shape[1]):
                local_ensemble[:, ens] = forecast_step_single(ensemble=local_ensemble[:,ens], **model_kwargs)

            # Avoid gather; update ensemble in place
            gathered_ensemble = comm.allgather(local_ensemble)

            if rank == 0:
                # print(f"Gathered ensemble shape: {np.hstack(gathered_ensemble).shape}")
                ensemble = np.hstack(gathered_ensemble)
            else:
                ensemble = None

            return ensemble
        
        # Using MPI split communicator for MPI application in forecast_step_single
        elif re.match(r"\AMPI_model\Z", self.parallel_flag, re.IGNORECASE):
            """
            if the numerical model ran in the forecast step is MPI parallelized
            - use parallelization with MPI split communicator (COMM_model) from the icesee_mpi_parallelization
            routine
            """
            nd = ensemble.shape[0]

            size = self.parallel_manager.size_model
            rank = self.parallel_manager.rank_model
            comm = self.parallel_manager.COMM_model
            # check the number of available processors
            print(f"\nranks: {rank}, size: {size}\n")

            for ens in range(ensemble.shape[1]):
                ensemble[:, ens] = forecast_step_single(ensemble=ensemble[:, ens], **model_kwargs)

            gathered_ensemble = self.parallel_manager.all_gather_data(comm, ensemble)
            #  initalize the ensemble to be returned
            ensemble_bcast = np.empty((nd,self.parameters["Nens"]))
            comm.Barrier()
            if rank == 0:
                ensemble_bcast = np.hstack(gathered_ensemble)
            # else:
                # ensemble = None # return None for other ranks since 
            self.parallel_manager.broadcast_data(comm, ensemble_bcast, root=0)
            # ensemble = comm.bcast(ensemble, root=0)
            return ensemble_bcast

            # gather the ensemble members
            # comm = self.parallel_manager.COMM_model
            # comm.Barrier()
            # gathered_ensemble = self.parallel_manager.all_gather_data(comm, ensemble)
            # gathered_ensemble = comm.allgather(ensemble)
            # print(f"[Rank {self.parallel_manager.rank_model}] Gathered shapes: {[arr.shape for arr in gathered_ensemble]}")

            # debug print
            # print(f"[Rank {self.parallel_manager.rank_model}] Ensemble shape after forecast: {np.hstack(gathered_ensemble).shape}")
            # print(f"[Rank {self.parallel_manager.rank_model}] Ensemble mean after forecast: {np.mean(np.hstack(gathered_ensemble), axis=1)}")

            # comm.Barrier()
            # if self.parallel_manager.rank_model == 0:
            #     ensemble = np.hstack(gathered_ensemble)
            #     # print(f"Ensemble shape after gathering: {ensemble.shape}")
            #     return ensemble
            # else:
                # return None
            

        # Parallel forecast step using Multiprocessing
        elif re.match(r"\AMultiprocessing\Z", self.parallel_flag, re.IGNORECASE):
            """
            Run the forecast_step_single function Nens times simultaneously using multiprocessing.
            """

            import multiprocessing as mp
            from multiprocessing import Pool

            nd, Nens = ensemble.shape
            print("before",ensemble[:10,1])
            # print(f"[DEBUG] EnKF Mean Before Update: {np.mean(ensemble, axis=1)}")

            # Prepare worker arguments (avoid copying full ensemble, just one member per worker)
            worker_args = [(ens_idx, ensemble, nd, Q_err, self.parameters, model_kwargs) for ens_idx in range(Nens)]

            # Run forecast_step_single in parallel
            # if __name__ == '__main__':
                # restrict the number of workers to the number of ensemble members

            
            nworkers = min(mp.cpu_count(), Nens)

            with Pool(nworkers) as pool:
                results = pool.map(worker, worker_args)

                # Convert the list of results back to a numpy array and transpose
                # ensemble = np.array(results).T
                # Convert results to numpy array
                results_array = np.array(results)

                # Debugging print for multiprocessing results
                # print(f"[DEBUG] Raw multiprocessing results shape: {results_array.shape}")

                # If results is 1D (should be 2D), reshape it
                # if results_array.ndim == 1:
                #     results_array = results_array.reshape(nd, Nens)

                # If results array is empty, raise an error
                if results_array.size == 0:
                    raise ValueError("Multiprocessing returned an empty ensemble!")

                # Transpose for expected shape
                ensemble = results_array.T  
            # print("After", ensemble[:10,1])
            # print(f"[DEBUG] EnKF Mean After Update: {np.mean(ensemble, axis=1)}")
            return ensemble


        # Parallel forecast step using Dask
        elif re.match(r"\ADask\Z", self.parallel_flag, re.IGNORECASE):
            import dask
            import dask.array as da
            from dask import compute, delayed
            import copy

            nd, Nens = ensemble.shape

            # Create delayed tasks for each ensemble member
            tasks = [
                    delayed(forecast_step_single)(ensemble[:, ens], **model_kwargs)
                    for ens in range(Nens)
                ]

            # Compute all tasks in parallel and collect the results
            results = compute(*tasks)

            # Update the ensemble matrix
            ensemble = np.array(results).T

            return ensemble

        # Parallel forecast step using Ray
        elif re.match(r"\ARay\Z", self.parallel_flag, re.IGNORECASE):
            import ray

            nd, Nens = ensemble.shape

            # Initialize Ray (Only needed once, preferably at the beginning of the script)
            if not ray.is_initialized():
                ray.init(ignore_reinit_error=True)

            @ray.remote
            def ray_worker(ensemble_member, nd, Q_err, parameters, model_kwargs):
                """
                Remote function to perform forecast step for a single ensemble member.
                This function will be executed in parallel by Ray workers.
                """
                return forecast_step_single(ensemble_member, **model_kwargs)

            # Launch tasks in parallel using Ray
            futures = [
                ray_worker.remote(ensemble[:, ens].copy())
                for ens in range(Nens)
            ]

            # Collect results from all Ray tasks
            results = ray.get(futures)

            # Convert the list of results back to a numpy array and transpose
            ensemble = np.array(results).T
            return ensemble

        
        # python openmp parallelization
        elif re.match(r"\APyomp\Z", self.parallel_flag, re.IGNORECASE):
            # check if the Pyomp module is installed
            try:
                # check python version: it should be [3.9 - 3.10]
                if sys.version_info[0] == 3 and (sys.version_info[1] > 9 or sys.version_info[1] <= 10):
                    from numba import njit
                    from numba.openmp import openmp_context as openmp
                    from numba.openmp import omp_get_thread_num, omp_get_num_threads
                    from collections import namedtuple

                    nd, Nens = ensemble.shape

                    #  get the maximum number of cores
                    MaxTHREADS = os.cpu_count()

                    ModelKwargs = namedtuple("ModelKwargs", model_kwargs.keys())
                    processed_args = ModelKwargs(**model_kwargs)

                    # create a wrapper for the forecast step function to use @njit
                    @njit
                    def forecast_wrapper():
                        partialEnsembles = np.zeros((nd,MaxTHREADS))
                        with openmp("parallel shared(partialEnsembles, numThrds) private(threadID,ens,localEnsemble)"):
                            threadID = omp_get_thread_num() #thread rank = 0...(numThreads-1)
                            with openmp("single"): #one thread does the work, others wait
                                numThrds = omp_get_num_threads() #get number of threads
                            for ens in range(threadID,Nens,numThrds):
                                partialEnsembles[:,threadID] = forecast_step_single(ensemble=ensemble[:,ens], **processed_args)
                                
                        return partialEnsembles

                    return forecast_wrapper()
                
            except ImportError:
                raise ImportError("Pyomp module not found. Please install it using 'conda install -c python-for-hpc -c conda-forge pyomp' of see https://github.com/Python-for-HPC/PyOMP for more details.")
            
        else:
            raise ValueError("Invalid parallel flag. Choose from 'serial', 'MPI', 'Dask', 'Ray', 'Multiprocessing'.")

    # Kalman gain computation
    def _compute_kalman_gain(self):
        """Compute the Kalman gain based on the Jacobian of the observation function."""
        # compute the mean of the forecast ensemble
        # ensemble_forecast_mean = np.mean(ensemble, axis=1)
    
        Jobs = self.Obs_Jacobian(self.Cov_model.shape[0])
        inv_matrix = np.linalg.inv(Jobs @ self.Cov_model @ Jobs.T + self.Cov_obs)
        KalGain = self.Cov_model @ Jobs.T @ inv_matrix
        return KalGain
    
    # Analysis steps
    def EnKF_Analysis(self, ensemble):
        """
        Stochastic Ensemble Kalman Filter (EnKF) analysis step.
        
        Parameters:
            ensemble: ndarray - Ensemble matrix (n x N).
        
        Returns:
            ensemble_analysis: updated ensemble matrix (n x N).
            analysis_error_cov: ndarray - Analysis error covariance matrix (n x n).
        """
        # Compute the Kalman gain
        KalGain = self._compute_kalman_gain()

        n,Nens = ensemble.shape
        m   =self.Observation_vec.shape[0] 

        # compute virtual observations and ensemble analysis
        virtual_observations  = np.zeros((m,Nens))
        ensemble_analysis     = np.zeros_like(ensemble)

        for ens in range(Nens):
            # Generate virtual observations
            virtual_observations[:,ens] =  self.Observation_vec+ multivariate_normal.rvs(mean=np.zeros(m), cov=self.Cov_obs)
            # Compute the analysis step
            ensemble_analysis[:,ens] = ensemble[:,ens] + KalGain @ (virtual_observations[:,ens] - self.Observation_function(ensemble[:,ens]))

        # compute the analysis error covariance
        # difference = ensemble_analysis - np.mean(ensemble_analysis, axis=1,keepdims=True)
        # analysis_error_cov =(1/(Nens-1)) * difference @ difference.T

        # Debug
        # print(f"[Debug] max of analysis error covariance: {np.max(analysis_error_cov[567:,567:])}")

        # return ensemble_analysis, analysis_error_cov
        return ensemble_analysis
       
    def DEnKF_Analysis(self, ensemble):
        """
        Deterministic Ensemble Kalman Filter (DEnKF) analysis step.
        
        Parameters:
            ensemble: ndarray - Ensemble matrix (n x N).
        
        Returns:
            ensemble_analysis: updated ensemble matrix (n x N).
            analysis_error_cov: ndarray - Analysis error covariance matrix (n x n).
        """
        # Compute the Kalman gain
        KalGain = self._compute_kalman_gain()

        n,N = ensemble.shape
        m   = self.Observation_vec.shape[0] 

        # compute ensemble mean
        ensemble_forecast_mean = np.mean(ensemble, axis=1)

        # compute the anlysis mean
        analysis_mean = ensemble_forecast_mean + KalGain @ (self.Observation_vec - self.Observation_function(ensemble_forecast_mean))

        # compute the forecast and analysis anomalies
        forecast_anomalies = np.zeros_like(ensemble)
        analysis_anomalies = np.zeros_like(ensemble)

        for i in range(N):
            forecast_anomalies[:,i] = ensemble[:,i] - ensemble_forecast_mean
            analysis_anomalies[:,i] = forecast_anomalies[:,i] - 0.5 * KalGain @ self.Observation_function(forecast_anomalies[:,i])

        # compute the analysis ensemble
        ensemble_analysis = analysis_anomalies + analysis_mean.reshape(-1,1)

        # compute the analysis error covariance
        # analysis_error_cov =(1/(N-1)) * analysis_anomalies @ analysis_anomalies.T

        # return ensemble_analysis, analysis_error_cov
        return ensemble_analysis

    def EnRSKF_Analysis(self, ensemble):
        """
        Deterministic Ensemble Square Root Filter (EnSRF) analysis step.
        
        Parameters:
            ensemble: ndarray - Ensemble matrix (n x N).
        
        Returns:
            ensemble_analysis: updated ensemble matrix (n x N).
            analysis_error_cov: ndarray - Analysis error covariance matrix (n x n).
        """
        _,N = ensemble.shape
        I = np.eye(N)

        # compute the mean of the forecast ensemble
        ensemble_forecast_mean = np.mean(ensemble, axis=1)

        obs_anomaly =self.Observation_vec.reshape(-1,1) - self.Observation_function(ensemble)

        V =  self.Observation_function(self.Cov_model)
        IN = self.Cov_obs + V@V.T

        # compute the analysis ensemble
        ensemble_analysis  = ensemble_forecast_mean.reshape(-1,1) + self.Cov_model @ V.T @ np.linalg.solve(IN, obs_anomaly)

        # #  singular value decomposition
        # U,S,Vt = np.linalg.svd(I - V.T @ np.linalg.solve(IN, V))

        # # compute the analysis error covariance
        # analysis_error_cov = ensemble_analysis + self.Cov_model@(U@np.diag(np.sqrt(S))@U.T)

        return ensemble_analysis
    
    def EnTKF_Analysis(self, ensemble):
        """
        Ensemble Transform Kalman Filter (EnTKF) analysis step.
        
        Parameters:
            ensemble: ndarray - Ensemble matrix (n x N).
        
        Returns:
            ensemble_analysis: updated ensemble matrix (n x N).
            analysis_error_cov: ndarray - Analysis error covariance matrix (n x n).
        """
        
        # compute the mean of the forecast ensemble
        ensemble_forecast_mean = np.mean(ensemble, axis=1)

        obs_anomaly =self.Observation_vec.reshape(-1,1) - self.Observation_function(ensemble)

        _obs_operator = self.Observation_function(self.Cov_model)

        # perform singular value decomposition
        U,S, _ = np.linalg.svd(_obs_operator.T @ np.linalg.solve(self.Cov_obs, _obs_operator))

        # compute the inverse of the analysis error covariance
        inv_analysis_error_cov = U @ np.diag(1/(S+1)) @ U.T

        # Right hand side of the analysis equation
        rhs = _obs_operator.T @ np.linalg.solve(self.Cov_obs, obs_anomaly)

        # compute the analysis ensemble
        ensemble_analysis = ensemble_forecast_mean.reshape(-1,1) + self.Cov_model@inv_analysis_error_cov@rhs

        # compute the analysis ensemble increment
        # ensemble_increment = self.Cov_model @ U @ np.diag(np.sqrt(1/(S+1))) @ U.T

        # # compute the analysis error covariance
        # analysis_error_cov = ensemble_analysis + ensemble_increment

        # return ensemble_analysis, analysis_error_cov
        return ensemble_analysis
    
    
# if __name__ == '__main__':
#     enkf = EnsembleKalmanFilter()
#     enkf.forecast_step()
