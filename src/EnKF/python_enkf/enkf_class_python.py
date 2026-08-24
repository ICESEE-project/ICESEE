# =============================================================================
# @author: Brian Kyanjo
# @date: 2024-09-24
# @description: Optimized and Robust Ensemble Kalman Filter (EnKF) class with
#               Parallelization.
#               - The class includes the analysis step of the EnKF with both
#                 Stochastic and Deterministic EnKF options.
# ==============================================================================

import numpy as np
from scipy.stats import multivariate_normal
from joblib import Parallel, delayed

class EnsembleKalmanFilter:
    def __init__(self, ObsFun, JObsFun, Cov_obs, Cov_model, taper, icesee_kwargs, ensemble_threshold=30, n_jobs=-1):
        """
        Initializes the Ensemble Kalman Filter (EnKF) with parallelization.

        Parameters:
        ObsFun: Callable - Observation function mapping state space to observation space.
        JObsFun: Callable - Jacobian of the observation function.
        Cov_obs: ndarray - Observation covariance matrix (m x m).
        Cov_model: ndarray - Model covariance matrix (n x n).
        taper: ndarray - Covariance taper matrix (n x n).
        icesee_kwargs: dict - Dictionary containing parameters like "m_obs" and others.
        ensemble_threshold: int - Threshold for choosing between Stochastic and Deterministic EnKF.
        n_jobs: int - Number of parallel jobs to run (-1 uses all available cores).
        """
        self.ObsFun = ObsFun
        self.JObsFun = JObsFun
        self.Cov_obs = Cov_obs
        self.Cov_model = Cov_model
        self.taper = taper
        self.icesee_kwargs = icesee_kwargs
        self.ensemble_threshold = ensemble_threshold
        self.n_jobs = n_jobs  # Number of parallel jobs for joblib
        self._validate_input()

    def _validate_input(self):
        """Validate input matrices for correct dimensions and consistency."""
        if not all(isinstance(matrix, np.ndarray) for matrix in [self.Cov_obs, self.Cov_model, self.taper]):
            raise ValueError("All matrices must be numpy arrays.")
        if self.Cov_obs.shape[0] != self.Cov_obs.shape[1]:
            raise ValueError("Observation covariance matrix must be square.")
        if self.Cov_model.shape[0] != self.Cov_model.shape[1]:
            raise ValueError("Model covariance matrix must be square.")
        if self.taper.shape != self.Cov_model.shape:
            raise ValueError("Taper matrix must have the same dimensions as the model covariance matrix.")

    def _compute_kalman_gain(self):
        """Compute the Kalman gain based on the Jacobian of the observation function."""
        Jobs = self.JObsFun(self.Cov_model.shape[0], self.icesee_kwargs["m_obs"])
        inv_matrix = np.linalg.inv(Jobs @ self.Cov_model @ Jobs.T + self.Cov_obs)
        KalGain = self.Cov_model @ Jobs.T @ inv_matrix
        return KalGain

    def _generate_virtual_observations(self, huxg_obs, num_samples):
        """Generate virtual observations based on a multivariate normal distribution."""
        return multivariate_normal.rvs(mean=huxg_obs, cov=self.Cov_obs, size=num_samples).T

    def _update_ensemble_member(self, member_idx, huxg_ens, obs_virtual, KalGain):
        """Update a single ensemble member using the stochastic analysis step."""
        return huxg_ens[:, member_idx] + KalGain @ (obs_virtual[:, member_idx] - self.ObsFun(huxg_ens[:, member_idx], self.icesee_kwargs["m_obs"]))

    def EnKF_analyze(self, huxg_ens, huxg_obs):
        """
        Perform the analysis step of the Ensemble Kalman Filter (Stochastic EnKF) with parallelization.

        Parameters:
        huxg_ens: ndarray (n x N) - Ensemble matrix of model states (n: state size, N: ensemble size).
        huxg_obs: ndarray (m,) - Observation vector (m: measurement size).

        Returns:
        analysis_ens: ndarray (n x N) - The updated ensemble after analysis.
        analysis_cov: ndarray (n x n) - The updated covariance after analysis.
        """
        self._check_input_dimensions(huxg_ens, huxg_obs)

        n, N = huxg_ens.shape
        KalGain = self._compute_kalman_gain()

        # Generate virtual observations for all ensemble members
        obs_virtual = self._generate_virtual_observations(huxg_obs, N)

        # Use parallel computing to update each ensemble member concurrently
        analysis_ens = np.array(Parallel(n_jobs=self.n_jobs)(
            delayed(self._update_ensemble_member)(i, huxg_ens, obs_virtual, KalGain) for i in range(N)
        )).T

        # Calculate analysis covariance matrix
        analysis_cov = self._compute_analysis_covariance(analysis_ens)

        return analysis_ens, analysis_cov

    def DEnKF_analyze(self, huxg_ens, huxg_obs):
        """
        Perform the analysis step of the Deterministic Ensemble Kalman Filter (DEnKF) with parallelization.

        Parameters:
        huxg_ens: ndarray (n x N) - Ensemble matrix of model states (n: state size, N: ensemble size).
        huxg_obs: ndarray (m,) - Observation vector (m: measurement size).

        Returns:
        analysis_ens: ndarray (n x N) - The updated ensemble after analysis.
        analysis_cov: ndarray (n x n) - The updated covariance after analysis.
        """
        self._check_input_dimensions(huxg_ens, huxg_obs)

        n, N = huxg_ens.shape
        KalGain = self._compute_kalman_gain()

        # Compute analysis of ensemble mean
        huxg_ens_mean = np.mean(huxg_ens, axis=1, keepdims=True)
        analysis_mean = huxg_ens_mean + KalGain @ (huxg_obs - self.ObsFun(huxg_ens_mean, self.icesee_kwargs["m_obs"]))

        # Create parallel computation for anomalies
        huxg_ens_anom = huxg_ens - huxg_ens_mean
        analysis_ens_anom = np.array(Parallel(n_jobs=self.n_jobs)(
            delayed(lambda i: huxg_ens_anom[:, i] - 0.5 * KalGain @ self.ObsFun(huxg_ens_anom[:, i], self.icesee_kwargs["m_obs"]))
            for i in range(N)
        )).T

        # Compute analysis ensemble and covariance matrix
        analysis_ens = analysis_ens_anom + analysis_mean
        analysis_cov = self._compute_analysis_covariance(analysis_ens)

        return analysis_ens, analysis_cov

    def analyze(self, huxg_ens, huxg_obs):
        """Automatically select between Stochastic and Deterministic EnKF based on ensemble size."""
        if huxg_ens.shape[1] < self.ensemble_threshold:
            print(f"Using Deterministic EnKF (DEnKF) with ensemble size {huxg_ens.shape[1]}.")
            return self.DEnKF_analyze(huxg_ens, huxg_obs)
        else:
            print(f"Using Stochastic EnKF with ensemble size {huxg_ens.shape[1]}.")
            return self.EnKF_analyze(huxg_ens, huxg_obs)

    def _compute_analysis_covariance(self, analysis_ens):
        """Compute the covariance of the analysis ensemble."""
        N = analysis_ens.shape[1]
        analysis_ens_mean = np.mean(analysis_ens, axis=1, keepdims=True)
        analysis_cov = (1 / (N - 1)) * (analysis_ens - analysis_ens_mean) @ (analysis_ens - analysis_ens_mean).T
        return analysis_cov * self.taper

    def _check_input_dimensions(self, huxg_ens, huxg_obs):
        """Check dimensions of ensemble and observation vectors."""
        if huxg_ens.ndim != 2 or huxg_obs.ndim != 1:
            raise ValueError("Ensemble matrix must be 2D and observation vector must be 1D.")

    def set_cov_model(self, Cov_model):
        """Update the model covariance matrix."""
        self._check_square_matrix(Cov_model)
        if Cov_model.shape != self.Cov_model.shape:
            raise ValueError("New model covariance matrix must have the same dimensions as the current one.")
        self.Cov_model = Cov_model

    def set_cov_obs(self, Cov_obs):
        """Update the observation covariance matrix."""
        self._check_square_matrix(Cov_obs)
        if Cov_obs.shape != self.Cov_obs.shape:
            raise ValueError("New observation covariance matrix must have the same dimensions as the current one.")
        self.Cov_obs = Cov_obs

    def _check_square_matrix(self, matrix):
        """Ensure the matrix is square."""
        if matrix.shape[0] != matrix.shape[1]:
            raise ValueError("Matrix must be square.")
