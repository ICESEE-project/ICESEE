import numpy as np
import pytest
from scipy.stats import multivariate_normal
from ICESEE.src.EnKF.python_enkf.EnKF import EnsembleKalmanFilter
from ICESEE.src.utils.localization import active_observation_std


def test_serial_covariance_uses_observed_rows_and_current_time():
    # Three observation times, six active rows, heterogeneous error scales.
    error = np.array([[1.,2.,3.,4.,5.,6.,7.,8.],
                      [2.,4.,6.,8.,10.,12.,14.,16.],
                      [3.,6.,9.,12.,15.,18.,21.,24.]])
    rows = np.array([0,1,2,4,6,7])
    covariance = np.diag(active_observation_std({"error_R":error},1,rows)**2)
    assert covariance.shape == (6,6)
    np.testing.assert_equal(np.diag(covariance),[4.,16.,36.,100.,196.,256.])


@pytest.mark.parametrize("n,m",[(1,1),(3,1),(1,3),(4,3)])
def test_batched_observation_draws_preserve_member_rng_sequence(monkeypatch,n,m):
    cov = np.diag(np.arange(1,m+1,dtype=float))
    observed = np.arange(m,dtype=float)
    members = np.zeros((m,n))
    analysis = EnsembleKalmanFilter(Observation_vec=observed,Cov_obs=cov,
                                   Observation_function=lambda x:x)
    monkeypatch.setattr(analysis,"_compute_kalman_gain",lambda:np.eye(m))
    np.random.seed(481)
    expected = np.column_stack([observed+np.asarray(
        multivariate_normal.rvs(mean=np.zeros(m),cov=cov)).reshape(m)
        for _ in range(n)])
    expected_state = np.random.get_state()
    np.random.seed(481)
    actual = analysis.EnKF_Analysis(members)
    actual_state = np.random.get_state()
    np.testing.assert_array_equal(actual,expected)
    np.testing.assert_array_equal(actual_state[1],expected_state[1])
    assert actual_state[2:] == expected_state[2:]
