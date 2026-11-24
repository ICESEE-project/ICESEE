"""Independent Kalman references for the parallel stochastic EnKF gain."""
import numpy as np
import pytest

from ICESEE.src.parallelization.distributed_analysis import (
    StochasticAnalysisProducts, ensemble_transform_from_products,
    prescribed_r_transform, resolve_gain_formulation,
)


@pytest.mark.parametrize("observations", [8, 32, 128, 1024])
@pytest.mark.parametrize("prior_sd", [1., .01, .001])
def test_observation_rich_scalar_gain_matches_known_R(observations, prior_sd):
    n = 32
    x = np.linspace(-1., 1., n)
    x = prior_sd * (x - x.mean()) / x.std(ddof=1)
    y = np.tile(x, (observations, 1))
    rng = np.random.default_rng(8)
    errors = .1 * rng.standard_normal(y.shape)
    errors -= errors.mean(axis=1, keepdims=True)
    products = StochasticAnalysisProducts.zeros(n, gain_formulation="prescribed_r")
    for start in range(0, observations, 17):
        rows = slice(start, start + 17)
        products.add_chunk(y[rows], np.ones(observations)[rows],
            error_mode="stochastic_R", observation_errors=errors[rows],
            observation_variances=np.full(observations, .01)[rows])
    transform = ensemble_transform_from_products(products)
    expected_gain = observations * prior_sd**2 / (.01 + observations * prior_sd**2)
    np.testing.assert_allclose((x @ transform).mean(), expected_gain, atol=2e-11, rtol=2e-11)
    np.testing.assert_allclose(transform.sum(axis=0), 1., atol=2e-12)


def fixture():
    rng = np.random.default_rng(19)
    x = rng.normal(size=(5, 32))
    h = rng.normal(size=(11, 5))
    d = rng.normal(size=11)
    sigma = np.linspace(.3, 1.2, 11)
    eta = rng.normal(size=(11, 32)) * sigma[:, None]
    eta -= eta.mean(axis=1, keepdims=True)
    y = h @ x
    p = np.cov(x)
    k = np.linalg.solve(h @ p @ h.T + np.diag(sigma**2), h @ p).T
    return x, h, d, sigma, eta, y, p, k


def test_memberwise_dense_Kalman_reference_and_permutations():
    x, h, d, sigma, eta, y, p, k = fixture()
    yp = y - y.mean(axis=1, keepdims=True)
    innovation = d[:, None] + eta - y
    transform = prescribed_r_transform(yp, innovation, sigma**2)
    np.testing.assert_allclose(x @ transform, x + k @ innovation, rtol=1e-12, atol=1e-12)
    rng = np.random.default_rng(4)
    members, rows = rng.permutation(32), rng.permutation(11)
    other = prescribed_r_transform(yp[rows][:, members], innovation[rows][:, members], sigma[rows]**2)
    np.testing.assert_allclose(x[:, members] @ other, (x @ transform)[:, members], atol=2e-12)


def test_stochastic_posterior_covariance_in_expectation():
    x, h, d, sigma, _, y, p, k = fixture()
    expected = p - k @ h @ p
    samples = []
    for seed in range(600):
        eta = np.random.default_rng(seed).normal(size=y.shape) * sigma[:, None]
        eta -= eta.mean(axis=1, keepdims=True)
        result = x @ prescribed_r_transform(y-y.mean(axis=1,keepdims=True), d[:,None]+eta-y, sigma**2)
        samples.append(np.cov(result))
    # Monte Carlo check of covariance, not a member-wise deterministic-filter assertion.
    relative = np.linalg.norm(np.mean(samples,axis=0)-expected) / np.linalg.norm(expected)
    assert relative < .035


@pytest.mark.parametrize("variance", [np.zeros(4), np.full(4,np.nan), np.ones((4,4)), np.ones(3)])
def test_invalid_covariance_is_rejected(variance):
    with pytest.raises(ValueError):
        prescribed_r_transform(np.ones((4,3)), np.ones((4,3)), variance)


def test_prescribed_R_retains_null_uncertainty_and_weak_observation_limit():
    np.testing.assert_array_equal(prescribed_r_transform(np.zeros((4,3)),np.ones((4,3)),np.ones(4)),np.eye(3))
    t = prescribed_r_transform(np.array([[-1.,0.,1.]]),np.array([[2.,1.,0.]]),np.array([1e20]))
    np.testing.assert_allclose(t,np.eye(3),atol=1e-18)
    assert resolve_gain_formulation({}) == "sampled_factor"
    with pytest.raises(ValueError):
        resolve_gain_formulation(dict(enkf_gain_formulation="prescribed_r"))
    with pytest.raises(NotImplementedError):
        resolve_gain_formulation(dict(enkf_gain_formulation="prescribed_r",enkf_observation_error_mode="stochastic_R",local_analysis=True))


def test_mode1_and_file_backed_mode2_use_same_prescribed_covariance(tmp_path):
    import h5py
    from mpi4py import MPI
    from ICESEE.src.parallelization._mpi_analysis_functions import EnKF_X5
    from ICESEE.src.parallelization.EnKF_parallel_io import EnKF_fully_parallel_IO
    rng=np.random.default_rng(6);x=rng.normal(size=(11,32));rows=np.array([9,2,0,8,5,7,4])
    sigma=np.linspace(.2,1.,7);d=rng.normal(size=7)
    options=dict(enkf_gain_formulation='prescribed_r',enkf_observation_error_mode='stochastic_R',
        Nens=32,nd=11,base_seed=42,k=0,km=0,local_analysis=False,localization_flag=False,
        time_analysis_file_writing=0.)
    full=np.zeros((11,1));full[rows,0]=d
    std=np.zeros(11);std[rows]=sigma
    class Utils:
        def __init__(self,**kw):pass
        def JObs_indices(self,nd):return rows
    mode1,_=EnKF_X5(0,x,32,full,dict(options,error_R=std),Utils)
    with h5py.File(tmp_path/'H_matrix.h5','w') as f:f['obs_indices']=rows
    with h5py.File(tmp_path/'synthetic_obs.h5','w') as f:
        f['hu_obs_compact']=d[:,None];f['obs_std']=sigma
    io=EnKF_fully_parallel_IO.__new__(EnKF_fully_parallel_IO)
    io.nens=32;io.nd=11;io.nt=3;io.rank=0;io.base_path=str(tmp_path);io.mpi_comm=MPI.COMM_SELF
    mode2,_,_=io.compute_X5_modified(state_dset=x,**options)
    np.testing.assert_allclose(x@mode1,x@mode2,atol=2e-12)
    eta=np.random.default_rng(42+1000003).standard_normal((7,32))*sigma[:,None]
    eta-=eta.mean(axis=1,keepdims=True)
    p=np.cov(x);gain=np.linalg.solve(p[np.ix_(rows,rows)]+np.diag(sigma**2),p[rows,:]).T
    np.testing.assert_allclose(x@mode2,x+gain@(d[:,None]+eta-x[rows]),atol=2e-12)
