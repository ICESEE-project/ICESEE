import h5py
import numpy as np

from ICESEE.src.parallelization.EnKF_parallel_io import (
    EnKF_fully_parallel_IO,
    _exclude_inversion_observations,
    _inversion_variable_info,
    _read_rows_preserve_order,
    ensemble_transform_from_products,
    hdf5_compression_options,
    legacy_transform_from_gram,
    normalize_hdf5_compression,
    resolve_row_chunk_size,
)
from ICESEE.src.parallelization._mpi_analysis_functions import (
    apply_analysis_controls_local,
)
from ICESEE.src.parallelization import _mpi_forecast_functions as forecast_module
from ICESEE.src.parallelization._mpi_forecast_functions import (
    add_member_process_noise,
    forecast_member_context,
    process_noise_is_due,
)
from ICESEE.src.parallelization import _mpi_ensemble_intialization as init_module
from ICESEE.src.parallelization._mpi_ensemble_intialization import (
    _member_initialization_context,
    generate_initial_member_increment,
)
from ICESEE.src.parallelization._parallel_i_o import (
    finalize_analysis_ensemble,
)
from ICESEE.src.utils.localization import (
    compute_X5_from_matrices,
    generated_observation_terms,
    iter_generated_observation_columns,
    stochastic_observation_terms,
)
from ICESEE.src.utils.tools import _list_sorted_files, icesee_get_index


class _SingleRankComm:
    def Get_rank(self):
        return 0

    def Get_size(self):
        return 1

    def Barrier(self):
        return None

    def bcast(self, value, root=0):
        return value

    def allgather(self, value):
        return [value]


def test_hdf5_compression_configuration_normalizes_yaml_none():
    assert normalize_hdf5_compression(None) is None
    assert normalize_hdf5_compression("None") is None
    assert normalize_hdf5_compression(" null ") is None
    assert normalize_hdf5_compression("LZF") == "lzf"
    assert normalize_hdf5_compression("gzip") == "gzip"
    assert hdf5_compression_options(None, 4) is None
    assert hdf5_compression_options("lzf", 4) is None
    assert hdf5_compression_options("gzip", 4) == 4


def test_generated_observation_terms_are_mode_and_order_invariant(monkeypatch):
    """Mode 1 and streamed mode 2 must construct identical generated_R."""
    from ICESEE.src.run_model_da import _error_generation

    def fake_field(**icesee_kwargs):
        return icesee_kwargs["rng"].standard_normal(
            int(icesee_kwargs["noise_dim"])
        )

    monkeypatch.setattr(_error_generation, "generate_enkf_field", fake_field)
    settings = {
        "base_seed": 731,
        "total_state_param_vars": 2,
        "num_state_vars": 2,
        "joint_estimation": False,
        "localization_flag": False,
        "sig_obs": [2.0, 5.0],
        "Lx": 10.0,
        "Ly": 4.0,
    }
    obs_indices = np.array([0, 3, 4, 7], dtype=np.int64)
    forecast_obs = np.arange(20, dtype=float).reshape(4, 5) / 7.0
    observations = np.array([1.0, -2.0, 0.5, 4.0])

    mode1 = generated_observation_terms(
        forecast_obs, observations, settings, obs_indices, 8, 3
    )
    raw_mode2 = np.empty_like(forecast_obs)
    for member, column in iter_generated_observation_columns(
        settings, obs_indices, forecast_obs.shape[1], 8, 3
    ):
        raw_mode2[:, member] = column
    raw_mode2 -= raw_mode2.mean(axis=1, keepdims=True)

    np.testing.assert_array_equal(mode1[1], raw_mode2)
    np.testing.assert_array_equal(mode1[0], forecast_obs - forecast_obs.mean(
        axis=1, keepdims=True
    ))
    np.testing.assert_array_equal(
        mode1[2], observations[:, None] + raw_mode2 - forecast_obs
    )
    np.testing.assert_allclose(
        raw_mode2.mean(axis=1), 0.0, rtol=0.0, atol=1.0e-15
    )

    repeated = generated_observation_terms(
        forecast_obs, observations, settings, obs_indices, 8, 3
    )
    for expected, actual in zip(mode1, repeated):
        np.testing.assert_array_equal(expected, actual)


class _DummyInverseModel:
    @staticmethod
    def inverse_step_single(ensemble=None, **icesee_kwargs):
        indices = _inversion_variable_info(
            ensemble.size, icesee_kwargs
        )["indices"]
        return {"coefficient": ensemble[indices] + 10.0}


def test_sorted_state_shards_exclude_aggregate_files(tmp_path):
    for name in (
        "icesee_enkf_ens_0010.h5",
        "icesee_enkf_ens_0002.h5",
        "icesee_enkf_ens_mean.h5",
    ):
        with h5py.File(tmp_path / name, "w") as handle:
            handle["states"] = np.zeros((2, 2))

    assert [path.rsplit("/", 1)[-1] for path in _list_sorted_files(tmp_path)] == [
        "icesee_enkf_ens_0002.h5",
        "icesee_enkf_ens_0010.h5",
    ]


class _FailingInverseModel:
    @staticmethod
    def inverse_step_single(ensemble=None, **icesee_kwargs):
        raise RuntimeError("synthetic inversion failure")


def test_equal_size_indexing_does_not_require_vec_inputs_old():
    _, index_map, local_size = icesee_get_index(
        vec_inputs=["x", "y", "z"],
        nd=3,
        dim_list=[3],
        default_run=False,
        comm_world=None,
        even_distribution=False,
    )

    np.testing.assert_array_equal(index_map["x"], [0])
    np.testing.assert_array_equal(index_map["y"], [1])
    np.testing.assert_array_equal(index_map["z"], [2])
    assert local_size == 3


def test_equal_size_indexing_recovers_missing_dim_list():
    _, index_map, local_size = icesee_get_index(
        vec_inputs=["x", "y", "z"],
        nd=3,
        default_run=False,
        comm_world=None,
        even_distribution=False,
    )

    np.testing.assert_array_equal(index_map["x"], [0])
    np.testing.assert_array_equal(index_map["y"], [1])
    np.testing.assert_array_equal(index_map["z"], [2])
    assert local_size == 3


def _legacy_observation_space_transform(forecast, observations, energy=0.999):
    anomalies = forecast - forecast.mean(axis=1, keepdims=True)
    innovation = observations[:, None] - forecast
    combined = 2.0 * anomalies
    u, singular, _ = np.linalg.svd(combined, full_matrices=False)
    eigenvalues = singular**2
    total = eigenvalues.sum()
    inverse = np.zeros_like(eigenvalues)
    subtotal = 0.0
    for index, value in enumerate(eigenvalues):
        if total > 0 and subtotal / total < energy:
            inverse[index] = 1.0 / value
            subtotal += value
        else:
            break
    x3 = u @ ((inverse[:, None] * u.T) @ innovation)
    return np.eye(forecast.shape[1]) + anomalies.T @ x3


def test_streamed_ensemble_transform_matches_observation_space_svd():
    rng = np.random.default_rng(42)
    forecast = rng.normal(size=(97, 12))
    observations = rng.normal(size=97)
    anomalies = forecast - forecast.mean(axis=1, keepdims=True)
    innovation = observations[:, None] - forecast
    gram = 4.0 * anomalies.T @ anomalies
    rhs = 2.0 * anomalies.T @ innovation

    expected = _legacy_observation_space_transform(forecast, observations)
    actual = legacy_transform_from_gram(gram, rhs)
    np.testing.assert_allclose(actual, expected, rtol=2e-11, atol=2e-11)


def test_general_streamed_transform_matches_legacy_observation_space_svd():
    rng = np.random.default_rng(52)
    forecast = rng.normal(size=(83, 10))
    observations = rng.normal(size=83)
    yprime = forecast - forecast.mean(axis=1, keepdims=True)
    eta = yprime.copy()
    innovation = observations[:, None] - forecast
    factor = yprime + eta

    expected = compute_X5_from_matrices(
        yprime, eta, innovation, forecast.shape[1]
    )
    actual = ensemble_transform_from_products(
        yprime.T @ factor,
        factor.T @ factor,
        factor.T @ innovation,
    )
    np.testing.assert_allclose(actual, expected, rtol=3e-11, atol=3e-11)


def test_general_streamed_transform_matches_stochastic_observation_space_svd():
    rng = np.random.default_rng(62)
    forecast = rng.normal(size=(89, 11))
    observations = rng.normal(size=89)
    sigma = np.linspace(0.1, 0.8, observations.size)
    yprime, eta, innovation = stochastic_observation_terms(
        forecast, observations, sigma, seed=112233
    )
    factor = yprime + eta

    expected = compute_X5_from_matrices(
        yprime, eta, innovation, forecast.shape[1]
    )
    actual = ensemble_transform_from_products(
        yprime.T @ factor,
        factor.T @ factor,
        factor.T @ innovation,
    )
    np.testing.assert_allclose(actual, expected, rtol=4e-11, atol=4e-11)


def test_streamed_products_are_chunk_invariant():
    rng = np.random.default_rng(9)
    forecast = rng.normal(size=(103, 8))
    observations = rng.normal(size=103)
    expected_gram = np.zeros((8, 8))
    expected_rhs = np.zeros((8, 8))
    for start in range(0, 103, 11):
        block = forecast[start:start + 11]
        anomaly = block - block.mean(axis=1, keepdims=True)
        innovation = observations[start:start + 11, None] - block
        expected_gram += 4.0 * anomaly.T @ anomaly
        expected_rhs += 2.0 * anomaly.T @ innovation
    full_anomaly = forecast - forecast.mean(axis=1, keepdims=True)
    np.testing.assert_allclose(expected_gram, 4.0 * full_anomaly.T @ full_anomaly)
    np.testing.assert_allclose(
        expected_rhs, 2.0 * full_anomaly.T @ (observations[:, None] - forecast)
    )


def test_streamed_analysis_controls_match_reference_analysis():
    """Mode 2 must reproduce mode-1 controls for every row partition."""
    rng = np.random.default_rng(17)
    nd, nens = 36, 9
    forecast = rng.normal(size=(nd, nens))
    transform = np.eye(nens) + 0.025 * rng.normal(size=(nens, nens))
    raw_analysis = forecast @ transform
    settings = {
        "nd": nd,
        "vec_inputs": ["thickness", "velocity", "bed"],
        "num_state_vars": 2,
        "analysis_relaxation_factor": 0.85,
        "analysis_relaxation_factors": {"velocity": 0.65},
        "analysis_increment_limits": {"thickness": 0.2, "bed": 0.15},
        "inflation_factor": 1.01,
        "state_inflation_factor": 1.02,
        "param_inflation_factor": 1.03,
        "bed_inflation_factor": 1.04,
        "frozen_analysis_vars": ["velocity"],
        "inference_plugin_enabled": True,
        "bed_update_mode": "snapshots",
        "bed_snap_cols": [1],
        "km": 0,
    }
    expected = apply_analysis_controls_local(
        raw_analysis.copy(), forecast.copy(), np.arange(nd), dict(settings)
    )

    actual = np.empty_like(expected)
    for start, stop in ((0, 5), (5, 19), (19, 31), (31, nd)):
        actual[start:stop] = apply_analysis_controls_local(
            raw_analysis[start:stop].copy(),
            forecast[start:stop].copy(),
            np.arange(start, stop),
            dict(settings),
        )

    np.testing.assert_allclose(actual, expected, rtol=0.0, atol=0.0)


def test_row_reader_preserves_order_and_duplicates(tmp_path):
    path = tmp_path / "rows.h5"
    values = np.arange(40).reshape(10, 4)
    with h5py.File(path, "w") as handle:
        handle["x"] = values
    with h5py.File(path, "r") as handle:
        actual = _read_rows_preserve_order(handle["x"], [7, 2, 7, 0])
    np.testing.assert_array_equal(actual, values[[7, 2, 7, 0]])


def test_chunk_size_honors_budget_and_explicit_override():
    assert resolve_row_chunk_size(10000, 40, budget_mb=1) < 10000
    assert resolve_row_chunk_size(1000, 40, explicit_rows=73) == 73
    assert resolve_row_chunk_size(20, 40, explicit_rows=73) == 20


def test_process_noise_schedule_is_execution_mode_independent():
    settings = {
        "process_noise_schedule": "observations",
        "obs_index": [2, 5, 9],
        "number_obs_instants": 3,
    }
    assert not process_noise_is_due(settings, 1)
    assert process_noise_is_due(settings, 2)
    assert not process_noise_is_due(settings, 3)
    assert process_noise_is_due(settings, 9)
    assert process_noise_is_due(
        {"process_noise_schedule": "every_step"}, 3
    )
    assert not process_noise_is_due(
        {"process_noise_schedule": "none"}, 2
    )


def test_member_process_noise_is_order_and_restart_invariant(tmp_path, monkeypatch):
    def fake_field(**kwargs):
        return np.random.standard_normal(int(kwargs["noise_dim"]))

    monkeypatch.setattr(forecast_module, "generate_enkf_field", fake_field)
    common = {
        "num_state_vars": 2,
        "total_state_param_vars": 3,
        "joint_estimation": True,
        "localization_flag": False,
        "sig_Q": [0.4, 0.2, 0.0],
        "alpha": 0.75,
        "rho": 1.1,
        "dt": 0.2,
        "Lx": 10.0,
        "Ly": 4.0,
        "base_seed": 927,
        "process_noise_schedule": "every_step",
        "sub_rank": 0,
    }

    first_dir = tmp_path / "forward_order"
    first = {}
    for member in (0, 1):
        settings = dict(common, _modelrun_datasets=str(first_dir), k=0)
        first[member] = add_member_process_noise(
            np.zeros(12), member, settings
        )
    continued = {}
    for member in (1, 0):
        settings = dict(common, _modelrun_datasets=str(first_dir), k=1)
        continued[member] = add_member_process_noise(
            first[member].copy(), member, settings
        )

    second_dir = tmp_path / "reverse_order"
    reverse_first = {}
    for member in (1, 0):
        settings = dict(common, _modelrun_datasets=str(second_dir), k=0)
        reverse_first[member] = add_member_process_noise(
            np.zeros(12), member, settings
        )
    reverse_continued = {}
    for member in (0, 1):
        settings = dict(common, _modelrun_datasets=str(second_dir), k=1)
        reverse_continued[member] = add_member_process_noise(
            reverse_first[member].copy(), member, settings
        )

    for member in (0, 1):
        np.testing.assert_array_equal(first[member], reverse_first[member])
        np.testing.assert_array_equal(
            continued[member], reverse_continued[member]
        )
    assert not np.array_equal(first[0], first[1])


def test_member_initialization_is_order_and_execution_mode_invariant(monkeypatch):
    """Dense and streamed initializers must use the same member fields."""
    def fake_field(**kwargs):
        return np.random.standard_normal(int(kwargs["noise_dim"]))

    monkeypatch.setattr(init_module, "generate_enkf_field", fake_field)
    settings = {
        "base_seed": 814,
        "total_state_param_vars": 3,
        "sig_Q": [0.4, 0.2, 0.1],
        "Lx": 10.0,
        "Ly": 4.0,
    }
    forward = {
        member: generate_initial_member_increment(4, settings, member, 12)[0]
        for member in (0, 1, 2)
    }
    reverse = {
        member: generate_initial_member_increment(4, settings, member, 12)[0]
        for member in (2, 1, 0)
    }
    for member in forward:
        np.testing.assert_array_equal(forward[member], reverse[member])
    assert not np.array_equal(forward[0], forward[1])

    # A reduced active state (for example hybrid analysis with frozen
    # friction) consumes only its active blocks, not configured padding.
    reduced, _ = generate_initial_member_increment(4, settings, 0, 8)
    np.testing.assert_array_equal(reduced, forward[0][:8])


def test_application_rng_context_is_member_keyed_and_order_invariant():
    settings = {"base_seed": 517, "k": 9}

    def sample(member):
        with forecast_member_context(settings, member) as member_kwargs:
            return (
                np.random.standard_normal(5),
                member_kwargs["rng"].standard_normal(5),
                member_kwargs["rank_seed"],
            )

    forward = {member: sample(member) for member in (0, 1, 2)}
    reverse = {member: sample(member) for member in (2, 1, 0)}
    for member in forward:
        np.testing.assert_array_equal(forward[member][0], reverse[member][0])
        np.testing.assert_array_equal(forward[member][1], reverse[member][1])
        assert forward[member][2] == reverse[member][2]
    assert not np.array_equal(forward[0][0], forward[1][0])


def test_model_initializer_context_is_member_keyed_and_order_invariant():
    settings = {"base_seed": 912}

    def sample(member):
        with _member_initialization_context(settings, member) as member_kwargs:
            return (
                np.random.standard_normal(5),
                member_kwargs["rng"].standard_normal(5),
            )

    forward = {member: sample(member) for member in (0, 1, 2)}
    reverse = {member: sample(member) for member in (2, 1, 0)}
    for member in forward:
        np.testing.assert_array_equal(forward[member][0], reverse[member][0])
        np.testing.assert_array_equal(forward[member][1], reverse[member][1])


def test_model_initializer_context_exposes_shape_without_dense_ensemble():
    settings = {"base_seed": 912, "nd": 2_000_000, "Nens": 400}
    with _member_initialization_context(settings, 7) as member_kwargs:
        proxy = member_kwargs["statevec_ens"]
        assert proxy.shape == (2_000_000, 400)
        assert not isinstance(proxy, np.ndarray)


def test_hybrid_cycles_exclude_friction_observations():
    settings = {
        "inversion_flag": True,
        "vec_inputs_old": ["thickness", "velocity", "coefficient"],
        "friction_idx": 2,
    }
    info = _inversion_variable_info(12, settings)
    assert (info["start"], info["stop"]) == (8, 12)
    observed_rows = np.array([0, 3, 8, 9, 11, 6])
    np.testing.assert_array_equal(
        _exclude_inversion_observations(observed_rows, settings, 12),
        [True, True, False, False, False, True],
    )


def _transaction_io(tmp_path):
    io = EnKF_fully_parallel_IO.__new__(EnKF_fully_parallel_IO)
    io.nd = 6
    io.nens = 3
    io.size = 1
    io.rank = 0
    io.mpi_comm = _SingleRankComm()
    io.base_path = str(tmp_path)
    io.file_prefix = "icesee_enkf_ens"
    io.storage_dtype = np.dtype("float64")
    io.files = []
    io.datasets = []
    io.current_batch_start = -1
    io._close_batch = lambda: None
    return io


def test_hybrid_inversion_is_applied_before_atomic_publish(tmp_path):
    io = _transaction_io(tmp_path)
    canonical = tmp_path / "icesee_enkf_ens_0001.h5"
    temporary = tmp_path / "icesee_enkf_ens_0001.h5.analysis.tmp"
    forecast = np.arange(18, dtype=float).reshape(6, 3)
    analysis = forecast.copy()
    analysis[:3] += 0.5
    with h5py.File(canonical, "w") as handle:
        handle["states"] = forecast
    with h5py.File(temporary, "w") as handle:
        handle["states"] = analysis
    settings = {
        "inversion_flag": True,
        "vec_inputs_old": ["thickness", "coefficient"],
        "friction_idx": 1,
        "model_module": _DummyInverseModel(),
    }

    io._publish_analysis_file(1, str(temporary), settings)

    with h5py.File(canonical, "r") as handle:
        published = np.asarray(handle["states"])
    np.testing.assert_array_equal(published[:3], analysis[:3])
    np.testing.assert_array_equal(published[3:], forecast[3:] + 10.0)
    assert not temporary.exists()


def test_failed_inversion_does_not_replace_last_valid_forecast(tmp_path):
    io = _transaction_io(tmp_path)
    canonical = tmp_path / "icesee_enkf_ens_0001.h5"
    temporary = tmp_path / "icesee_enkf_ens_0001.h5.analysis.tmp"
    forecast = np.arange(18, dtype=float).reshape(6, 3)
    with h5py.File(canonical, "w") as handle:
        handle["states"] = forecast
    with h5py.File(temporary, "w") as handle:
        handle["states"] = forecast + 100.0
    settings = {
        "inversion_flag": True,
        "vec_inputs_old": ["thickness", "coefficient"],
        "friction_idx": 1,
        "model_module": _FailingInverseModel(),
    }

    with np.testing.assert_raises_regex(RuntimeError, "synthetic inversion failure"):
        io._publish_analysis_file(1, str(temporary), settings)

    with h5py.File(canonical, "r") as handle:
        np.testing.assert_array_equal(handle["states"][:], forecast)
    assert temporary.exists()


def test_member_slab_issm_finalization_matches_mode1_full_ensemble(tmp_path):
    """Bounded mode-2 ISSM finalization must be numerically identical."""
    io = _transaction_io(tmp_path)
    io.nd = 12
    io.nens = 5
    canonical = tmp_path / "icesee_enkf_ens_0002.h5"
    temporary = tmp_path / "icesee_enkf_ens_0002.h5.analysis.tmp"
    rng = np.random.default_rng(91)
    forecast = np.empty((12, 5))
    forecast[0:4] = rng.uniform(150.0, 500.0, size=(4, 5))
    forecast[4:8] = rng.uniform(-150.0, 350.0, size=(4, 5))
    forecast[8:12] = rng.uniform(-900.0, -100.0, size=(4, 5))
    analysis = forecast + rng.normal(0.0, 35.0, size=forecast.shape)
    settings = {
        "model_name": "issm",
        "vec_inputs": ["thickness", "surface", "bed"],
        "dt": 0.2,
        "di": 0.893,
        "rho_ice": 917.0,
        "rho_sw": 1028.0,
        "geometry_projection_mode": "preserve_thickness",
        "analysis_finalize_memory_budget_mb": 1.0,
        "analysis_finalize_member_chunk_size": 1,
    }
    expected = finalize_analysis_ensemble(
        analysis.copy(), forecast.copy(), 2, dict(settings)
    )
    with h5py.File(canonical, "w") as handle:
        handle["states"] = forecast
    with h5py.File(temporary, "w") as handle:
        handle["states"] = analysis

    io._publish_analysis_file(2, str(temporary), settings)

    with h5py.File(canonical, "r") as handle:
        actual = np.asarray(handle["states"])
    np.testing.assert_allclose(actual, expected, rtol=0.0, atol=0.0)


def test_hybrid_transaction_matches_mode1_reduced_state_semantics():
    """Mode 2's fixed-shape hybrid update must equal mode 1 reinjection."""
    rng = np.random.default_rng(314)
    hdim, nens = 4, 6
    full_names = ["thickness", "surface", "bed", "coefficient"]
    reduced_names = full_names[:-1]
    forecast = np.empty((hdim * len(full_names), nens))
    forecast[0:hdim] = rng.uniform(150.0, 500.0, (hdim, nens))
    forecast[hdim:2 * hdim] = rng.uniform(-100.0, 300.0, (hdim, nens))
    forecast[2 * hdim:3 * hdim] = rng.uniform(-800.0, -50.0, (hdim, nens))
    forecast[3 * hdim:4 * hdim] = rng.uniform(1800.0, 4200.0, (hdim, nens))
    transform = np.eye(nens) + 0.015 * rng.normal(size=(nens, nens))
    common = {
        "model_name": "issm",
        "num_state_vars": 2,
        "dt": 0.2,
        "di": 0.893,
        "rho_ice": 917.0,
        "rho_sw": 1028.0,
        "geometry_projection_mode": "preserve_thickness",
        "analysis_relaxation_factor": 0.9,
        "analysis_increment_limits": {"bed": 12.0},
        "state_inflation_factor": 1.01,
        "param_inflation_factor": 1.02,
        "bed_inflation_factor": 1.03,
    }

    # Mode 1 removes friction before EnKF analysis, finalizes the reduced
    # state, reinjects it into the full forecast, then runs inversion.
    reduced_forecast = forecast[:3 * hdim].copy()
    reduced_settings = dict(
        common, nd=3 * hdim, vec_inputs=reduced_names,
        total_state_param_vars=3,
    )
    mode1 = forecast.copy()
    mode1[:3 * hdim] = apply_analysis_controls_local(
        reduced_forecast @ transform,
        reduced_forecast,
        np.arange(3 * hdim),
        reduced_settings,
    )
    mode1[:3 * hdim] = finalize_analysis_ensemble(
        mode1[:3 * hdim], reduced_forecast, 1, dict(reduced_settings)
    )
    mode1[3 * hdim:] = forecast[3 * hdim:] + 10.0

    # Mode 2 keeps the full on-disk shape, freezes friction for the EnKF
    # update, uses the same finalizer, then applies the same inversion result.
    full_settings = dict(
        common, nd=4 * hdim, vec_inputs=full_names,
        total_state_param_vars=4, frozen_analysis_vars=["coefficient"],
    )
    mode2 = apply_analysis_controls_local(
        forecast @ transform,
        forecast,
        np.arange(4 * hdim),
        full_settings,
    )
    mode2 = finalize_analysis_ensemble(
        mode2, forecast, 1, dict(full_settings)
    )
    mode2[3 * hdim:] = forecast[3 * hdim:] + 10.0

    np.testing.assert_allclose(mode2, mode1, rtol=0.0, atol=0.0)


def test_auto_history_selects_rolling_when_full_history_cannot_fit(tmp_path):
    comm = _SingleRankComm()
    io = EnKF_fully_parallel_IO(
        "icesee_enkf_ens",
        nd=10**10,
        nens=40,
        nt=100,
        subcomm=comm,
        mpi_comm=comm,
        icesee_kwargs={
            "ensemble_history_mode": "auto",
            "ensemble_storage_dtype": "float64",
            "restart_enabled": True,
            "fail_on_insufficient_storage": False,
        },
        serial_file_creation=True,
        base_path=str(tmp_path),
        batch_size=2,
    )
    assert io.nt == 101
    assert io.history_mode == "rolling"


def test_rolling_history_prunes_numeric_state_shards(tmp_path):
    comm = _SingleRankComm()
    io = EnKF_fully_parallel_IO.__new__(EnKF_fully_parallel_IO)
    io.history_mode = "rolling"
    io.files = []
    io.datasets = []
    io.current_batch_start = 0
    io.base_path = str(tmp_path)
    io.file_prefix = "icesee_enkf_ens"
    io.mpi_comm = comm
    io.rank = 0
    for timestep in range(3):
        with h5py.File(tmp_path / f"icesee_enkf_ens_{timestep:04d}.h5", "w") as handle:
            handle["states"] = np.zeros((2, 2))
    with h5py.File(tmp_path / "icesee_enkf_ens_mean.h5", "w") as handle:
        handle["mean"] = np.zeros((2, 3))

    io.prune_history(2)

    assert not (tmp_path / "icesee_enkf_ens_0000.h5").exists()
    assert not (tmp_path / "icesee_enkf_ens_0001.h5").exists()
    assert (tmp_path / "icesee_enkf_ens_0002.h5").exists()
    assert (tmp_path / "icesee_enkf_ens_mean.h5").exists()
