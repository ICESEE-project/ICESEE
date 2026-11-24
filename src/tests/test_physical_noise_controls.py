"""Physical time, variable selection and graph-coordinate noise contracts."""
import numpy as np
import pytest
from ICESEE.src.utils.process_noise import process_noise_event_coefficients, process_noise_coefficients
from ICESEE.src.parallelization import _mpi_forecast_functions as canonical
from ICESEE.src.run_model_da._error_generation import physical_graph_field, generate_pseudo_random_field_1d


def options(**changes):
    result = dict(dt=.05,nt=1640,process_noise_normalization="physical_ar1",
        process_noise_correlation_time=2.,process_noise_start="stationary",
        obs_index=[300,340,380],number_obs_instants=3)
    result.update(changes)
    return result


def test_elapsed_time_between_observation_events_and_stationary_first_draw():
    first, rho = process_noise_event_coefficients(options(),299)
    assert first == 0
    assert rho*np.sqrt(.05) == pytest.approx(1)
    assert process_noise_event_coefficients(options(),339)[0] == pytest.approx(np.exp(-1))
    assert process_noise_event_coefficients(options(process_noise_start="zero"),299)[0] == pytest.approx(np.exp(-15/2))
    assert process_noise_event_coefficients(options(process_noise_schedule="every_step"),1)[0] == pytest.approx(np.exp(-.05/2))
    assert process_noise_event_coefficients(options(nt=4000),339) == process_noise_event_coefficients(options(),339)


@pytest.mark.parametrize("changes", [dict(process_noise_correlation_time=-1),dict(process_noise_correlation_time=np.nan),dict(process_noise_start="guess")])
def test_physical_noise_requires_explicit_valid_controls(changes):
    with pytest.raises(ValueError): process_noise_coefficients(options(**changes))


def test_only_selected_variables_receive_restartable_noise(monkeypatch,tmp_path):
    monkeypatch.setattr(canonical,"generate_enkf_field",lambda **kw:kw['rng'].standard_normal(kw['noise_dim']))
    kw=options(nd=8,num_state_vars=2,total_state_param_vars=2,vec_inputs=['h','s'],
        process_noise_vars=['h'],sig_Q=[.8,.08],base_seed=42,joint_estimation=False,
        localization_flag=False,data_path=str(tmp_path),sub_rank=0,k=299)
    first=canonical.add_member_process_noise(np.zeros(8),3,kw)
    np.testing.assert_array_equal(first[4:],np.zeros(4))
    from pathlib import Path
    path=Path(canonical._process_noise_state_path(kw,3))
    saved=path.read_bytes()
    kw['k']=339
    continued=canonical.add_member_process_noise(first.copy(),3,kw)
    path.write_bytes(saved)
    repeated=canonical.add_member_process_noise(first.copy(),3,dict(kw,k=339))
    np.testing.assert_array_equal(continued,repeated)
    assert np.isfinite(continued).all() and np.any(first[:4])
    with pytest.raises(RuntimeError,match='controls'):
        canonical.add_member_process_noise(continued,3,dict(kw,k=379,process_noise_correlation_time=3.))


def test_physical_graph_length_changes_covariance_and_preserves_coordinate_identity():
    coords=np.arange(200.)[:,None]
    short=physical_graph_field(coords,2.,42,max_neighbors=64)
    long=physical_graph_field(coords,10.,42,max_neighbors=64)
    assert np.corrcoef(short[:-5],short[5:])[0,1] < np.corrcoef(long[:-5],long[5:])[0,1]
    order=np.random.default_rng(8).permutation(len(coords))
    reordered=physical_graph_field(coords[order],10.,42,max_neighbors=64)
    np.testing.assert_allclose(reordered,long[order],atol=2e-13)
    configured=generate_pseudo_random_field_1d(N=200,coords=coords,rh=10.,seed=42,
        method='graph',graph_physical_length=True,graph_max_neighbors=64)
    np.testing.assert_array_equal(configured,long)
    with pytest.raises(MemoryError,match='neighbors'):
        physical_graph_field(coords,100.,42,max_neighbors=10)


def test_fork_noise_origin_is_explicit_and_does_not_replay_new_noise_before_fork():
    kw=options(process_noise_schedule='every_step',process_noise_start_step=299)
    assert not canonical.process_noise_is_due(kw,298)
    assert canonical.process_noise_is_due(kw,299)
    assert process_noise_event_coefficients(kw,299)[0]==0.
    assert process_noise_event_coefficients(kw,300)[0]==pytest.approx(np.exp(-.05/2))


def test_scalar_control_noise_is_spatially_coherent_and_member_keyed(tmp_path):
    kw=options(nd=8,num_state_vars=2,total_state_param_vars=2,vec_inputs=['h','basal_melt_field'],
        process_noise_vars=['basal_melt_field'],process_noise_spatial_policy={'basal_melt_field':'global_scalar'},
        sig_Q=[0,.1],base_seed=42,joint_estimation=False,localization_flag=False,
        data_path=str(tmp_path),sub_rank=0,k=299)
    values=[]
    for mid in (0,1):
        vector=canonical.add_member_process_noise(np.zeros(8),mid,kw)
        np.testing.assert_array_equal(vector[:4],np.zeros(4))
        np.testing.assert_array_equal(vector[4:],np.repeat(vector[4],4))
        values.append(vector[4])
    assert values[0]!=values[1]
