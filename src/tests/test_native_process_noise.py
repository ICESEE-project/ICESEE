from types import SimpleNamespace
import numpy as np
import pytest
from mpi4py import MPI

from ICESEE.src.parallelization.native_process_noise import NativeProcessNoiseAdapter
from ICESEE.src.parallelization import _mpi_forecast_functions as canonical
from ICESEE.src.parallelization.distributed_native_runtime import NativeDistributedMember
from ICESEE.src.parallelization.distributed_whole_member_adapter import build_whole_member_field_registry


def setup_adapter(monkeypatch):
    monkeypatch.setattr(canonical, 'generate_enkf_field',
                        lambda **kw: kw['rng'].standard_normal(kw['noise_dim']))
    options = dict(nd=9, nt=6, dt=.05, num_state_vars=2,
                   total_state_param_vars=3, joint_estimation=True,
                   localization_flag=False, Lx=100., Ly=20.,
                   sig_Q=[.8,.4,.2], base_seed=42,
                   obs_index=np.array([1,3,5]), number_obs_instants=3)
    topology = SimpleNamespace(spatial_comm=MPI.COMM_SELF)
    underlying = SimpleNamespace(forecast_native_member=lambda *a,**kw:None,
                                 restore_native_checkpoint=object())
    rows = np.array([2,0,1,5,3,4,8,6,7])
    wrapped = NativeProcessNoiseAdapter(
        underlying, topology=topology, canonical_rows=rows,
        reference_coords=np.arange(3.)[:,None], icesee_kwargs=options)
    return wrapped, topology, options, rows


def member(member_id=7, values=None):
    values = np.zeros(9) if values is None else values
    fields = build_whole_member_field_registry(
        {name:values[i*3:(i+1)*3].copy()
         for i,name in enumerate(('h','u','parameter'))})
    return NativeDistributedMember(member_id,fields)


def test_native_noise_matches_canonical_fields_and_preserves_parameters(monkeypatch,tmp_path):
    wrapped, topology, options, rows = setup_adapter(monkeypatch)
    physical = member()
    expected = np.zeros(9)
    baseline_options = dict(wrapped.options, data_path=str(tmp_path),
                            _modelrun_datasets=str(tmp_path))
    try:
        assert wrapped.restore_native_checkpoint is wrapped.adapter.restore_native_checkpoint
        for k in range(6):
            baseline_options['k'] = k
            expected = canonical.add_member_process_noise(expected,7,baseline_options)
            wrapped.forecast_native_member(physical,k,topology=topology,
                                           icesee_kwargs=options)
            np.testing.assert_array_equal(physical.pack_owned(),expected[rows])
            np.testing.assert_array_equal(physical.pack_owned()[6:],np.zeros(3))
        assert np.any(expected[:6])
    finally:
        wrapped.close()


def test_checkpoint_resume_replays_ar_history_without_perturbing_initial_fields(monkeypatch):
    first,topology,options,rows = setup_adapter(monkeypatch)
    uninterrupted = member()
    try:
        np.testing.assert_array_equal(uninterrupted.pack_owned(),np.zeros(9))
        for k in range(2):
            first.forecast_native_member(uninterrupted,k,topology=topology,icesee_kwargs=options)
        saved = uninterrupted.pack_owned().copy()
        for k in range(2,6):
            first.forecast_native_member(uninterrupted,k,topology=topology,icesee_kwargs=options)
        restored,_,_,_ = setup_adapter(monkeypatch)
        resumed = member(values=saved)
        try:
            # Constructing/restoring the policy introduces no initial noise.
            np.testing.assert_array_equal(resumed.pack_owned(),saved)
            for k in range(2,6):
                restored.forecast_native_member(resumed,k,topology=topology,icesee_kwargs=options)
            np.testing.assert_array_equal(resumed.pack_owned(),uninterrupted.pack_owned())
        finally:
            restored.close()
    finally:
        first.close()


def test_workspace_limit_refuses_before_allocating_a_global_increment():
    with pytest.raises(ValueError,match='workspace budget'):
        NativeProcessNoiseAdapter(object(),topology=SimpleNamespace(spatial_comm=MPI.COMM_SELF),
                                  canonical_rows=np.arange(3),reference_coords=np.zeros((1,1)),
                                  icesee_kwargs={'nd':3},workspace_budget_mb=0)


def test_ownership_map_cannot_duplicate_or_omit_a_reference_row():
    with pytest.raises(ValueError,match='exactly once'):
        NativeProcessNoiseAdapter(object(),topology=SimpleNamespace(spatial_comm=MPI.COMM_SELF),
                                  canonical_rows=np.array([0,0,2]),reference_coords=np.zeros((1,1)),
                                  icesee_kwargs={'nd':3})


def test_configured_noise_precedes_first_analysis_without_initial_perturbation(monkeypatch):
    wrapped, topology, options, rows = setup_adapter(monkeypatch)
    options.update(obs_index=np.array([10,15,20]), number_obs_instants=3)
    members = [member(0),member(1)]
    try:
        for k in range(9):
            for physical in members:
                wrapped.forecast_native_member(physical,k,topology=topology,icesee_kwargs=options)
        for physical in members:
            np.testing.assert_array_equal(physical.pack_owned(),np.zeros(9))
        for physical in members:
            wrapped.forecast_native_member(physical,9,topology=topology,icesee_kwargs=options)
        assert np.any(members[0].pack_owned()[:6] != members[1].pack_owned()[:6])
        assert [k for k in range(20) if canonical.process_noise_is_due(options,k)] == [9,14,19]
    finally:
        wrapped.close()


def test_noise_checkpoint_identity_includes_time_contract(monkeypatch):
    from ICESEE.src.parallelization import native_process_noise as native
    original=native.canonical_noise_identity({})
    assert original.startswith('canonical-noise-v2-')
    monkeypatch.setattr(native,'TIME_GRID_CONTRACT','other-clock')
    assert native.canonical_noise_identity({}) != original
