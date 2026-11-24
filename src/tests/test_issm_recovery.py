import json
from pathlib import Path
import re
from types import SimpleNamespace as NS
import h5py
import numpy as np
import pytest

from ICESEE.applications.issm_model.issm_utils.recovery import (
    MatlabMemberTransaction, checked_command, NUMERICAL_REJECTION,
    publish_scalar_diagnostics, commit_backend, restore_backend_member,
    recovery_identity, prune_backend,
)
from ICESEE.src.parallelization.forecast_recovery import RecoverableForecastError


class Server:
    def __init__(self):
        self.rng = 10
        self.error = None
    def send_command(self, wrapped, timeout):
        match = re.match(r"issm_recovery_command\('((?:[^']|'')*)', '((?:[^']|'')*)'\);",wrapped)
        command,errorfile = [v.replace("''", "'") for v in match.groups()]
        if self.error:
            Path(errorfile).write_text(json.dumps(dict(identifier=self.error,message="test error")))
            return False
        rng = re.match(r"issm_recovery_rng\('(save|restore)', '((?:[^']|'')*)'\);",command)
        if rng:
            action,filename = rng.groups()
            filename = Path(filename.replace("''", "'"))
            if action == "save":filename.write_text(str(self.rng))
            else:self.rng = int(filename.read_text())
        return True


@pytest.fixture
def options(tmp_path):
    workspace=tmp_path/'workspace';workspace.mkdir()
    data=tmp_path/'data';data.mkdir()
    return dict(issm_examples_dir=str(workspace),icesee_path=str(tmp_path),data_path='data',
                server=Server(),Nens=2,nt=2,scalar_inputs=['IceVolume'])


def test_full_model_outputs_and_matlab_rng_are_restored(options):
    transaction=MatlabMemberTransaction(options,0)
    transaction.models.mkdir(parents=True)
    model=transaction.models/'enkf_state.mat';model.write_bytes(b'complete model before step')
    output=transaction.data/'ensemble_output_0.h5';output.write_bytes(b'previous output')
    snapshot=transaction.capture(transaction.data/'snapshot')
    model.write_bytes(b'failed model');output.write_bytes(b'failed output')
    (transaction.models/'inverse_state.mat').write_bytes(b'failed auxiliary state')
    options['server'].rng=991
    transaction.restore(snapshot)
    assert model.read_bytes()==b'complete model before step'
    assert output.read_bytes()==b'previous output'
    assert not (transaction.models/'inverse_state.mat').exists()
    assert options['server'].rng==10


def test_corrupt_snapshot_is_refused_before_replacing_live_files(options):
    transaction=MatlabMemberTransaction(options,0)
    transaction.models.mkdir(parents=True)
    model=transaction.models/'enkf_state.mat';model.write_bytes(b'old')
    snapshot=transaction.capture(transaction.data/'snapshot')
    (snapshot/'models/enkf_state.mat').write_bytes(b'corrupt')
    model.write_bytes(b'live')
    with pytest.raises(ValueError,match='size mismatch'):
        transaction.restore(snapshot)
    assert model.read_bytes()==b'live'


def test_snapshot_budget_and_symlinks_fail_before_trial(options):
    transaction=MatlabMemberTransaction(dict(options,issm_recovery_snapshot_mb=.000001),0)
    transaction.models.mkdir(parents=True)
    (transaction.models/'enkf_state.mat').write_bytes(b'x'*30)
    with pytest.raises(MemoryError):transaction.capture(transaction.data/'snapshot')
    assert not (transaction.data/'snapshot').exists()
    (transaction.models/'enkf_state.mat').unlink()
    (transaction.models/'enkf_state.mat').symlink_to(transaction.data/'external.mat')
    with pytest.raises(ValueError,match='symlink'):transaction.capture(transaction.data/'snapshot')


def test_only_explicit_numerical_rejection_is_retryable(options):
    options['server'].error=NUMERICAL_REJECTION
    with pytest.raises(RecoverableForecastError):checked_command(options['server'],'solve_model();',options)
    options['server'].error='MATLAB:UndefinedFunction'
    with pytest.raises(RuntimeError,match='non-retryable'):checked_command(options['server'],'solve_model();',options)


def test_scalar_means_wait_for_all_accepted_members_without_mpi(options):
    data=Path(options['icesee_path'])/'data'
    for member,value in enumerate((10.,20.)):
        with h5py.File(data/f'ensemble_out_scalar_{member}.h5','w') as h:
            h['IceVolume']=[value,value+1]
    publish_scalar_diagnostics(dict(options,ens_id=0),0)
    assert not (data/'ensemble_scalar_output.h5').exists()
    publish_scalar_diagnostics(dict(options,ens_id=1),0)
    with h5py.File(data/'ensemble_scalar_output.h5') as h:
        assert h['IceVolume'][0]==15
        assert np.isnan(h['IceVolume'][1])


def test_backend_manifest_commit_restores_model_cache(options):
    class Comm:
        def bcast(self,value,root):return value
        def allgather(self,value):return [value]
        def Get_size(self):return 1
    topology=NS(world=Comm(),world_rank=0)
    pool=NS(member_ids=[0,1])
    for member in pool.member_ids:
        model=Path(options['issm_examples_dir'])/f'Models/ens_id_{member}'
        model.mkdir(parents=True);(model/'enkf_state.mat').write_bytes(f'member {member}'.encode())
    descriptor=commit_backend(pool,topology,options,0)
    changed=Path(options['issm_examples_dir'])/'Models/ens_id_1/enkf_state.mat'
    changed.write_bytes(b'newer uncommitted model')
    checkpoint=NS(metadata={'issm_backend':descriptor})
    restore_backend_member(1,checkpoint,topology,options)
    assert changed.read_bytes()==b'member 1'
    topology.world.Get_size=lambda:2
    with pytest.raises(ValueError,match='rank count'):
        restore_backend_member(1,checkpoint,topology,options)


def test_backend_identity_preserves_science_and_ignores_resume_controls(options):
    science = dict(options, dt=.2, seed=np.int64(42), sig_Q=np.array([1.,2.]))
    original = recovery_identity(science)
    assert recovery_identity(dict(science, resume_from_checkpoint=True,
        initialize_ensemble=False, run_segment_cycles=1, rank=3)) == original
    assert recovery_identity(dict(science, dt=.1)) != original
    assert recovery_identity(dict(science, seed=np.int64(43))) != original
    assert recovery_identity(dict(science, sig_Q=np.array([1.,3.]))) != original


def test_backend_retention_preserves_referenced_and_incomplete_bundles(options):
    data = Path(options['icesee_path'])/'data'
    initial,steps = data/'initial', data/'steps'
    initial.mkdir();steps.mkdir()
    bundles = data/'_issm_backend'
    names = ['step_00000001_kept','step_00000002_orphan','step_00000003_incomplete']
    for name in names:
        path = bundles/name;path.mkdir(parents=True)
        if not name.endswith('incomplete'):
            (path/'manifest.json').write_text(json.dumps(
                dict(format='icesee-issm-backend-v1',complete=True)))
    checkpoint = steps/'checkpoint_00000000';checkpoint.mkdir()
    (checkpoint/'manifest.json').write_text(json.dumps(dict(metadata=dict(
        issm_backend=dict(path='_issm_backend/'+names[0])))))
    manager = NS(initial_root=initial,steps_root=steps,policy=NS(keeps_all=True))
    prune_backend(manager,NS(world_rank=0),options)
    assert (bundles/names[1]).exists()
    manager.policy.keeps_all=False
    prune_backend(manager,NS(world_rank=0),options)
    assert (bundles/names[0]).exists()
    assert not (bundles/names[1]).exists()
    assert (bundles/names[2]).exists()
