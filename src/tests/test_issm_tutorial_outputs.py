"""The teaching launcher must not accept a caught application error as success."""
from pathlib import Path
import importlib.util
import json
import sys
import types

import h5py
import pytest


@pytest.fixture(params=['basal_friction_variation', 'ISMIP_Choi'])
def tutorial(request):
    source = (Path(__file__).resolve().parents[2] / 'applications/issm_model/examples' /
              request.param / 'tutorials/short_assimilation/run_case.py')
    spec = importlib.util.spec_from_file_location('issm_tutorial_' + request.param, source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_missing_output_fails_even_when_child_process_returns_zero(tutorial, tmp_path, monkeypatch):
    info = dict(variant='fixed', command=['fake-child'], output_path=str(tmp_path / 'output'),
                expected_state_columns=11)
    monkeypatch.setattr(tutorial, 'repository_root', lambda: tmp_path)
    monkeypatch.setattr(tutorial, 'inspect_case', lambda *args: (dict(info), []))
    monkeypatch.setattr(tutorial.subprocess, 'run', lambda *args, **kwargs:
                        types.SimpleNamespace(returncode=0))
    monkeypatch.setattr(sys, 'argv', ['run_case.py', '--variant', 'fixed',
                                    '--run-id', 'missing', '--record-root', str(tmp_path / 'records')])
    assert tutorial.main() == 1
    record = json.loads((tmp_path / 'records/missing/fixed/execution-status.json').read_text())
    assert record['simulation_exit_code'] == 0
    assert record['tutorial_exit_code'] == 1
    assert record['output_verification']['status'] == 'failed'
    assert 'Missing expected output: icesee_ensemble_data.h5' in record['output_verification']['issues']


def test_existing_wrong_shape_is_rejected(tutorial, tmp_path):
    with h5py.File(tmp_path / 'icesee_ensemble_data.h5', 'w') as f:
        f.create_dataset('ensemble', shape=(2, 4, 3))
        f.create_dataset('ensemble_mean', shape=(2, 3))
    (tmp_path / 'synthetic_obs.h5').touch()
    result = tutorial.verify_output_artifacts(dict(output_path=str(tmp_path), expected_state_columns=11))
    assert result['status'] == 'failed'
    assert result['issues']


def test_complete_shapes_are_only_structural_evidence(tutorial, tmp_path):
    with h5py.File(tmp_path / 'icesee_ensemble_data.h5', 'w') as f:
        f.create_dataset('ensemble', shape=(2, 4, 11))
        f.create_dataset('ensemble_mean', shape=(2, 11))
    (tmp_path / 'synthetic_obs.h5').touch()
    result = tutorial.verify_output_artifacts(dict(output_path=str(tmp_path), expected_state_columns=11))
    assert result['status'] == 'passed structural completion checks'
    assert result['issues'] == []
    assert 'Does not validate physical time' in result['limitation']
