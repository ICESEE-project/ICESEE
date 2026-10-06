# ==============================================================================
# @des: Regression test for the missing `import h5py` in
# applications/icepack_model/examples/synthetic_ice_stream/_icepack_enkf.py.
#
# The file references h5py.Dataset (in generate_true_state's
# isinstance(statevec_true, h5py.Dataset) check, used to skip returning the
# full in-memory trajectory when the state is already file-backed) without
# importing h5py itself -- a NameError waiting to happen the first time that
# branch is actually reached. This test exercises that exact branch with a
# real h5py.Dataset, using lightweight duck-typed stand-ins for the
# Firedrake Function objects (h0/u0) so it does not require a real Icepack
# solver/mesh, and setting nt=0 to skip the simulation loop entirely --
# only the pre-loop initial-state write and the isinstance/flush/return-None
# path are exercised.
# ==============================================================================
import os
import sys
import tempfile

import h5py
import numpy as np

# _icepack_model.py (imported by _icepack_enkf.py) imports
# config/_utility_imports.py, which parses sys.argv, loads a params.yaml
# file, and deletes/recreates icesee_kwargs['data_path'] at module import
# time -- see test_enkf_serial_process_noise.py for the same pattern. Point
# -F at Lorenz96's example config purely to satisfy the import-time read
# (nothing in this test depends on its contents), and --data_path at a
# disposable temp directory so the auto-clean guard never touches the
# repository's own _modelrun_datasets.
_lorenz96_params = os.path.join(
    os.path.dirname(__file__),
    "..", "..", "applications", "lorenz_model", "examples", "lorenz96", "params.yaml",
)
_safe_data_path = tempfile.mkdtemp(prefix="icesee_synthetic_ice_stream_test_")
_saved_argv = sys.argv
sys.argv = [sys.argv[0], "-F", _lorenz96_params, "--data_path", _safe_data_path]
try:
    from ICESEE.applications.icepack_model.examples.synthetic_ice_stream._icepack_enkf import (
        generate_true_state,
    )
finally:
    sys.argv = _saved_argv


class _FakeDat:
    def __init__(self, data):
        self.data_ro = data


class _FakeFiredrakeFunction:
    """Duck-typed stand-in for a Firedrake Function: only .dat.data_ro and
    .copy(deepcopy=True) are used by generate_true_state before its
    simulation loop, which this test disables via nt=0."""

    def __init__(self, data):
        self.dat = _FakeDat(np.asarray(data))

    def copy(self, deepcopy=True):
        return _FakeFiredrakeFunction(np.array(self.dat.data_ro, copy=True))


def test_generate_true_state_returns_none_for_file_backed_h5py_dataset(tmp_path):
    hdim = 3
    vec_inputs = ["h", "u", "v"]
    nd = hdim * len(vec_inputs)

    h0 = _FakeFiredrakeFunction(np.array([1.0, 2.0, 3.0]))
    u0 = _FakeFiredrakeFunction(np.array([[4.0, 7.0], [5.0, 8.0], [6.0, 9.0]]))

    h5_path = tmp_path / "true_state.h5"
    with h5py.File(h5_path, "w") as f:
        statevec_true = f.create_dataset("true_state", (nd, 1), dtype="f8")

        result = generate_true_state(
            statevec_true=statevec_true,
            h0=h0,
            u0=u0,
            a=None,
            b=None,
            dt=1.0,
            A=None,
            C=None,
            Q=None,
            V=None,
            solver=None,
            nt=0,
            joint_estimation=False,
            vec_inputs=vec_inputs,
            nd=nd,
            default_run=True,
            even_distribution=False,
        )

        # The isinstance(statevec_true, h5py.Dataset) branch must be taken
        # (this is the exact line that used to NameError without the
        # missing `import h5py`) and must return None rather than trying
        # to read the whole file-backed trajectory back into memory.
        assert result is None

        # The initial state was still written through to the dataset.
        np.testing.assert_array_equal(statevec_true[0:3, 0], [1.0, 2.0, 3.0])
        np.testing.assert_array_equal(statevec_true[3:6, 0], [4.0, 5.0, 6.0])
        np.testing.assert_array_equal(statevec_true[6:9, 0], [7.0, 8.0, 9.0])
