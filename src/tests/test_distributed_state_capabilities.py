# ==============================================================================
# @des: Level-1 unit tests for
# src/parallelization/distributed_state_capabilities.py (Stage 4D.3).
# Pure Python -- no MPI, no Firedrake, no real HDF5 file required.
# ==============================================================================
from ICESEE.src.parallelization.distributed_state_capabilities import (
    StateBackendCapabilities,
    detect_state_backend_capabilities,
)


class _FakeComm:
    def __init__(self, size):
        self._size = size

    def Get_size(self):
        return self._size


def test_detects_real_local_h5py_mpi_flag():
    # This is a real, unmocked probe -- not asserting a specific value
    # (that depends on how h5py was built on whatever machine runs this
    # test), only that it returns a bool and does not raise.
    caps = detect_state_backend_capabilities()
    assert isinstance(caps.h5py_mpi_enabled, bool)


def test_no_comm_argument_falls_back_to_mpi_comm_world_if_importable():
    caps = detect_state_backend_capabilities()
    assert isinstance(caps.mpi_available, bool)
    assert caps.mpi_world_size >= 1


def test_explicit_comm_overrides_world_size():
    caps = detect_state_backend_capabilities(comm=_FakeComm(8))
    assert caps.mpi_available is True
    assert caps.mpi_world_size == 8


def test_explicit_comm_that_raises_reports_unavailable_not_a_crash():
    class _BrokenComm:
        def Get_size(self):
            raise RuntimeError("no communicator")

    caps = detect_state_backend_capabilities(comm=_BrokenComm())
    assert caps.mpi_available is False
    assert caps.mpi_world_size == 1


def test_supports_collective_hdf5_requires_both_mpi_h5py_and_multiple_ranks():
    single_rank = StateBackendCapabilities(
        h5py_mpi_enabled=True, mpi_available=True, mpi_world_size=1,
        node_local_scratch=None,
    )
    assert single_rank.supports_collective_hdf5() is False

    serial_h5py = StateBackendCapabilities(
        h5py_mpi_enabled=False, mpi_available=True, mpi_world_size=4,
        node_local_scratch=None,
    )
    assert serial_h5py.supports_collective_hdf5() is False

    no_mpi = StateBackendCapabilities(
        h5py_mpi_enabled=True, mpi_available=False, mpi_world_size=1,
        node_local_scratch=None,
    )
    assert no_mpi.supports_collective_hdf5() is False

    genuinely_ready = StateBackendCapabilities(
        h5py_mpi_enabled=True, mpi_available=True, mpi_world_size=4,
        node_local_scratch=None,
    )
    assert genuinely_ready.supports_collective_hdf5() is True


def test_node_local_scratch_only_reported_if_env_var_points_at_a_real_directory(
    monkeypatch, tmp_path
):
    monkeypatch.delenv("ICESEE_NODE_LOCAL_SCRATCH", raising=False)
    monkeypatch.delenv("TMPDIR", raising=False)
    monkeypatch.delenv("SLURM_TMPDIR", raising=False)
    caps = detect_state_backend_capabilities()
    assert caps.node_local_scratch is None

    monkeypatch.setenv("ICESEE_NODE_LOCAL_SCRATCH", str(tmp_path))
    caps = detect_state_backend_capabilities()
    assert caps.node_local_scratch == str(tmp_path)


def test_node_local_scratch_env_var_pointing_nowhere_is_ignored(monkeypatch):
    monkeypatch.setenv("ICESEE_NODE_LOCAL_SCRATCH", "/no/such/directory/xyz")
    monkeypatch.delenv("TMPDIR", raising=False)
    monkeypatch.delenv("SLURM_TMPDIR", raising=False)
    caps = detect_state_backend_capabilities()
    assert caps.node_local_scratch is None
