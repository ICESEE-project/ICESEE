# ==============================================================================
# @des: Level-1 unit tests for src/utils/state_ownership.py.
#
# These test the ownership-resolution arithmetic and the scalar/array
# normalization helper directly, without launching real MPI ranks. A
# lightweight in-process fake communicator stands in for a multi-rank
# ``subcomm`` where the arithmetic needs to be exercised across several
# simulated ranks (Level 2 concerns -- real collectives across real
# processes -- belong in the MPI-launched tests instead).
# ==============================================================================
import numpy as np
import pytest
from mpi4py import MPI

from ICESEE.src.utils.state_ownership import (
    StateOwnership,
    resolve_state_ownership,
    combine_member_state,
    ensure_state_array,
)


class _FakeSubcomm:
    """Minimal stand-in for one rank's view of a multi-rank subcomm.

    Simulates ``allgather``/``gather`` by returning the pre-supplied
    per-rank values directly, so ownership arithmetic can be exercised
    for several simulated ranks without any real MPI communication.
    """

    def __init__(self, rank, per_rank_values):
        self._rank = rank
        self._per_rank_values = per_rank_values

    def Get_rank(self):
        return self._rank

    def Get_size(self):
        return len(self._per_rank_values)

    def allgather(self, value):
        return list(self._per_rank_values)

    def gather(self, value, root=0):
        # Real mpi4py only returns the assembled list on `root`; every
        # other rank gets None. Emulate that so a test cannot accidentally
        # depend on non-root return values, matching real behavior.
        if self._rank != root:
            return None
        return list(self._per_rank_values)


# ---------------------------------------------------------------------------
# resolve_state_ownership
# ---------------------------------------------------------------------------

def test_replicated_state_ownership_holds_regardless_of_subcomm_size():
    """global_size == local_size for replicated state, for every subcomm size."""
    for sub_size in (1, 2, 5):
        subcomm = _FakeSubcomm(rank=0, per_rank_values=[7] * sub_size)
        icesee_kwargs = {"nd": 7, "state_distribution": "replicated"}
        ownership = resolve_state_ownership(icesee_kwargs, subcomm)
        assert ownership.distribution == "replicated"
        assert ownership.global_size == 7
        assert ownership.local_size == 7
        assert ownership.local_offset == 0


def test_default_distribution_is_distributed_not_replicated():
    """The default must preserve pre-existing (sum/concatenate) arithmetic."""
    subcomm = _FakeSubcomm(rank=0, per_rank_values=[3, 3])
    icesee_kwargs = {"nd": 3}  # no state_distribution key at all
    ownership = resolve_state_ownership(icesee_kwargs, subcomm)
    assert ownership.distribution == "distributed"
    assert ownership.global_size == 6


def test_distributed_ownership_sums_disjoint_local_sizes():
    """global_size == sum(local_size_r) for disjoint partitions."""
    local_sizes = [4, 6, 2, 8]
    for rank, expected_offset in enumerate([0, 4, 10, 12]):
        subcomm = _FakeSubcomm(rank=rank, per_rank_values=local_sizes)
        icesee_kwargs = {"nd": local_sizes[rank], "state_distribution": "distributed"}
        ownership = resolve_state_ownership(icesee_kwargs, subcomm)
        assert ownership.distribution == "distributed"
        assert ownership.global_size == sum(local_sizes)
        assert ownership.local_size == local_sizes[rank]
        assert ownership.local_offset == expected_offset


def test_single_rank_subcomm_is_replicated_regardless_of_config():
    """No ownership ambiguity exists at subcomm size 1 (today's only tested
    configuration for every shipped application): local == global either
    way, so this must not require an allgather at all."""
    subcomm = _FakeSubcomm(rank=0, per_rank_values=[9])
    for requested in ("replicated", "distributed"):
        icesee_kwargs = {"nd": 9, "state_distribution": requested}
        ownership = resolve_state_ownership(icesee_kwargs, subcomm)
        assert ownership.distribution == "replicated"
        assert ownership.global_size == 9
        assert ownership.local_size == 9


def test_invalid_distribution_value_raises():
    subcomm = _FakeSubcomm(rank=0, per_rank_values=[1, 1])
    icesee_kwargs = {"nd": 1, "state_distribution": "sharded"}
    with pytest.raises(ValueError):
        resolve_state_ownership(icesee_kwargs, subcomm)


def test_state_ownership_rejects_invalid_distribution_directly():
    with pytest.raises(ValueError):
        StateOwnership(distribution="sharded", global_size=1, local_size=1)


# ---------------------------------------------------------------------------
# combine_member_state
# ---------------------------------------------------------------------------

def test_combine_replicated_state_is_a_no_op():
    """No communication for replicated state: data returned unchanged."""
    ownership = StateOwnership(distribution="replicated", global_size=3, local_size=3)
    data = {"x": np.array([1.0, 2.0, 3.0])}
    subcomm = _FakeSubcomm(rank=0, per_rank_values=[data])
    result = combine_member_state(ownership, subcomm, data, root=0)
    assert result is data  # identity, not merely equal -- confirms no copy/gather happened


def test_combine_distributed_state_concatenates_dict_pieces():
    piece_rank0 = {"x": np.array([1.0, 2.0])}
    piece_rank1 = {"x": np.array([3.0, 4.0, 5.0])}
    ownership = StateOwnership(distribution="distributed", global_size=5, local_size=2)
    subcomm = _FakeSubcomm(rank=0, per_rank_values=[piece_rank0["x"], piece_rank1["x"]])
    result = combine_member_state(ownership, subcomm, piece_rank0, root=0)
    np.testing.assert_array_equal(result["x"], np.array([1.0, 2.0, 3.0, 4.0, 5.0]))


def test_combine_distributed_state_concatenates_array_pieces_not_stacks():
    """np.concatenate, not np.vstack: 1-D pieces must join end-to-end into
    one longer 1-D vector, not be promoted into rows of a new 2-D array."""
    ownership = StateOwnership(distribution="distributed", global_size=4, local_size=2)
    subcomm = _FakeSubcomm(
        rank=0, per_rank_values=[np.array([1.0, 2.0]), np.array([3.0, 4.0])]
    )
    result = combine_member_state(ownership, subcomm, np.array([1.0, 2.0]), root=0)
    assert result.shape == (4,)
    np.testing.assert_array_equal(result, np.array([1.0, 2.0, 3.0, 4.0]))


def test_combine_distributed_state_returns_none_off_root():
    ownership = StateOwnership(distribution="distributed", global_size=4, local_size=2)
    subcomm = _FakeSubcomm(
        rank=1, per_rank_values=[np.array([1.0, 2.0]), np.array([3.0, 4.0])]
    )
    result = combine_member_state(ownership, subcomm, np.array([3.0, 4.0]), root=0)
    assert result is None


# ---------------------------------------------------------------------------
# ensure_state_array
# ---------------------------------------------------------------------------

def test_ensure_state_array_promotes_bare_scalar():
    """The confirmed Stage-3 bug: initialize_ensemble returning u0b[0] (a
    numpy.float64 scalar) must become a proper one-element array."""
    scalar = np.array([2.0, 3.0, 4.0])[0]  # numpy.float64, matches u0b[0]
    assert np.isscalar(scalar) or scalar.ndim == 0
    result = ensure_state_array(scalar)
    assert result.shape == (1,)
    assert result[0] == 2.0


def test_ensure_state_array_promotes_python_float():
    result = ensure_state_array(5.0)
    assert result.shape == (1,)
    assert result[0] == 5.0


def test_ensure_state_array_promotes_0d_array():
    result = ensure_state_array(np.array(6.0))
    assert result.shape == (1,)
    assert result[0] == 6.0


def test_ensure_state_array_leaves_1d_array_unchanged():
    arr = np.array([1.0, 2.0, 3.0])
    result = ensure_state_array(arr)
    assert result.shape == (3,)
    np.testing.assert_array_equal(result, arr)


def test_ensure_state_array_leaves_multidimensional_array_unchanged():
    arr = np.zeros((4, 7))
    result = ensure_state_array(arr)
    assert result.shape == (4, 7)


def test_ensure_state_array_distinguishes_scalar_from_one_element_vector():
    """A scalar state component and a genuine one-element state vector must
    normalize to the same shape -- they are the same logical quantity."""
    from_scalar = ensure_state_array(np.float64(9.0))
    from_one_element_vector = ensure_state_array(np.array([9.0]))
    assert from_scalar.shape == from_one_element_vector.shape == (1,)
    np.testing.assert_array_equal(from_scalar, from_one_element_vector)
