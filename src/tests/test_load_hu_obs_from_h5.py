# ==============================================================================
# @des: Tests for load_hu_obs_from_h5 (src/utils/tools.py), ported from
# test/ICESEE during Pass 2 reconciliation. It transparently reconstructs
# the dense (state_dim x m_obs) synthetic-observation matrix from either
# storage layout ICESEE writes to synthetic_obs.h5:
#   - dense: a full "hu_obs" dataset (modes 0/1 always; mode 2 when
#     synthetic_observation_storage is 'dense' or 'both').
#   - compact: "hu_obs_compact" (nobs, m_obs) + "obs_indices" (nobs,) sparse
#     rows + a "state_dimension" attr, scattered into a zero-filled dense
#     array (mode 2's default, for bounded per-rank memory).
# Small in-memory HDF5 fixtures (h5py.File(..., driver="core", backing_store
# =False)) stand in for real scientific datasets.
# ==============================================================================
import h5py
import numpy as np
import pytest

from ICESEE.src.utils.tools import load_hu_obs_from_h5


def _memory_h5():
    return h5py.File("in_memory_test.h5", "w", driver="core", backing_store=False)


def test_dense_representation_returns_full_matrix():
    nd, m_obs = 6, 3
    expected = np.arange(nd * m_obs, dtype=np.float64).reshape(nd, m_obs)
    with _memory_h5() as f:
        f.create_dataset("hu_obs", data=expected)
        result = load_hu_obs_from_h5(f)
    np.testing.assert_array_equal(result, expected)
    assert result.shape == (nd, m_obs)
    assert result.dtype == np.float64


def test_compact_representation_scatters_into_dense_zero_filled_array():
    nd, m_obs = 6, 2
    obs_indices = np.array([1, 3, 4], dtype=np.int64)
    compact = np.array([[10.0, 11.0], [20.0, 21.0], [30.0, 31.0]])
    with _memory_h5() as f:
        f.create_dataset("hu_obs_compact", data=compact)
        f.create_dataset("obs_indices", data=obs_indices)
        f.attrs["state_dimension"] = nd
        result = load_hu_obs_from_h5(f)

    assert result.shape == (nd, m_obs)
    expected = np.zeros((nd, m_obs))
    expected[obs_indices, :] = compact
    np.testing.assert_array_equal(result, expected)
    # Unobserved rows stay exactly zero.
    for row in range(nd):
        if row not in obs_indices:
            np.testing.assert_array_equal(result[row], np.zeros(m_obs))


def test_compact_representation_without_state_dimension_attr_infers_from_max_index():
    obs_indices = np.array([0, 2], dtype=np.int64)
    compact = np.array([[1.0], [2.0]])
    with _memory_h5() as f:
        f.create_dataset("hu_obs_compact", data=compact)
        f.create_dataset("obs_indices", data=obs_indices)
        # no state_dimension attr set
        result = load_hu_obs_from_h5(f)
    # Inferred nd = max(obs_indices) + 1 = 3
    assert result.shape == (3, 1)
    np.testing.assert_array_equal(result[[0, 2], :], compact)
    np.testing.assert_array_equal(result[1], [0.0])


def test_dense_takes_precedence_when_both_present():
    # ICESEE writes both when synthetic_observation_storage == 'both'; the
    # dense dataset is authoritative and cheaper to read directly.
    dense = np.ones((4, 2))
    with _memory_h5() as f:
        f.create_dataset("hu_obs", data=dense)
        f.create_dataset("hu_obs_compact", data=np.zeros((1, 2)))
        f.create_dataset("obs_indices", data=np.array([0], dtype=np.int64))
        result = load_hu_obs_from_h5(f)
    np.testing.assert_array_equal(result, dense)


def test_missing_both_representations_raises_key_error():
    with _memory_h5() as f:
        f.create_dataset("unrelated_dataset", data=np.zeros(3))
        with pytest.raises(KeyError):
            load_hu_obs_from_h5(f)


def test_missing_obs_indices_with_only_compact_raises_key_error():
    with _memory_h5() as f:
        f.create_dataset("hu_obs_compact", data=np.zeros((2, 2)))
        # obs_indices intentionally omitted
        with pytest.raises(KeyError):
            load_hu_obs_from_h5(f)
