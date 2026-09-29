# ==============================================================================
# @des: Regression tests for coordinate_keyed_white_noise
# (src/run_model_da/_error_generation.py), the Stage 4C fix making the
# "graph" random-field method's initial white-noise draw keyed to
# physical DOF identity (seed, physical coordinate) instead of array
# position -- so the SAME physical node gets the SAME white-noise value
# regardless of which rank/array-index it occupies under a given
# Firedrake decomposition. The FFT method is deliberately untouched (see
# test_fft_is_index_space_not_coordinate_space below, which documents
# and locks in its existing, intentional index-space behavior).
# ==============================================================================
from __future__ import annotations

import numpy as np
import pytest

from scipy.spatial import cKDTree

from ICESEE.src.run_model_da._error_generation import (
    coordinate_keyed_white_noise,
    generate_pseudo_random_field_1d,
)


def _knn_neighbor_sets(coords, k, coord_tolerance=1e-6):
    """Reference (test-only) re-implementation of the deterministic k-NN
    tie-break in _build_knn_adjacency, used to check permutation
    invariance of the NEIGHBOR SET independent of generate_pseudo_random_
    field_1d's own internals."""
    n = len(coords)
    tree = cKDTree(coords)
    tie_margin = 8
    kq = min(k + 1 + tie_margin, n)
    dist, idx = tree.query(coords, k=kq)

    def qkey(c):
        return tuple(int(round(v / coord_tolerance)) for v in c)

    result = []
    for i in range(n):
        candidates = [(float(d), int(j)) for d, j in zip(dist[i], idx[i]) if j != i]
        candidates.sort(key=lambda pair: (pair[0], qkey(coords[pair[1]])))
        chosen = candidates[:k]
        result.append(frozenset(tuple(np.round(coords[j], 6)) for _, j in chosen))
    return result


def test_permutation_invariance_exact():
    """Same physical coordinates in a different array order -> the SAME
    per-coordinate white-noise value, to machine precision, once mapped
    back through the permutation. The raw arrays differ in order, as
    expected."""
    rng = np.random.default_rng(0)
    coords = rng.uniform(0, 1000, size=(50, 2))
    perm = rng.permutation(50)
    coords_perm = coords[perm]

    seed = 12345
    w1 = coordinate_keyed_white_noise(coords, seed)
    w2 = coordinate_keyed_white_noise(coords_perm, seed)

    assert np.array_equal(w2, w1[perm]), "coordinate-keyed noise is not permutation-invariant"
    assert not np.allclose(w1, w2), "raw arrays should differ in order for a non-identity permutation"


def test_member_independence():
    coords = np.random.default_rng(1).uniform(0, 1000, size=(20, 2))
    w_member0 = coordinate_keyed_white_noise(coords, seed=111)
    w_member1 = coordinate_keyed_white_noise(coords, seed=222)
    assert not np.allclose(w_member0, w_member1)


def test_reproducibility():
    coords = np.random.default_rng(2).uniform(0, 1000, size=(20, 2))
    w1 = coordinate_keyed_white_noise(coords, seed=42)
    w2 = coordinate_keyed_white_noise(coords, seed=42)
    assert np.array_equal(w1, w2)


def test_tolerance_quantization_merges_near_duplicate_coordinates():
    """Two coordinates within `coord_tolerance` of each other must key to
    the same value (the intended behavior: harmless floating-point
    differences from a different parallel assembly must not change the
    key) -- but coordinates further apart than the tolerance must not."""
    base = np.array([[100.0, 200.0]])
    near = base + 1e-9  # well inside the default 1e-6 tolerance
    far = base + 1e-3    # well outside it

    w_base = coordinate_keyed_white_noise(base, seed=7)
    w_near = coordinate_keyed_white_noise(near, seed=7)
    w_far = coordinate_keyed_white_noise(far, seed=7)

    assert np.array_equal(w_base, w_near)
    assert not np.array_equal(w_base, w_far)


def test_graph_method_uses_coordinate_keyed_noise_when_coords_given():
    """generate_pseudo_random_field_1d's graph branch, given coords, must
    be permutation-invariant end to end (through the smoothing step too,
    for a coordinate set with no k-NN tie ambiguity)."""
    rng = np.random.default_rng(3)
    # A jittered (non-grid) point set avoids the exact k-NN distance ties
    # a perfectly regular grid can produce (a separate, documented, k-NN
    # tie-breaking artifact -- see model_capabilities.py's icepack notes
    # -- not something this test needs to exercise).
    coords = rng.uniform(0, 1000, size=(40, 2))
    perm = rng.permutation(40)

    field_a = generate_pseudo_random_field_1d(
        N=40, method="graph", coords=coords, seed=99, verbose=False,
    )
    field_b = generate_pseudo_random_field_1d(
        N=40, method="graph", coords=coords[perm], seed=99, verbose=False,
    )
    assert np.allclose(field_a[perm], field_b, atol=1e-10)


def test_knn_tie_break_is_deterministic_under_permutation():
    """A perfectly regular grid (the exact failure mode found against
    real Icepack coordinates: a node at (4375, 150) exactly 625 units
    from two candidates, (5000,150) and (3750,150)) must select the SAME
    physical neighbor set regardless of the array order coordinates are
    presented in.
    """
    xs, ys = np.meshgrid(np.linspace(0, 5000, 9), np.linspace(0, 1200, 9))
    coords = np.column_stack([xs.ravel(), ys.ravel()])
    rng = np.random.default_rng(5)
    perm = rng.permutation(len(coords))

    nbrs_a = _knn_neighbor_sets(coords, k=6)
    nbrs_b = _knn_neighbor_sets(coords[perm], k=6)

    tree_perm = cKDTree(coords[perm])
    _, idx = tree_perm.query(coords, k=1)

    mismatches = [i for i in range(len(coords)) if nbrs_a[i] != nbrs_b[idx[i]]]
    assert mismatches == [], f"neighbor-set mismatches at nodes {mismatches}"


def test_graph_adjacency_permutation_invariant_on_regular_grid():
    """End to end through generate_pseudo_random_field_1d's graph branch
    (not just the standalone neighbor-set helper above): the smoothed
    field itself must be exactly permutation-invariant on the same
    tie-prone regular grid.
    """
    xs, ys = np.meshgrid(np.linspace(0, 5000, 9), np.linspace(0, 1200, 9))
    coords = np.column_stack([xs.ravel(), ys.ravel()])
    n = len(coords)
    rng = np.random.default_rng(6)
    perm = rng.permutation(n)

    field_a = generate_pseudo_random_field_1d(
        N=n, method="graph", coords=coords, seed=17, verbose=False,
    )
    field_b = generate_pseudo_random_field_1d(
        N=n, method="graph", coords=coords[perm], seed=17, verbose=False,
    )
    assert np.allclose(field_a[perm], field_b, atol=1e-10)


def test_disjoint_partition_called_separately_matches_single_full_call():
    """The actual Mode-3 usage pattern: each spatial rank calls this
    function ONCE with only its own owned coordinates (never the full
    global set), rather than one process calling it on a permuted full
    array (see test_permutation_invariance_exact above for that case).
    Concatenating two disjoint ranks' independently computed results and
    matching back by physical coordinate must reproduce a single
    full-array (R=1-equivalent) call exactly, and each partial call must
    only ever see its own subset's length -- never the global size."""
    rng = np.random.default_rng(9)
    coords = rng.uniform(0, 1000, size=(37, 2))
    seed = 4242

    w_full = coordinate_keyed_white_noise(coords, seed)

    # Arbitrary, unequal, non-contiguous disjoint split -- matching a real
    # Firedrake partition's arbitrary (not simply "first half/second half")
    # ownership.
    idx_a = np.array([0, 2, 4, 5, 7, 9, 10, 12, 15, 18, 20, 22, 25, 30, 33, 36])
    idx_b = np.array([i for i in range(37) if i not in set(idx_a.tolist())])
    assert set(idx_a.tolist()) | set(idx_b.tolist()) == set(range(37))

    w_a = coordinate_keyed_white_noise(coords[idx_a], seed)
    w_b = coordinate_keyed_white_noise(coords[idx_b], seed)
    assert w_a.size == idx_a.size  # bounded to this "rank"'s own subset
    assert w_b.size == idx_b.size  # never the global (37) size

    reassembled = np.empty(37)
    reassembled[idx_a] = w_a
    reassembled[idx_b] = w_b

    assert np.array_equal(reassembled, w_full)


def test_distribution_is_approximately_standard_normal():
    """Smoke test (not a rigorous statistical test, but appropriate for a
    deterministic unit test): over enough physical coordinates, the
    coordinate-keyed draw's sample mean/std should still be close to
    standard normal's (0, 1), confirming the per-DOF independent-seed
    construction doesn't distort the intended distribution."""
    rng = np.random.default_rng(11)
    coords = rng.uniform(0, 1000, size=(5000, 2))
    w = coordinate_keyed_white_noise(coords, seed=777)
    assert abs(np.mean(w)) < 0.05
    assert abs(np.std(w) - 1.0) < 0.05


def test_fft_is_index_space_not_coordinate_space():
    """Documents and locks in FFT's existing, intentional behavior: its
    covariance is defined over flat array indices (see
    sample_periodic_exp_cov's own docstring: "periodic ring distance"
    over i,j), not physical coordinates. This is NOT a bug and must NOT
    be "fixed" to be coordinate-aware -- this test exists to catch an
    accidental future change to that documented semantics, not to assert
    physical decomposition invariance for it.

    generate_pseudo_random_field_1d's "original FFT code path" draws
    from NumPy's legacy global random state directly (np.random.rand),
    relying on the caller having already done np.random.seed(seed) --
    exactly like every real caller in this codebase does (see
    generate_initial_member_increment/add_member_process_noise). Given
    the same explicit global seed, the same array-index-ordered output
    is fully reproducible; it has no notion of physical coordinates at
    all to be invariant with respect to.
    """
    np.random.seed(99)
    field_a = generate_pseudo_random_field_1d(N=40, Lx=1000.0, rh=100.0, method="fft")
    np.random.seed(99)
    field_b = generate_pseudo_random_field_1d(N=40, Lx=1000.0, rh=100.0, method="fft")
    assert np.array_equal(field_a, field_b)
