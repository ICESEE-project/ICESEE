"""Serial-reference to distributed-layout node numbering (execution mode 3).

PACE job 14548042 (two members, model_nprocs=2) and job 14053366 (twenty
members, model_nprocs=2) stored member rows whose values are a permutation
of the serial truth rows: elementwise |member - truth| reached 2.9e3 m in
thickness while the sorted values agreed, and the basal-melt field matched
exactly as a multiset but differed by 20 m/yr elementwise. Observations
generated on the serial mesh were therefore assimilated at the wrong nodes.
"""

import numpy as np
import pytest

from ICESEE.applications.icepack_model.icepack_utils._reference_numbering import (
    match_owned_nodes_to_reference,
    reference_to_layout_permutation,
    translate_observation_ids,
)


class _GroupComm:
    """allgather over a simulated spatial communicator of several ranks."""

    def __init__(self, pieces):
        self.pieces = pieces

    def allgather(self, value):
        return list(self.pieces)


def _reference_coords(n=7, seed=0):
    return np.random.default_rng(seed).uniform(-1e6, 1e6, size=(n, 2))


def test_two_rank_renumbered_layout_maps_every_observation_to_its_node():
    reference = _reference_coords()
    # Distributed numbering: rank 0 owns rows 0..2, rank 1 rows 3..6, each a
    # different physical node than the serial row with the same index.
    layout_order = np.array([5, 0, 3, 6, 1, 4, 2])
    owned = [(0, reference[layout_order[:3]]), (3, reference[layout_order[3:]])]
    pieces = [
        (start, match_owned_nodes_to_reference(reference, coords))
        for start, coords in owned
    ]
    perm = reference_to_layout_permutation(
        pieces[0][0], pieces[0][1], _GroupComm(pieces), reference.shape[0]
    )
    # perm[reference node] is the distributed row holding that node.
    np.testing.assert_array_equal(layout_order[perm], np.arange(7))

    n = reference.shape[0]
    ids = np.array([3 * n + 0, 3 * n + 6, 0 * n + 2])  # two s rows, one h row
    layout_ids = translate_observation_ids(ids, perm)
    np.testing.assert_array_equal(layout_ids // n, ids // n)
    np.testing.assert_array_equal(layout_order[layout_ids % n], ids % n)


def test_identical_numbering_is_the_identity():
    reference = _reference_coords()
    index = match_owned_nodes_to_reference(reference, reference)
    perm = reference_to_layout_permutation(
        0, index, _GroupComm([(0, index)]), reference.shape[0]
    )
    np.testing.assert_array_equal(perm, np.arange(reference.shape[0]))


def test_unmatched_node_is_refused():
    reference = _reference_coords()
    owned = reference[:3].copy()
    owned[1] += 10.0
    with pytest.raises(ValueError, match="no matching reference node"):
        match_owned_nodes_to_reference(reference, owned)


def test_incomplete_or_overlapping_ownership_is_refused():
    reference = _reference_coords()
    first = match_owned_nodes_to_reference(reference, reference[:3])
    with pytest.raises(ValueError, match="not owned"):
        reference_to_layout_permutation(0, first, _GroupComm([(0, first)]), 7)
    overlap = [(0, first), (3, match_owned_nodes_to_reference(reference, reference[2:]))]
    with pytest.raises(ValueError, match="more than one"):
        reference_to_layout_permutation(0, first, _GroupComm(overlap), 7)
