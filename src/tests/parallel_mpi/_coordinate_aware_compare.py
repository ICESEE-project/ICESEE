# ==============================================================================
# @des: Reusable, physical-coordinate-aware comparator for two
# .npz dumps produced by _icepack_coordinate_noise_worker.py (or any
# other producer following the same {coords, increment, hdim,
# total_state_param_vars, vec_inputs} contract).
#
# Unlike np.sort(values), this keeps each value associated with its
# physical mesh coordinate: it builds a nearest-neighbor coordinate
# mapping from run B's nodes onto run A's nodes (via a KDTree, since both
# runs cover the identical physical mesh but Firedrake's own partition-
# dependent DOF ordering differs between them), reorders run B's values
# into run A's coordinate order, and only then compares value to value.
#
# Not part of any production code path: this is validation/test tooling,
# run offline, never imported by the DA driver.
# ==============================================================================
from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree


def check_coordinate_uniqueness(coords, atol=1e-9):
    """Return (is_unique, min_pairwise_distance) for a coordinate array.

    A duplicate/near-duplicate pair (distance below atol) would mean
    coordinates alone cannot key a DOF uniquely -- the caller would need
    to fall back to component/variable/global-DOF-id disambiguation
    instead of a pure nearest-neighbor coordinate match.
    """
    tree = cKDTree(coords)
    # k=2: nearest neighbor other than the point itself.
    dist, _ = tree.query(coords, k=2)
    min_dist = float(np.min(dist[:, 1]))
    return min_dist > atol, min_dist


def match_by_coordinate(coords_a, coords_b, atol=1e-6):
    """Return an index array `perm` such that coords_b[perm] ~= coords_a.

    Raises if any nearest-neighbor match exceeds `atol` (a genuine
    mismatch, not the expected exact-same-physical-mesh correspondence),
    or if the match is not a bijection (would indicate coordinates do not
    uniquely key each DOF -- see check_coordinate_uniqueness).
    """
    if coords_a.shape != coords_b.shape:
        raise ValueError(f"coordinate array shapes differ: {coords_a.shape} vs {coords_b.shape}")
    tree_b = cKDTree(coords_b)
    dist, idx = tree_b.query(coords_a, k=1)
    max_dist = float(np.max(dist))
    if max_dist > atol:
        raise ValueError(
            f"coordinate matching exceeded tolerance: max nearest-neighbor "
            f"distance {max_dist:.3e} > atol {atol:.3e} -- run A and run B "
            f"do not appear to cover the identical physical mesh"
        )
    if len(set(idx.tolist())) != len(idx):
        raise ValueError(
            "coordinate matching is not a bijection -- some coords_b node "
            "was matched by more than one coords_a node; coordinates alone "
            "do not uniquely key each DOF here"
        )
    return idx, max_dist


def compare_variable_blocks(npz_a_path, npz_b_path, coord_atol=1e-6, value_atol=1e-6, value_rtol=1e-6):
    """Coordinate-aware, per-variable comparison of two noise/state dumps.

    Returns a dict keyed by variable name, each entry containing:
      max_abs_diff, rms_diff, correlation, within_tolerance,
      sorted_only_would_match (True if a naive np.sort comparison would
      have reported equality despite genuine per-physical-DOF differences
      -- exactly the failure mode np.sort alone cannot detect).
    """
    a = np.load(npz_a_path, allow_pickle=True)
    b = np.load(npz_b_path, allow_pickle=True)

    coords_a = np.asarray(a["coords"])
    coords_b = np.asarray(b["coords"])
    hdim = int(a["hdim"])
    assert hdim == int(b["hdim"]), "hdim mismatch between the two runs"
    vec_inputs = [str(v) for v in a["vec_inputs"]]
    assert vec_inputs == [str(v) for v in b["vec_inputs"]]

    is_unique_a, min_dist_a = check_coordinate_uniqueness(coords_a)
    is_unique_b, min_dist_b = check_coordinate_uniqueness(coords_b)

    perm, coord_match_max_dist = match_by_coordinate(coords_a, coords_b, atol=coord_atol)

    inc_a = np.asarray(a["increment"])
    inc_b = np.asarray(b["increment"])

    results = {
        "_coordinate_uniqueness": {
            "run_a_unique": is_unique_a, "run_a_min_pairwise_dist": min_dist_a,
            "run_b_unique": is_unique_b, "run_b_min_pairwise_dist": min_dist_b,
            "coord_match_max_dist": coord_match_max_dist,
        }
    }
    for i, name in enumerate(vec_inputs):
        block_a = inc_a[i * hdim:(i + 1) * hdim]
        block_b_raw = inc_b[i * hdim:(i + 1) * hdim]
        block_b_reordered = block_b_raw[perm]

        diff = block_a - block_b_reordered
        max_abs = float(np.max(np.abs(diff)))
        rms = float(np.sqrt(np.mean(diff ** 2)))
        if np.std(block_a) > 0 and np.std(block_b_reordered) > 0:
            corr = float(np.corrcoef(block_a, block_b_reordered)[0, 1])
        else:
            corr = float("nan")
        within_tol = bool(
            np.allclose(block_a, block_b_reordered, atol=value_atol, rtol=value_rtol)
        )
        sorted_only_match = bool(
            np.allclose(np.sort(block_a), np.sort(block_b_raw), atol=value_atol, rtol=value_rtol)
        )
        results[name] = {
            "max_abs_diff": max_abs,
            "rms_diff": rms,
            "correlation": corr,
            "within_tolerance": within_tol,
            "sorted_only_would_match": sorted_only_match,
            "genuinely_different_field": sorted_only_match and not within_tol,
        }
    return results


if __name__ == "__main__":
    import sys
    import json

    out = compare_variable_blocks(sys.argv[1], sys.argv[2])
    print(json.dumps(out, indent=2, default=str))
