#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# =============================================================================
# @author: Brian Kyanjo
# @date: December 2025
# @description:
#   Generate a kriging-based conditional random bed ensemble for ICESEE/ISSM.
#
#   Strategy:
#     bed_ensemble = bed_background + kriged(residuals) + conditional_noise
#
#   This preserves the large-scale bed structure while perturbing only
#   unresolved components.
#
#   Output:
#       bed_ens shape = (npts, Ne)
# =============================================================================

import argparse
import os
import sys
from typing import Optional

import h5py
import numpy as np
import gstools as gs
from gstools import krige
from scipy.spatial import cKDTree

try:
    from tqdm import trange
    HAVE_TQDM = True
except ImportError:
    HAVE_TQDM = False


# ---------------------------------------------------------------------
# Utility
# ---------------------------------------------------------------------
def _rmse(a, b):
    return float(np.sqrt(np.mean((np.asarray(a) - np.asarray(b)) ** 2)))


def smooth_background_2d(
    field,
    x,
    y,
    corr_len_m=40_000.0,
    truncation_sigma=3.0,
):
    """Return a coordinate-aware smooth bed background on an unstructured mesh.

    The previous one-dimensional convolution followed ISSM vertex numbering,
    which is not a spatial coordinate and can mix unrelated parts of the bed.
    This Gaussian-distance average is invariant to mesh-node ordering and uses
    both horizontal coordinates.
    """
    values = np.asarray(field, dtype=float).ravel()
    x = np.asarray(x, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()
    if not (values.size == x.size == y.size):
        raise ValueError("field, x, and y must have the same length")
    if values.size == 0:
        return values.copy()

    corr_len_m = float(corr_len_m)
    truncation_sigma = float(truncation_sigma)
    if corr_len_m <= 0.0:
        raise ValueError("corr_len_m must be positive")
    if truncation_sigma <= 0.0:
        raise ValueError("truncation_sigma must be positive")
    if not np.all(np.isfinite(values)) or not np.all(np.isfinite(x + y)):
        raise ValueError("field and coordinates must be finite")

    coords = np.column_stack((x, y))
    tree = cKDTree(coords)
    radius = truncation_sigma * corr_len_m
    background = np.empty_like(values)
    for row, neighbors in enumerate(tree.query_ball_point(coords, r=radius)):
        index = np.asarray(neighbors, dtype=int)
        # KD-tree result order is an implementation detail. Coordinate sorting
        # makes floating-point summation deterministic under vertex permutation.
        order = np.lexsort((y[index], x[index]))
        index = index[order]
        distance = np.hypot(x[index] - x[row], y[index] - y[row])
        weights = np.exp(-0.5 * (distance / corr_len_m) ** 2)
        background[row] = np.dot(weights, values[index]) / np.sum(weights)
    return background


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
def generate_bed_kriging(
    icesee_path: str = "./",
    data_path: str = "_modelrun_datasets",
    Ne: int = 60,
    stride_km: float = 10.0,
    track_half_width_m: float = 2_500.0,
    grounded_observations_only: bool = True,
    snap_idx: int = 0,
    sigma_noise: float = 10.0,
    sill_bed: float = 5000.0,
    range_bed: float = 40000.0,
    nugget_bed: float = 50.0,
    ensemble_inflation: float = 1.5,
    background_length_m: float = 40_000.0,
    background_truncation_sigma: float = 3.0,
    enforce_obs: bool = False,
    use_floating_mask: bool = False,
    seed_base: int = 1234,
    output_file: Optional[str] = None,
):

    base_dir = os.path.join(icesee_path, data_path)
    mesh_file = os.path.join(base_dir, "mesh_idxy_0.h5")
    true_file = os.path.join(base_dir, "true_nurged_states.h5")

    if output_file is None:
        output_file = os.path.join(base_dir, "bed_kriging_results.h5")

    print("[bed-kriging] ================================================")
    print("[bed-kriging] Residual kriging mode")
    print("[bed-kriging] ================================================")

    # ------------------------------------------------------------------
    # Read mesh
    # ------------------------------------------------------------------
    with h5py.File(mesh_file, "r") as f:
        x = f["/fric_x"][:].ravel()
        y = f["/fric_y"][:].ravel()

    npts = x.size

    # ------------------------------------------------------------------
    # Read true state
    # ------------------------------------------------------------------
    with h5py.File(true_file, "r") as f:
        true_state = f["true_state"][:]

    ndim = true_state.shape[0] // 6
    nt = true_state.shape[1]
    
    H_true   = true_state[0:ndim, snap_idx]
    bed_true = true_state[4 * ndim:5 * ndim, snap_idx]

    # ------------------------------------------------------------------
    # Large-scale background
    # ------------------------------------------------------------------
    bed_bg = smooth_background_2d(
        bed_true,
        x,
        y,
        corr_len_m=background_length_m,
        truncation_sigma=background_truncation_sigma,
    )

    # residual field
    bed_resid = bed_true - bed_bg

    # ------------------------------------------------------------------
    # Radar-track mask
    # ------------------------------------------------------------------
    x_km = x / 1000.0
    track_xs = np.arange(x_km.min(), x_km.max(), stride_km)

    obs_mask = np.zeros(npts, dtype=bool)

    band = float(track_half_width_m) / 1000.0
    if stride_km <= 0.0 or band < 0.0:
        raise ValueError("track stride must be positive and width non-negative")
    for tx in track_xs:
        obs_mask |= np.abs(x_km - tx) <= band

    # Match the EnKF bed-observation experiment: radar returns constrain only
    # the true grounded region, not the floating shelf bed.
    density_ratio = 917.0 / 1028.0
    grounded_mask = H_true + bed_true / density_ratio > 0.0
    if grounded_observations_only:
        obs_mask &= grounded_mask

    # x_min, x_max = float(x_km.min()), float(x_km.max())

    # track_xs = np.arange(x_min, x_max + 1.0e-12, float(stride_km))

    # if track_xs.size > 1:
    #     dx_nom = (x_max - x_min) / (track_xs.size - 1)
    # else:
    #     dx_nom = float(stride_km)

    # band = min(0.25, 0.15 * dx_nom)

    # # Hydrostatic grounding criterion
    # di = 917.0 / 1028.0
    # ocean_levelset = H_true + bed_true / di

    # grounded_mask = ocean_levelset > 0

    # # Near-grounding-line band
    # near_gl_mask = np.abs(ocean_levelset) < 50.0   # meters

    # dynamic_mask = grounded_mask | near_gl_mask

    # obs_mask = np.zeros(npts, dtype=bool)

    # for x_line in track_xs:
    #     track_mask = np.abs(x_km - x_line) <= band
    #     obs_mask |= track_mask & dynamic_mask

    n_obs = obs_mask.sum()

    print(f"[bed-kriging] Radar-bed obs points: {n_obs} / {npts}")

    # ------------------------------------------------------------------
    # Synthetic observations on residuals
    # ------------------------------------------------------------------
    rng = np.random.default_rng(seed_base)

    x_obs = x[obs_mask]
    y_obs = y[obs_mask]

    z_obs = (
        bed_resid[obs_mask]
        + sigma_noise * rng.standard_normal(n_obs)
    )

    # ------------------------------------------------------------------
    # Optional floating-prior mask
    # ------------------------------------------------------------------
    floating_mask = np.zeros(npts, dtype=bool)

    if use_floating_mask:
        # Legacy runs used the hidden true bed both to define this mask and
        # to overwrite every masked ensemble value.  That made the floating
        # bed artificially perfect and imprinted a truth-shaped edge in the
        # prior.  If this optional guard is requested, diagnose it from the
        # smooth background and retain that background instead.
        floating_mask = bed_bg < -900.0
        print(
            "[bed-kriging] Deep-background points held at the smooth prior: "
            f"{floating_mask.sum()}"
        )

    # ------------------------------------------------------------------
    # Residual kriging
    # ------------------------------------------------------------------
    model = gs.Exponential(
        dim=2,
        var=sill_bed,
        len_scale=range_bed,
        nugget=nugget_bed,
    )

    ok = krige.Ordinary(
        model,
        cond_pos=[x_obs, y_obs],
        cond_val=z_obs,
    )

    resid_kriged, resid_var = ok([x, y], return_var=True)

    resid_kriged = resid_kriged.ravel()
    resid_var = resid_var.ravel()

    # reconstructed mean bed
    bed_kriged = bed_bg + resid_kriged
    if use_floating_mask:
        bed_kriged[floating_mask] = bed_bg[floating_mask]
        resid_var[floating_mask] = 0.0

    # ------------------------------------------------------------------
    # Conditional residual ensemble
    # ------------------------------------------------------------------
    cond_srf = gs.CondSRF(ok)

    bed_ens = np.zeros((npts, Ne))

    iterator = trange(Ne, desc="[bed-kriging] ensemble") if HAVE_TQDM else range(Ne)

    for e in iterator:
        resid_field = cond_srf([x, y], seed=seed_base + e)
        resid_field = np.asarray(resid_field).ravel()

        resid_field = (
            resid_kriged
            + ensemble_inflation * (resid_field - resid_kriged)
        )

        field = bed_bg + resid_field

        if enforce_obs:
            field[obs_mask] = bed_true[obs_mask]

        if use_floating_mask:
            field[floating_mask] = bed_bg[floating_mask]

        bed_ens[:, e] = field

    # ------------------------------------------------------------------
    # Diagnostics
    # ------------------------------------------------------------------
    bed_mean = np.mean(bed_ens, axis=1)
    bed_spread = np.std(bed_ens, axis=1)

    diagnostics = {
        "rmse_kriged_all": _rmse(bed_kriged, bed_true),
        "rmse_ensmean_all": _rmse(bed_mean, bed_true),
        "spread_all": float(np.mean(bed_spread)),
        "n_obs": int(n_obs),
        "grounded_coverage_percent": float(
            100.0 * np.count_nonzero(obs_mask & grounded_mask)
            / max(np.count_nonzero(grounded_mask), 1)
        ),
    }

    print("[bed-kriging] Diagnostics:")
    for k, v in diagnostics.items():
        print(f"   {k}: {v}")

    # ------------------------------------------------------------------
    # Write output
    # ------------------------------------------------------------------
    with h5py.File(output_file, "w") as f:
        f["x"] = x
        f["y"] = y
        f["bed_true"] = bed_true
        f["bed_background"] = bed_bg
        f["bed_kriged"] = bed_kriged
        f["bed_var"] = resid_var
        f["obs_mask"] = obs_mask.astype(np.uint8)
        f["bed_ens"] = bed_ens
        f.attrs["track_stride_km"] = float(stride_km)
        f.attrs["track_half_width_m"] = float(track_half_width_m)
        f.attrs["grounded_observations_only"] = bool(
            grounded_observations_only
        )
        f.attrs["background_length_m"] = float(background_length_m)
        f.attrs["background_truncation_sigma"] = float(
            background_truncation_sigma
        )

    print("[bed-kriging] Done.")
    print(output_file)


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------
def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--icesee-path", default="./")
    parser.add_argument("--data-path", default="_modelrun_datasets")
    parser.add_argument("--Ne", type=int, default=60)
    parser.add_argument("--stride-km", type=float, default=10.0)
    parser.add_argument(
        "--track-half-width-km",
        type=float,
        default=2.5,
        help="Half-width of each constant-x radar corridor",
    )
    parser.add_argument(
        "--include-floating-observations",
        action="store_true",
        help="Also condition the kriging on floating-shelf bed observations",
    )
    parser.add_argument("--snap-idx", type=int, default=0)
    parser.add_argument("--sigma-noise", type=float, default=10.0)
    parser.add_argument("--sill-bed", type=float, default=5000.0)
    parser.add_argument("--range-bed", type=float, default=40000.0)
    parser.add_argument("--nugget-bed", type=float, default=50.0)
    parser.add_argument("--ensemble-inflation", type=float, default=1.5)
    parser.add_argument(
        "--background-length-km",
        type=float,
        default=40.0,
        help="Gaussian length scale for the 2-D bed background",
    )
    parser.add_argument(
        "--background-truncation-sigma",
        type=float,
        default=3.0,
        help="Truncate Gaussian smoothing beyond this many length scales",
    )
    parser.add_argument("--enforce-obs", action="store_true")
    parser.add_argument(
        "--use-floating-mask",
        action="store_true",
        help=(
            "Hold deep-bed points at the smooth background prior; never "
            "replace them with the hidden true bed"
        ),
    )
    parser.add_argument("--seed-base", type=int, default=1234)
    parser.add_argument(
        "--output-file",
        default=None,
        help="Output HDF5 path; run_da_issm.py expects bed_kriging_results.h5 in this directory",
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    generate_bed_kriging(
        icesee_path=args.icesee_path,
        data_path=args.data_path,
        Ne=args.Ne,
        stride_km=args.stride_km,
        track_half_width_m=1000.0 * args.track_half_width_km,
        grounded_observations_only=not args.include_floating_observations,
        snap_idx=args.snap_idx,
        sigma_noise=args.sigma_noise,
        sill_bed=args.sill_bed,
        range_bed=args.range_bed,
        nugget_bed=args.nugget_bed,
        ensemble_inflation=args.ensemble_inflation,
        background_length_m=1000.0 * args.background_length_km,
        background_truncation_sigma=args.background_truncation_sigma,
        enforce_obs=args.enforce_obs,
        use_floating_mask=args.use_floating_mask,
        seed_base=args.seed_base,
        output_file=args.output_file,
    )


# sample run command:
# python3 generate_bed_kringing.py \
#   --Ne 80 \
#   --stride-km 5 \
#   --sigma-noise 10 \
#   --sill-bed 12000 \
#   --range-bed 60000 \
#   --nugget-bed 50 \
#   --ensemble-inflation 2.5 \
#   --use-floating-mask
