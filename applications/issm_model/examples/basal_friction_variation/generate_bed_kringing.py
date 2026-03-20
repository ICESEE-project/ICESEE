#!/usr/bin/env python3
# -*- @author: brian kyanjo -*-
# -*- @date: December 2025 -*-
"""
Generate a kriging-based conditional random bed field ensemble.

- Reads:
    {icesee_path}/{data_path}/mesh_idxy_0.h5        -> /fric_x, /fric_y
    {icesee_path}/{data_path}/true_nurged_states.h5 -> /true_state

- Assumes:
    true_state has shape (6 * ndim, nt)
    bed_true_all = true_state[4*ndim : 5*ndim, :] (bedrock)

- Produces:
    {icesee_path}/{data_path}/bed_kriging_results.h5 with datasets:
        x, y, bed_true, obs_mask, z_obs,
        bed_kriged, bed_var, bed_ens
"""

import argparse
import os
import sys
from typing import Optional

import h5py
import numpy as np
import gstools as gs
from gstools import krige

# Try to use tqdm if available, otherwise fall back to simple logging
try:
    from tqdm import trange
    HAVE_TQDM = True
except ImportError:
    HAVE_TQDM = False


def generate_bed_kriging(
    icesee_path: str = "./",
    data_path: str = "_modelrun_datasets",
    Ne: int = 60,
    stride_km: float = 10.0,
    snap_idx: int = 0,
    sigma_noise: float = 15.0,
    sill_bed: float = 5000.0,
    range_bed: float = 30000.0,
    nugget_bed: float = 100.0,
    seed_base: int = 1234,
    output_file: Optional[str] = None,
) -> None:
    # ------------------------------------------------------------------
    # Paths
    # ------------------------------------------------------------------
    base_dir = os.path.join(icesee_path, data_path)
    mesh_file = os.path.join(base_dir, "mesh_idxy_0.h5")
    true_file = os.path.join(base_dir, "true_nurged_states.h5")

    if output_file is None:
        output_file = os.path.join(base_dir, "bed_kriging_results.h5")

    print("[bed-kriging] ================================================")
    print("[bed-kriging] Starting kriging-based bed ensemble generation")
    print("[bed-kriging] -----------------------------------------------")
    print(f"[bed-kriging] icesee_path   : {icesee_path}")
    print(f"[bed-kriging] data_path     : {data_path}")
    print(f"[bed-kriging] mesh_file     : {mesh_file}")
    print(f"[bed-kriging] true_file     : {true_file}")
    print(f"[bed-kriging] output_file   : {output_file}")
    print(f"[bed-kriging] Ne            : {Ne}")
    print(f"[bed-kriging] stride_km     : {stride_km}")
    print(f"[bed-kriging] snap_idx      : {snap_idx}")
    print("[bed-kriging] ================================================")
    sys.stdout.flush()

    # ------------------------------------------------------------------
    # Read mesh coordinates (node positions, in meters)
    # ------------------------------------------------------------------
    print("[bed-kriging] [1/6] Reading mesh coordinates ...", end="", flush=True)
    with h5py.File(mesh_file, "r") as f:
        x_param = f["/fric_x"][:]  # shape (fdim,)
        y_param = f["/fric_y"][:]

    x = np.asarray(x_param).ravel()  # meters
    y = np.asarray(y_param).ravel()  # meters
    npts = x.size
    print(" done.")
    print(f"[bed-kriging]        Number of mesh points: {npts}")
    sys.stdout.flush()

    # ------------------------------------------------------------------
    # Read true state and extract bed at one time snapshot
    # ------------------------------------------------------------------
    print("[bed-kriging] [2/6] Reading true state and extracting bed ...", end="", flush=True)
    with h5py.File(true_file, "r") as f:
        w = f["true_state"][:]  # shape (6*ndim, nt)

    ndim = w.shape[0] // 6
    nt = w.shape[1]

    if ndim != npts:
        raise RuntimeError(
            f"[bed-kriging] ndim from true_state ({ndim}) "
            f"!= number of mesh points ({npts})"
        )

    if snap_idx < 0 or snap_idx >= nt:
        raise ValueError(
            f"[bed-kriging] snap_idx={snap_idx} out of range [0, {nt-1}]"
        )

    bed_true_all = w[4 * ndim : 5 * ndim, :]  # (ndim, nt)
    bed_true = bed_true_all[:, snap_idx]      # (ndim,)
    print(" done.")
    print(f"[bed-kriging]        nt snapshots available: {nt}")
    sys.stdout.flush()

    # ------------------------------------------------------------------
    # Radar-track mask: tracks perpendicular to flow along x
    # ------------------------------------------------------------------
    print("[bed-kriging] [3/6] Building radar-track observation mask ...", end="", flush=True)
    # Convert to km for spacing logic only; kriging still uses meters
    x_km = x / 1000.0

    stride = float(stride_km)
    x_min, x_max = x_km.min(), x_km.max()
    track_xs = np.arange(x_min, x_max + 1e-6, stride)

    if track_xs.size > 1:
        dx_nom = (x_max - x_min) / max(track_xs.size - 1, 1)
    else:
        dx_nom = stride
    band = 0.5 * dx_nom  # km

    obs_mask = np.zeros(npts, dtype=bool)
    for x_line in track_xs:
        obs_mask |= np.abs(x_km - x_line) <= band

    n_obs = int(obs_mask.sum())
    print(" done.")
    print(f"[bed-kriging]        Radar-bed obs points: {n_obs} / {npts}")
    if n_obs == 0:
        raise RuntimeError("[bed-kriging] No observation points selected; check stride_km / domain.")
    sys.stdout.flush()

    # ------------------------------------------------------------------
    # Synthetic radar observations = true bed + noise
    # ------------------------------------------------------------------
    print("[bed-kriging] [4/6] Generating synthetic radar observations ...", end="", flush=True)
    rng = np.random.default_rng(seed_base)

    x_obs = x[obs_mask]
    y_obs = y[obs_mask]
    z_obs = bed_true[obs_mask] + sigma_noise * rng.standard_normal(n_obs)
    print(" done.")
    sys.stdout.flush()

    # ------------------------------------------------------------------
    # Kriging model (exponential covariance in *meters*)
    # ------------------------------------------------------------------
    print("[bed-kriging] [5/6] Fitting kriging model and computing mean field ...", end="", flush=True)
    model_bed = gs.Exponential(
        dim=2,
        var=sill_bed,
        len_scale=range_bed,  # in meters, same units as x,y
        nugget=nugget_bed,
    )

    # Ordinary kriging conditioned on radar tracks
    ok_bed = krige.Ordinary(
        model_bed,
        cond_pos=[x_obs, y_obs],
        cond_val=z_obs,
    )

    # Kriged mean + variance at all nodes
    bed_kriged, bed_var = ok_bed([x, y], return_var=True)
    bed_kriged = np.asarray(bed_kriged)
    bed_var = np.asarray(bed_var)
    print(" done.")
    sys.stdout.flush()

    # ------------------------------------------------------------------
    # Conditional random field ensemble
    # ------------------------------------------------------------------
    print("[bed-kriging] [6/6] Drawing conditional ensemble realizations ...")
    cond_srf_bed = gs.CondSRF(ok_bed)

    bed_ens = np.zeros((Ne, npts))

    if HAVE_TQDM:
        # Fancy progress bar if tqdm is installed
        for e in trange(Ne, desc="[bed-kriging]   ensemble", unit="ens"):
            seed = seed_base + e
            field = cond_srf_bed([x, y], seed=seed)
            bed_ens[e, :] = np.asarray(field)
    else:
        # Simple textual progress if tqdm isn't available
        for e in range(Ne):
            seed = seed_base + e
            field = cond_srf_bed([x, y], seed=seed)
            bed_ens[e, :] = np.asarray(field)

            # Print occasional progress updates
            if Ne <= 10 or e == Ne - 1 or (e + 1) % max(1, Ne // 10) == 0:
                pct = 100.0 * (e + 1) / Ne
                print(f"[bed-kriging]   ensemble {e+1}/{Ne} ({pct:5.1f}%)", flush=True)

    # ------------------------------------------------------------------
    # Write results
    # ------------------------------------------------------------------
    print("[bed-kriging] Writing results to HDF5 ...", end="", flush=True)
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with h5py.File(output_file, "w") as f:
        f.create_dataset("x", data=x)
        f.create_dataset("y", data=y)
        f.create_dataset("bed_true", data=bed_true)
        f.create_dataset("obs_mask", data=obs_mask)
        f.create_dataset("z_obs", data=z_obs)
        f.create_dataset("bed_kriged", data=bed_kriged)
        f.create_dataset("bed_var", data=bed_var)
        f.create_dataset("bed_ens", data=bed_ens)
    print(" done.")
    print("[bed-kriging] ================================================")
    print("[bed-kriging] Finished kriging-based bed ensemble generation")
    print("[bed-kriging] ================================================")
    sys.stdout.flush()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate kriging-based conditional random bed ensemble "
                    "from radar-like tracks."
    )
    parser.add_argument(
        "--icesee-path",
        type=str,
        default="./",
        help="Root ICESEE application path "
             "(e.g. /path/to/ICESEE/applications/issm_model/examples/ISMIP_Choi)",
    )
    parser.add_argument(
        "--data-path",
        type=str,
        default="_modelrun_datasets",
        help="Relative data directory under icesee_path (default: _modelrun_datasets)",
    )
    parser.add_argument(
        "--Ne",
        type=int,
        default=60,
        help="Number of ensemble realizations to generate (default: 60)",
    )
    parser.add_argument(
        "--stride-km",
        type=float,
        default=7.0,
        help="Radar track spacing in km (default: 7.0)",
    )
    parser.add_argument(
        "--snap-idx",
        type=int,
        default=0,
        help="Time snapshot index into true_state (0-based, default: 0)",
    )
    parser.add_argument(
        "--sigma-noise",
        type=float,
        default=10.0,
        help="Standard deviation of radar measurement noise in meters (default: 10.0)",
    )
    parser.add_argument(
        "--sill-bed",
        type=float,
        default=5000.0,
        help="Covariance sill (variance) for bed in m^2 (default: 5000.0)",
    )
    parser.add_argument(
        "--range-bed",
        type=float,
        default=30000.0,
        help="Covariance range (correlation length) for bed in meters (default: 30000.0)",
    )
    parser.add_argument(
        "--nugget-bed",
        type=float,
        default=100.0,
        help="Nugget (measurement error / microscale variance) in m^2 (default: 100.0)",
    )
    parser.add_argument(
        "--seed-base",
        type=int,
        default=1234,
        help="Base random seed for noise and ensemble (default: 1234)",
    )
    parser.add_argument(
        "--output-file",
        type=str,
        default=None,
        help=("Output HDF5 file path. "
              "Default: {icesee_path}/{data_path}/bed_kriging_results.h5"),
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    generate_bed_kriging(
        icesee_path=args.icesee_path,
        data_path=args.data_path,
        Ne=args.Ne,
        stride_km=args.stride_km,
        snap_idx=args.snap_idx,
        sigma_noise=args.sigma_noise,
        sill_bed=args.sill_bed,
        range_bed=args.range_bed,
        nugget_bed=args.nugget_bed,
        seed_base=args.seed_base,
        output_file=args.output_file,
    )