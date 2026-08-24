#!/usr/bin/env python3
"""
bed_roughness.py

Generate 2-D midpoint-displacement (diamond-square) roughness on a regular grid
(0..Lx, 0..Ly) at dx resolution, then interpolate it onto ISSM coordinates stored in:
  /fric_x, /fric_y  (1D arrays)

Writes an HDF5 file containing:
  /bed_roughness  (1D, float32) aligned with the ordering of the ISSM points.

Usage:
  python bed_roughness.py \
    --mesh_h5 _modelrun_datasets/mesh_idxy_0.h5 \
    --out_h5 bed_roughness_flat.h5 \
    --sigma0_m 500 --H 0.7 --h 0.7 --nrec 10 --dx_m 100 --seed 1
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import numpy as np

try:
    import h5py  # type: ignore
except Exception:
    h5py = None


@dataclass
class Params:
    Lx_m: int = 640_000
    Ly_m: int = 80_000
    dx_m: int = 100
    nrec: int = 10
    sigma0_m: float = 500.0
    h: float = 0.7
    H: float = 0.7
    seed: int = 1234


def _diamond_square_rect(z: np.ndarray, step: int, sigma: float, rng: np.random.Generator) -> None:
    """
    One rectangular diamond-square refinement step (in-place).

    Critical detail: the square step MUST average only finite neighbors; otherwise NaNs
    propagate and the grid never fills (which is what you observed).
    """
    ny, nx = z.shape
    half = step // 2

    # ---------------- Diamond step: set centers ----------------
    for y0 in range(0, ny - 1, step):
        y1 = y0 + step
        yc = y0 + half
        for x0 in range(0, nx - 1, step):
            x1 = x0 + step
            xc = x0 + half

            c00 = z[y0, x0]
            c10 = z[y0, x1]
            c01 = z[y1, x0]
            c11 = z[y1, x1]

            # corners should already be finite at this level; but guard anyway
            corners = np.array([c00, c10, c01, c11], dtype=float)
            good = np.isfinite(corners)
            if not np.any(good):
                continue

            mean_corners = corners[good].mean()
            z[yc, xc] = mean_corners + rng.normal(0.0, sigma)

    # ---------------- Square step: set edge midpoints ----------------
    for y in range(0, ny, half):
        x_start = half if (y // half) % 2 == 0 else 0
        for x in range(x_start, nx, step):
            if not np.isnan(z[y, x]):
                continue

            vals = []
            if y - half >= 0 and np.isfinite(z[y - half, x]):
                vals.append(z[y - half, x])
            if y + half < ny and np.isfinite(z[y + half, x]):
                vals.append(z[y + half, x])
            if x - half >= 0 and np.isfinite(z[y, x - half]):
                vals.append(z[y, x - half])
            if x + half < nx and np.isfinite(z[y, x + half]):
                vals.append(z[y, x + half])

            if len(vals) == 0:
                continue

            z[y, x] = (np.mean(vals) + rng.normal(0.0, sigma))


def generate_midpoint_displacement_2d(p: Params) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Generate roughness br(y,x) over [0,Lx]x[0,Ly] at dx_m resolution.

    sigma_k = sigma0 * H / (2**(h*k)), k = 0..nrec-1
    """
    rng = np.random.default_rng(p.seed)

    nx_t = p.Lx_m // p.dx_m + 1
    ny_t = p.Ly_m // p.dx_m + 1

    step0 = 2 ** p.nrec
    nx_work = (int(np.ceil((nx_t - 1) / step0)) * step0) + 1
    ny_work = (int(np.ceil((ny_t - 1) / step0)) * step0) + 1

    z = np.full((ny_work, nx_work), np.nan, dtype=float)

    # Initialize corners
    z[0, 0] = rng.normal(0.0, p.sigma0_m * p.H)
    z[0, -1] = rng.normal(0.0, p.sigma0_m * p.H)
    z[-1, 0] = rng.normal(0.0, p.sigma0_m * p.H)
    z[-1, -1] = rng.normal(0.0, p.sigma0_m * p.H)

    # Refine
    step = step0
    for k in range(p.nrec):
        sigma_k = (p.sigma0_m * p.H) / (2.0 ** (p.h * k))
        _diamond_square_rect(z, step=step, sigma=sigma_k, rng=rng)
        step //= 2

    # Crop
    br = z[:ny_t, :nx_t].copy()

    # Coordinates
    x = np.arange(nx_t, dtype=float) * p.dx_m
    y = np.arange(ny_t, dtype=float) * p.dx_m

    # If any NaNs remain, set to 0 (should be very small now)
    n_nan = int(np.isnan(br).sum())
    n_inf = int(np.isinf(br).sum())
    if n_nan or n_inf:
        print(f"[WARN] br contains NaN/Inf after crop: NaN={n_nan}, Inf={n_inf} -> replacing with 0")
    br = np.nan_to_num(br, nan=0.0, posinf=0.0, neginf=0.0)

    # Zero-mean
    br -= br.mean()

    return br, x, y


def bilinear_interp_regular_grid(
    br: np.ndarray, x: np.ndarray, y: np.ndarray, xq: np.ndarray, yq: np.ndarray
) -> np.ndarray:
    """Vectorized bilinear interpolation of br(y,x) to query points (xq,yq)."""
    dx = float(x[1] - x[0])
    dy = float(y[1] - y[0])

    xq_c = np.clip(xq, x[0], x[-1])
    yq_c = np.clip(yq, y[0], y[-1])

    ix = (xq_c - x[0]) / dx
    iy = (yq_c - y[0]) / dy

    i0 = np.floor(ix).astype(int)
    j0 = np.floor(iy).astype(int)

    i0 = np.clip(i0, 0, br.shape[1] - 2)
    j0 = np.clip(j0, 0, br.shape[0] - 2)

    tx = ix - i0
    ty = iy - j0

    f00 = br[j0, i0]
    f10 = br[j0, i0 + 1]
    f01 = br[j0 + 1, i0]
    f11 = br[j0 + 1, i0 + 1]

    out = (1 - tx) * (1 - ty) * f00 + tx * (1 - ty) * f10 + (1 - tx) * ty * f01 + tx * ty * f11
    return out.astype(np.float32)


def load_fric_xy(mesh_h5: Path) -> tuple[np.ndarray, np.ndarray]:
    if h5py is None:
        raise RuntimeError("h5py is not installed. Install with: pip install h5py")

    with h5py.File(mesh_h5, "r") as f:
        if "/fric_x" not in f or "/fric_y" not in f:
            raise KeyError(f"Expected datasets '/fric_x' and '/fric_y' in {mesh_h5}")
        x_param = np.asarray(f["/fric_x"][:], dtype=float).ravel()
        y_param = np.asarray(f["/fric_y"][:], dtype=float).ravel()

    if x_param.size != y_param.size:
        raise ValueError(f"fric_x and fric_y size mismatch: {x_param.size} vs {y_param.size}")
    return x_param, y_param


def save_bed_roughness_h5(out_h5: Path, bed_roughness: np.ndarray, p: Params, src_mesh_h5: Path) -> None:
    if h5py is None:
        raise RuntimeError("h5py is not installed. Install with: pip install h5py")

    out_h5.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(out_h5, "w") as f:
        f.create_dataset("bed_roughness", data=bed_roughness.astype(np.float32),
                         compression="gzip", compression_opts=4)

        f.attrs["source_mesh_h5"] = str(src_mesh_h5)
        f.attrs["npoints"] = int(bed_roughness.size)
        f.attrs["dx_m"] = int(p.dx_m)
        f.attrs["Lx_m"] = int(p.Lx_m)
        f.attrs["Ly_m"] = int(p.Ly_m)
        f.attrs["nrec"] = int(p.nrec)
        f.attrs["sigma0_m"] = float(p.sigma0_m)
        f.attrs["h"] = float(p.h)
        f.attrs["H"] = float(p.H)
        f.attrs["seed"] = int(p.seed)


def main() -> None:
    ap = argparse.ArgumentParser()

    ap.add_argument("--Lx_km", type=float, default=640.0)
    ap.add_argument("--Ly_km", type=float, default=80.0)
    ap.add_argument("--dx_m", type=int, default=100)
    ap.add_argument("--nrec", type=int, default=10)
    ap.add_argument("--sigma0_m", type=float, default=500.0)
    ap.add_argument("--h", type=float, default=0.7)
    ap.add_argument("--H", type=float, default=0.7)
    ap.add_argument("--seed", type=int, default=1234)

    ap.add_argument("--mesh_h5", type=str, required=True,
                    help="Path to mesh_idxy_0.h5 containing /fric_x and /fric_y")
    ap.add_argument("--out_h5", type=str, default="bed_roughness_flat.h5")

    args = ap.parse_args()

    p = Params(
        Lx_m=int(round(args.Lx_km * 1000)),
        Ly_m=int(round(args.Ly_km * 1000)),
        dx_m=int(args.dx_m),
        nrec=int(args.nrec),
        sigma0_m=float(args.sigma0_m),
        h=float(args.h),
        H=float(args.H),
        seed=int(args.seed),
    )

    mesh_h5 = Path(args.mesh_h5)
    out_h5 = Path(args.out_h5)

    br, x, y = generate_midpoint_displacement_2d(p)

    x_param, y_param = load_fric_xy(mesh_h5)

    print("fric_x: nan=", int(np.isnan(x_param).sum()), "inf=", int(np.isinf(x_param).sum()),
          "min=", float(np.nanmin(x_param)), "max=", float(np.nanmax(x_param)))
    print("fric_y: nan=", int(np.isnan(y_param).sum()), "inf=", int(np.isinf(y_param).sum()),
          "min=", float(np.nanmin(y_param)), "max=", float(np.nanmax(y_param)))

    good = np.isfinite(x_param) & np.isfinite(y_param)
    if not np.all(good):
        print(f"[WARN] Dropping {good.size - good.sum()} non-finite (fric_x, fric_y) points")
    x_param = x_param[good]
    y_param = y_param[good]

    print(f"Grid x range: [{x[0]}, {x[-1]}] m ; Grid y range: [{y[0]}, {y[-1]}] m")

    bed_roughness_1d = bilinear_interp_regular_grid(br, x, y, x_param, y_param)
    bed_roughness_1d = np.nan_to_num(bed_roughness_1d, nan=0.0, posinf=0.0, neginf=0.0)
    bed_roughness_1d -= float(bed_roughness_1d.mean())

    print("Roughness after interpolation:")
    print(" std =", float(bed_roughness_1d.std(ddof=1)), "m")
    print(" min =", float(bed_roughness_1d.min()), "m")
    print(" max =", float(bed_roughness_1d.max()), "m")

    save_bed_roughness_h5(out_h5, bed_roughness_1d, p, src_mesh_h5=mesh_h5)

    print(f"Wrote bed roughness aligned to ISSM points: {out_h5}")
    print(f"npoints = {bed_roughness_1d.size}")


if __name__ == "__main__":
    main()