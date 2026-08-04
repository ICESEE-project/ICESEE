#!/usr/bin/env python3
"""Audit the initial prior before spending time on a full DA experiment."""

from __future__ import annotations

import argparse
import os
import tempfile
from pathlib import Path

import h5py
os.environ["MPLCONFIGDIR"] = str(Path(tempfile.gettempdir()) / "icesee-mpl")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np
from scipy.interpolate import LinearNDInterpolator


def read_vector(path: Path, key: str) -> np.ndarray:
    with h5py.File(path, "r") as handle:
        return np.asarray(handle[key][...], dtype=float).reshape(-1)


def summarize_error(label: str, delta: np.ndarray, mask: np.ndarray) -> dict[str, float]:
    values = np.abs(delta[mask])
    metrics = {
        "rmse": float(np.sqrt(np.mean(delta[mask] ** 2))),
        "median": float(np.median(values)),
        "p95": float(np.percentile(values, 95)),
        "maximum": float(np.max(values)),
        "exact": float(np.mean(values == 0.0)),
        "within_5m": float(np.mean(values <= 5.0)),
    }
    print(
        f"{label}: RMSE={metrics['rmse']:.2f} m, "
        f"median |error|={metrics['median']:.2f} m, "
        f"p95={metrics['p95']:.2f} m, max={metrics['maximum']:.2f} m, "
        f"exact truth={100 * metrics['exact']:.1f}%, "
        f"within 5 m={100 * metrics['within_5m']:.1f}%"
    )
    return metrics


def centerline_grounding_x(
    x: np.ndarray,
    y: np.ndarray,
    phi: np.ndarray,
    y_center: float,
) -> float:
    """Return the grounded-to-floating zero crossing nearest the domain center."""
    x_line = np.linspace(np.min(x), np.max(x), 5000)
    y_line = np.full_like(x_line, y_center)
    values = np.asarray(
        LinearNDInterpolator(np.column_stack((x, y)), phi)(
            np.column_stack((x_line, y_line))
        ),
        dtype=float,
    )
    valid = np.isfinite(values)
    crossing = np.flatnonzero(
        valid[:-1]
        & valid[1:]
        & np.isfinite(values[:-1])
        & np.isfinite(values[1:])
        & (values[:-1] >= 0.0)
        & (values[1:] < 0.0)
    )
    if not crossing.size:
        return float("nan")
    candidates = []
    for index in crossing:
        left = values[index]
        right = values[index + 1]
        fraction = left / max(left - right, np.finfo(float).eps)
        candidates.append(x_line[index] + fraction * (x_line[index + 1] - x_line[index]))
    domain_center = 0.5 * (np.min(x) + np.max(x))
    return float(min(candidates, key=lambda value: abs(value - domain_center)))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "data_dir",
        nargs="?",
        default="_modelrun_datasets_rebuttal_ic_stretched_seed_bed_v6",
        type=Path,
    )
    parser.add_argument(
        "--expect-seaward-gl",
        action="store_true",
        help=(
            "Require the prior centerline grounding line to start seaward "
            "(larger x) of truth instead of preserving the true topology"
        ),
    )
    args = parser.parse_args()
    data_dir = args.data_dir

    true_file = data_dir / "ensemble_true_state_0.h5"
    prior_file = data_dir / "ensemble_nurged_state_0.h5"
    mesh_file = data_dir / "mesh_idxy_0.h5"
    for path in (true_file, prior_file, mesh_file):
        if not path.is_file():
            parser.error(f"missing preflight output: {path}")

    fields = {}
    for name, key in (
        ("thickness", "Thickness_1"),
        ("bed", "bed_1"),
        ("friction", "coefficient_1"),
    ):
        fields[("true", name)] = read_vector(true_file, key)
        fields[("prior", name)] = read_vector(prior_file, key)

    x = read_vector(mesh_file, "fric_x")
    y = read_vector(mesh_file, "fric_y")
    h_true = fields[("true", "thickness")]
    h_prior = fields[("prior", "thickness")]
    b_true = fields[("true", "bed")]
    b_prior = fields[("prior", "bed")]
    density_ratio = 917.0 / 1023.0
    true_ice = h_true > 1.0
    true_phi = h_true + b_true / density_ratio
    prior_phi = h_prior + b_prior / density_ratio
    true_grounded = true_ice & (true_phi >= 0.0)
    true_floating = true_ice & ~true_grounded
    prior_grounded = true_ice & (prior_phi >= 0.0)
    prior_floating = true_ice & ~prior_grounded

    bed_delta = b_prior - b_true
    thickness_delta = h_prior - h_true
    grounded_bed = summarize_error("grounded bed", bed_delta, true_grounded)
    floating_bed = summarize_error("floating bed", bed_delta, true_floating)
    summarize_error("grounded thickness", thickness_delta, true_grounded)
    summarize_error("floating thickness", thickness_delta, true_floating)

    floating_to_grounded = true_floating & prior_grounded
    grounded_to_floating = true_grounded & prior_floating
    topology_mismatch = true_ice & (true_grounded != prior_grounded)
    print(
        "flotation topology: "
        f"truth grounded={np.count_nonzero(true_grounded)}, "
        f"prior grounded={np.count_nonzero(prior_grounded)}, "
        f"floating->grounded={np.count_nonzero(floating_to_grounded)}, "
        f"grounded->floating={np.count_nonzero(grounded_to_floating)}, "
        f"mismatch={100 * np.mean(topology_mismatch[true_ice]):.2f}%"
    )
    y_center = 0.5 * (float(np.min(y)) + float(np.max(y)))
    true_gl_x = centerline_grounding_x(x, y, true_phi, y_center)
    prior_gl_x = centerline_grounding_x(x, y, prior_phi, y_center)
    gl_offset_km = (prior_gl_x - true_gl_x) / 1000.0
    print(
        "centerline grounding line: "
        f"truth x={true_gl_x / 1000.0:.2f} km, "
        f"prior x={prior_gl_x / 1000.0:.2f} km, "
        f"prior-truth={gl_offset_km:+.2f} km"
    )

    fig, axes = plt.subplots(2, 3, figsize=(16, 6.2), constrained_layout=True)
    panels = (
        (fields[("true", "bed")], "True initial bed (m)", "viridis", None),
        (
            bed_delta,
            "Prior - truth bed (m)",
            "seismic",
            "symmetric",
        ),
        (
            thickness_delta,
            "Prior - truth thickness (m)",
            "seismic",
            "symmetric",
        ),
        (
            fields[("true", "friction")],
            "True initial friction coefficient",
            "viridis",
            None,
        ),
        (
            fields[("prior", "friction")] - fields[("true", "friction")],
            "Prior - truth friction coefficient",
            "seismic",
            "symmetric",
        ),
        (
            prior_phi,
            "Prior flotation function (m)",
            "seismic",
            "symmetric",
        ),
    )
    for ax, (values, title, cmap, scaling) in zip(axes.flat, panels):
        kwargs = {}
        if scaling == "symmetric":
            limit = np.nanpercentile(np.abs(values), 99)
            kwargs = {"vmin": -limit, "vmax": limit}
        nx, ny = 650, 120
        grid_x = np.linspace(np.min(x), np.max(x), nx)
        grid_y = np.linspace(np.min(y), np.max(y), ny)
        mesh_x, mesh_y = np.meshgrid(grid_x, grid_y)
        sample_points = np.column_stack((mesh_x.ravel(), mesh_y.ravel()))
        points = np.column_stack((x, y))
        grid_values = LinearNDInterpolator(points, values)(sample_points).reshape(ny, nx)
        grid_true_phi = LinearNDInterpolator(points, true_phi)(sample_points).reshape(ny, nx)
        grid_prior_phi = LinearNDInterpolator(points, prior_phi)(sample_points).reshape(ny, nx)
        image = ax.pcolormesh(
            mesh_x, mesh_y, grid_values, shading="auto", cmap=cmap, **kwargs
        )
        ax.contour(mesh_x, mesh_y, grid_true_phi, levels=[0], colors="k", linewidths=1.2)
        ax.contour(
            mesh_x,
            mesh_y,
            grid_prior_phi,
            levels=[0],
            colors="cyan",
            linewidths=1.2,
            linestyles="--",
        )
        ax.set_title(title)
        ax.set_xlabel("x (km)")
        ax.set_ylabel("y (km)")
        km_formatter = FuncFormatter(lambda value, _position: f"{value / 1000:g}")
        ax.xaxis.set_major_formatter(km_formatter)
        ax.yaxis.set_major_formatter(km_formatter)
        fig.colorbar(image, ax=ax)

    output = data_dir / "initial_prior_audit.png"
    fig.savefig(output, dpi=180)
    print(f"wrote {output}")

    failures = []
    if floating_bed["exact"] > 0.05:
        failures.append("more than 5% of the floating bed is copied exactly from truth")
    if floating_bed["rmse"] > 100.0:
        failures.append("floating-bed RMSE exceeds the 100 m design ceiling")
    floating_bed_maximum = 200.0 if args.expect_seaward_gl else 125.0
    if floating_bed["maximum"] > floating_bed_maximum:
        failures.append(
            f"floating-bed maximum error exceeds {floating_bed_maximum:g} m"
        )
    if args.expect_seaward_gl:
        if not np.isfinite(gl_offset_km) or gl_offset_km <= 0.0:
            failures.append("the prior centerline GL is not seaward of truth")
        elif gl_offset_km > 40.0:
            failures.append(
                "the prior centerline GL is more than 40 km seaward of truth"
            )
        if np.count_nonzero(floating_to_grounded) <= np.count_nonzero(grounded_to_floating):
            failures.append(
                "the prior does not have a net seaward grounding-topology bias"
            )
    elif np.any(floating_to_grounded):
        failures.append("the prior grounds vertices that are floating in truth")
    if grounded_bed["rmse"] < 20.0:
        failures.append("grounded-bed prior is too close to truth to challenge the filter")

    if failures:
        for failure in failures:
            print(f"FAIL: {failure}")
        return 2
    if args.expect_seaward_gl:
        print(
            "PASS: independent bed error is moderate and the initial "
            "centerline grounding line is intentionally seaward of truth"
        )
    else:
        print(
            "PASS: independent tapered floating-bed error is moderate and "
            "the true floating domain remains floating"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
