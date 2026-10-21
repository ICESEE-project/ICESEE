"""Icepack example views
(applications/icepack_model/examples/synthetic_ice_stream/read_results.ipynb).

State layout (as in the notebook): four blocks of ``hdim = nd // 4`` rows --
thickness ``h``, velocity components ``u`` and ``v``, surface mass balance
``smb`` -- each in the scalar function space's degree-of-freedom ordering.

The notebook draws maps with ``firedrake.tripcolor`` on the Firedrake mesh.
Firedrake is not available where the launcher renders, so maps here use
matplotlib's ``tripcolor`` (Gouraud shading) on the exported degree-of-
freedom coordinates (``dof_coordinates.h5``, written in the Icepack
environment by :mod:`ICESEE.src.postprocessing.icepack_mesh`). The plotted
values, colour maps, colour limits, labels and units are the notebook's; the
triangulation between nodes is matplotlib's.
"""
from __future__ import annotations

import string

import numpy as np
from matplotlib import gridspec
from matplotlib.figure import Figure
from matplotlib.ticker import MaxNLocator, StrMethodFormatter

from .. import diagnostics as dx
from ..outputs import ResultsError, load_dof_coordinates
from . import Option, ViewContext, ViewResult, ViewSpec, register

#: icepack.constants.ice_density / water_density (917 / 1024)
RHO_I_OVER_RHO_W = 917.0 / 1024.0

#: the notebook's figure years for the difference maps and error tables
DESIRED_YEARS = [5, 15, 30, 60, 85]


def _hdim(ctx: ViewContext) -> int:
    return ctx.outputs.n_state // 4


def _grid(ctx: ViewContext) -> tuple[float, float, int, int]:
    meta = ctx.outputs.metadata
    if "Lxy" not in meta or "nxy" not in meta:
        raise ResultsError("missing", "true-wrong-icepack.h5 has no Lxy/nxy (domain size and mesh)")
    Lx, Ly = meta["Lxy"][0], meta["Lxy"][1]
    nx, ny = int(meta["nxy"][0]), int(meta["nxy"][1])
    return Lx, Ly, nx, ny


def _middle_index(ctx: ViewContext) -> int:
    Lx, Ly, nx, ny = _grid(ctx)
    return dx.grid_point_index(nx, ny, _hdim(ctx), "middle")


def _marker_xy(ctx: ViewContext, h_indx: int) -> tuple[float, float]:
    Lx, Ly, nx, ny = _grid(ctx)
    return dx.grid_point_xy(h_indx, nx, ny, _hdim(ctx), Lx, Ly)


def _max_index(ctx: ViewContext) -> int:
    return _hdim(ctx) - 1


def _max_step(ctx: ViewContext) -> int:
    o = ctx.outputs
    return min(o.t.size, o.ensemble_mean.shape[1]) - 1


def _n(ctx: ViewContext) -> int:
    o = ctx.outputs
    return min(o.t.size, o.true_state.shape[1], o.ensemble_mean.shape[1])


# ─────────────────────────────────────────────────────────────────────────────
# cells 5 & 6: truth vs. wrong at one grid point
# ─────────────────────────────────────────────────────────────────────────────
def truth_vs_wrong(ctx: ViewContext, opts: dict) -> ViewResult:
    o = ctx.outputs
    hdim = _hdim(ctx)
    h_indx = opts["h_indx"]
    n = _n(ctx)
    t = o.t[:n]
    fig = Figure(figsize=(10, 5))
    ax = fig.subplots()
    if opts["variable"] == "thickness":
        ax.plot(t, o.true_state[:hdim][h_indx, :n], "r", label="True")
        ax.plot(t, o.nurged_state[:hdim][h_indx, :n], "g", label="Wrong")
        obs = o.observations[h_indx, :]
        ax.plot(o.t[o.obs_index], obs, "kx", label="Observations")
        ax.set_xlabel("Time (years)")
        ax.set_ylabel("Ice thickness (m)")
        ax.legend()
        ax.grid()
        ax.set_title(f"middle profile at h_indx = {h_indx}")
    else:
        ax.plot(t, o.true_state[3 * hdim:][h_indx, :n], "r", label="True")
        ax.plot(t, o.nurged_state[3 * hdim:][h_indx, :n], "g", label="Wrong")
        ax.set_xlabel("Time (years)")
        ax.set_ylabel("SMB (m/yr)")
        ax.legend()
        ax.grid()
    return ViewResult(figure=fig, notes=list(o.warnings))


# ─────────────────────────────────────────────────────────────────────────────
# cell 14: thickness / velocity / SMB profiles with the ensemble
# ─────────────────────────────────────────────────────────────────────────────
def ensemble_profiles(ctx: ViewContext, opts: dict) -> ViewResult:
    o = ctx.outputs
    hdim = _hdim(ctx)
    h_indx = opts["h_indx"]
    filter_type = o.filter_type
    t = o.t
    n = _n(ctx)
    _t = t[:n]
    w = o.observations
    ind_m = o.obs_index

    fig = Figure(figsize=(14, 12))
    ax = fig.subplots(3, 1)

    # (a) thickness
    obs_h = w[h_indx, :]
    h_true = o.true_state[:hdim][h_indx, :n]
    h_nurged = o.nurged_state[:hdim][h_indx, :n]
    h_ens_mean = o.ensemble_mean[:hdim][h_indx, :n]
    members = o.ensemble_file is not None
    if members:
        h_ens_mem = o.ensemble_members(h_indx)[:n]
        for ens in range(h_ens_mem.shape[1]):
            ax[0].plot(_t, h_ens_mem[:, ens], lw=0.5, color="gray", label=f"{filter_type}" if ens == 0 else "")
    ax[0].plot(t[ind_m], obs_h, "kx", label="Observations", markersize=8)
    ax[0].plot(_t, h_true, "--r", label="True", linewidth=2)
    ax[0].plot(_t, h_nurged, "--g", label="Wrong", linewidth=2)
    ax[0].plot(_t, h_ens_mean, "--b", label=f"{filter_type} mean", linewidth=2)
    ax[0].legend(prop={"size": 12, "weight": "bold"})
    ax[0].set_xlim(opts["xlim"])
    ax[0].grid(True, linestyle="--", alpha=0.6)
    ax[0].set_ylabel("Thickness (m)", fontdict={"fontsize": 16, "fontweight": "bold"})

    # (b) velocity magnitude
    u_indx, v_indx = hdim + h_indx, 2 * hdim + h_indx
    uobs = dx.magnitude(w[u_indx, :], w[v_indx, :])
    utrue = dx.magnitude(o.true_state[u_indx, :n], o.true_state[v_indx, :n])
    unurge = dx.magnitude(o.nurged_state[u_indx, :n], o.nurged_state[v_indx, :n])
    uens_mean = dx.magnitude(o.ensemble_mean[u_indx, :n], o.ensemble_mean[v_indx, :n])
    if members:
        uens_mem = dx.magnitude(o.ensemble_members(u_indx)[:n], o.ensemble_members(v_indx)[:n])
        for ens in range(uens_mem.shape[1]):
            ax[1].plot(_t, uens_mem[:, ens], lw=0.5, color="gray", label=f"{filter_type}" if ens == 0 else "")
    ax[1].plot(_t, utrue, "--r", label="True", linewidth=2)
    ax[1].plot(_t, unurge, "--g", label="Wrong", linewidth=2)
    ax[1].plot(_t, uens_mean, "--b", label=f"{filter_type} mean", linewidth=2)
    ax[1].plot(t[ind_m], uobs, "kx", label="Observations", markersize=8)
    ax[1].legend(prop={"size": 12, "weight": "bold"})
    ax[1].set_ylabel("Velocity (m/year)", fontdict={"fontsize": 16, "fontweight": "bold"})
    ax[1].grid(True, linestyle="--", alpha=0.6)
    ax[1].set_xlim(opts["xlim"])

    # (c) surface mass balance
    smb_idx = 3 * hdim + h_indx
    if members:
        smb_ens = o.ensemble_members(smb_idx)[:n]
        for ens in range(smb_ens.shape[1]):
            ax[2].plot(_t, smb_ens[:, ens], lw=0.5, color="gray", label=f"{filter_type}" if ens == 0 else "")
    ax[2].plot(_t, o.ensemble_mean[smb_idx, :n], "-b", label=f"{filter_type} mean", linewidth=2)
    ax[2].plot(_t, o.true_state[smb_idx, :n], "--r", label="True", linewidth=2)
    ax[2].plot(_t, o.nurged_state[smb_idx, :n], "--g", label="Wrong", linewidth=2)
    ax[2].plot(t[ind_m], w[smb_idx, :], "kx", label="Observations", markersize=8)
    ax[2].set_xlabel("Time (years)", fontdict={"fontsize": 16, "fontweight": "bold"})
    ax[2].set_ylabel("SMB (m/year)", fontdict={"fontsize": 16, "fontweight": "bold"})
    ax[2].grid(True, linestyle="--", alpha=0.6)
    ax[2].set_xlim(opts["xlim"])
    ax[2].legend(prop={"size": 12, "weight": "bold"})

    for a in ax:
        a.tick_params(labelsize=14)
        for lab in (*a.get_xticklabels(), *a.get_yticklabels()):
            lab.set_fontweight("bold")
    for a, letter in zip(ax, "abc"):
        a.text(0.01, 0.9, f"({letter})", transform=a.transAxes, fontsize=16, fontweight="bold")
    notes = list(o.warnings)
    notes.append(f"Grid point h_indx = {h_indx} (state rows {h_indx}, {u_indx}, {v_indx}, {smb_idx}); "
                 f"{o.n_ens} ensemble members.")
    return ViewResult(figure=fig, notes=notes)


# ─────────────────────────────────────────────────────────────────────────────
# maps (cells 10, 16, 17) -- need the degree-of-freedom coordinates
# ─────────────────────────────────────────────────────────────────────────────
def _coords(ctx: ViewContext):
    xy = load_dof_coordinates(ctx.outputs.data_dir if ctx.data_dir is None else ctx.data_dir)
    hdim = _hdim(ctx)
    if xy.shape[0] != hdim:
        raise ResultsError("incomplete", f"dof_coordinates.h5 has {xy.shape[0]} nodes, the state {hdim} per field")
    from matplotlib.tri import Triangulation

    return xy, Triangulation(xy[:, 0], xy[:, 1])


def _fields(o, hdim: int, step: int):
    """get_firedrake_functions() of the notebook, as plain arrays."""
    h = [o.true_state[:hdim, step], o.nurged_state[:hdim, step], o.ensemble_mean[:hdim, step]]
    u = [dx.magnitude(s[hdim:2 * hdim, step], s[2 * hdim:3 * hdim, step])
         for s in (o.true_state, o.nurged_state, o.ensemble_mean)]
    smb = [o.true_state[3 * hdim:4 * hdim, step], o.nurged_state[3 * hdim:4 * hdim, step],
           o.ensemble_mean[3 * hdim:4 * hdim, step]]
    return h, u, smb


def snapshot_maps(ctx: ViewContext, opts: dict) -> ViewResult:
    """cell 10: thickness and velocity of the true, wrong and assimilated
    simulations at one time step, on common colour scales."""
    o = ctx.outputs
    hdim = _hdim(ctx)
    xy, tri = _coords(ctx)
    step = opts["step"]
    h, u, _ = _fields(o, hdim, step)
    h_min, h_max = h[0].min(), h[0].max()
    u_min, u_max = u[0].min(), u[0].max()
    x_h, y_h = _marker_xy(ctx, _middle_index(ctx))

    fig = Figure(figsize=(22.5, 5.2))
    axes = fig.subplots(2, 3)
    titles = ["True simulation ", "Wrong simulation", "Assimilated  simulation "]
    for col in range(3):
        ch = axes[0][col].tripcolor(tri, h[col], shading="gouraud", vmin=h_min, vmax=h_max)
        cu = axes[1][col].tripcolor(tri, u[col], shading="gouraud", vmin=u_min, vmax=u_max)
        axes[0][col].set_title(titles[col], fontsize=18, fontweight="bold")
    for row, colors, label in ((0, ch, "Thickness (m)"), (1, cu, "Velocity (m/year)")):
        cbar = fig.colorbar(colors, ax=axes[row][2], fraction=0.012, pad=0.04, label=label)
        cbar.set_label(label, fontsize=16, fontweight="bold")
        for lab in cbar.ax.yaxis.get_ticklabels():
            lab.set_fontweight("bold")
            lab.set_fontsize(18)
    for row in range(2):
        for col in range(3):
            ax = axes[row][col]
            ax.set_aspect("equal")
            ax.yaxis.set_major_formatter(StrMethodFormatter("{x:.0f}"))
            ax.tick_params(axis="both", which="major", labelsize=18)
            if col == 0:
                ax.set_ylabel("Y (m)", fontsize=18, fontweight="bold")
            else:
                ax.set_yticklabels([])
            if row == 1:
                ax.set_xlabel("X (m)", fontsize=18, fontweight="bold")
            else:
                ax.set_xticklabels([])
            ax.plot(x_h, y_h, marker="o", markersize=8, color="red", markeredgecolor="black", zorder=10)
            for lab in (*ax.get_xticklabels(), *ax.get_yticklabels()):
                lab.set_fontweight("bold")
    fig.tight_layout()
    notes = list(o.warnings) + [f"Time step {step} (t = {o.t[step]:g} years). Red dot: the profile point."]
    return ViewResult(figure=fig, notes=notes)


def difference_maps(ctx: ViewContext, opts: dict) -> ViewResult:
    """cell 16: true thickness / speed / SMB, then signed errors of the
    wrong run and of the assimilation at chosen years, on one symmetric
    colour scale per column."""
    o = ctx.outputs
    hdim = _hdim(ctx)
    xy, tri = _coords(ctx)
    t = o.t[:_n(ctx)]
    desired_years = list(opts["years"])
    steps = dx.nearest_steps(t, desired_years)
    no_of_rows = len(steps) + 2
    x_h, y_h = _marker_xy(ctx, _middle_index(ctx))

    h0, u0, s0 = _fields(o, hdim, 0)
    diffs = [[h0[1] - h0[0]], [u0[1] - u0[0]], [s0[1] - s0[0]]]
    per_step = []
    for step in steps:
        h, u, s = _fields(o, hdim, step)
        d = (h[2] - h[0], u[2] - u[0], s[2] - s[0])
        per_step.append(d)
        for col in range(3):
            diffs[col].append(d[col])
    lim = [max(np.max(np.abs(d)) for d in diffs[col]) for col in range(3)]

    fig = Figure(figsize=(30, 13), constrained_layout=True)
    gs = gridspec.GridSpec(no_of_rows, 3, figure=fig, height_ratios=[1] * no_of_rows,
                           hspace=0.001, wspace=0.001)
    ax = np.array([[fig.add_subplot(gs[i, j]) for j in range(3)] for i in range(no_of_rows)])

    true_fields = (h0[0], u0[0], s0[0])
    titles = ("True thickness", "True velocity", "True smb")
    labels = ("Thickness (m)", "Velocity (m/yr)", "SMB (m/yr)")
    for col in range(3):
        ax[0, col].set_title(titles[col], fontsize=18, fontweight="bold", loc="center")
        im = ax[0, col].tripcolor(tri, true_fields[col], shading="gouraud", cmap="viridis",
                                  vmin=true_fields[col].min(), vmax=true_fields[col].max())
        cb = fig.colorbar(im, ax=ax[0, col], fraction=0.025, pad=0.02)
        cb.set_label(labels[col], fontsize=18, fontweight="bold")
        cb.ax.tick_params(labelsize=16, width=1.4, direction="out", length=6)
        if col == 1:
            cb.locator = MaxNLocator(nbins=5)
            cb.update_ticks()

    last = [None, None, None]
    for col in range(3):
        ax[1, col].set_title("no assimilation", fontsize=18, fontweight="bold", loc="center")
        last[col] = ax[1, col].tripcolor(tri, diffs[col][0], shading="gouraud", cmap="bwr",
                                         vmin=-lim[col], vmax=lim[col])
    for ii, (year, d) in enumerate(zip(desired_years, per_step)):
        for col in range(3):
            ax[ii + 2, col].set_title(f"after {year:g} years of assimilation", fontsize=18,
                                      fontweight="bold", loc="center")
            last[col] = ax[ii + 2, col].tripcolor(tri, d[col], shading="gouraud", cmap="bwr",
                                                  vmin=-lim[col], vmax=lim[col])

    diff_labels = (r"$\Delta$ Thickness (m)", r"$\Delta |u|$ (m/yr)", r"$\Delta$ SMB (m/yr)")
    for col in range(3):
        cb = fig.colorbar(last[col], ax=ax[1:, col].ravel().tolist(), fraction=0.018, pad=0.02)
        cb.set_label(diff_labels[col], fontsize=18, fontweight="bold")
        cb.ax.tick_params(labelsize=14, width=1.4, length=6)
        cb.locator = MaxNLocator(nbins=6, symmetric=True)
        cb.update_ticks()

    for i in range(no_of_rows):
        for j in range(3):
            a = ax[i, j]
            for sp in a.spines.values():
                sp.set_linewidth(1.4)
            a.plot(x_h, y_h, "o", color="magenta", markeredgecolor="black", markersize=10)
            a.tick_params(axis="both", labelsize=17, width=1.2)
            if i == no_of_rows - 1:
                a.set_xlabel("X (km)", fontsize=18, fontweight="bold")
                a.xaxis.set_major_formatter(_km_formatter())
            else:
                a.set_xticklabels([])
            if j == 0:
                a.set_ylabel("Y (km)", fontsize=18, fontweight="bold")
                a.yaxis.set_major_formatter(_km_formatter())
            else:
                a.set_yticklabels([])
            for lab in (*a.get_xticklabels(), *a.get_yticklabels()):
                lab.set_fontweight("bold")
            a.text(0.019, 0.81, f"({string.ascii_lowercase[i]}$_{{{j + 1}}}$)", transform=a.transAxes,
                   fontsize=18, fontweight="bold",
                   bbox=dict(facecolor="white", edgecolor="black", boxstyle="round,pad=0.2"))
    notes = list(o.warnings) + [
        "Years → time steps: " + ", ".join(f"{y:g}→{s} (t={t[s]:g})" for y, s in zip(desired_years, steps)),
        "Each column shares one symmetric colour scale for all difference rows.",
    ]
    return ViewResult(figure=fig, notes=notes)


def _km_formatter():
    from matplotlib.ticker import FuncFormatter

    # the notebook relabels the metre ticks as str(float(tick / 1000))
    return FuncFormatter(lambda v, _pos: f"{float(v / 1000)}")


def transects(ctx: ViewContext, opts: dict) -> ViewResult:
    """cell 13: bed, ice base and surface along the y = 0 centreline for the
    true, wrong and assimilated (EnKF) runs at one step, and the thickness."""
    o = ctx.outputs
    hdim = _hdim(ctx)
    xy, _tri = _coords(ctx)
    Lx, Ly, nx, ny = _grid(ctx)
    meta = o.metadata
    if "b_io" not in meta:
        raise ResultsError("missing", "true-wrong-icepack.h5 has no b_io (bed elevation at inflow/outflow)")
    b_in, b_out = float(meta["b_io"][0]), float(meta["b_io"][1])
    step = opts["step"]
    xs = np.array([(Lx * k / nx, 0) for k in range(nx + 1)], dtype=float)
    nodes = []
    for x, y in xs:
        d = np.hypot(xy[:, 0] - x, xy[:, 1] - y)
        k = int(np.argmin(d))
        if d[k] > 1e-6 * max(Lx, Ly):
            raise ResultsError("incomplete", f"no mesh node at transect point ({x:g}, {y:g})")
        nodes.append(k)
    nodes = np.array(nodes)
    bs = b_in - (b_in - b_out) * xs[:, 0] / Lx

    def surf(h):
        return np.maximum(h + bs, (1 - RHO_I_OVER_RHO_W) * h)

    h = [o.true_state[:hdim, step][nodes], o.nurged_state[:hdim, step][nodes], o.ensemble_mean[:hdim, step][nodes]]
    s = [surf(v) for v in h]
    fig = Figure(figsize=(8, 7))
    axes = fig.subplots(2, 1)
    color_bed, color_ice_base, color_surface = "#3C3C3C", "#1B9E77", "#D95F02"
    xk = xs[:, 0] / 1e3
    for (name, ls), hv, sv in zip((("True", "-"), ("Wrong", ":"), ("EnKF", "--")), h, s):
        axes[0].plot(xk, bs, color=color_bed, linestyle=ls, label=f"Bed ({name})")
        axes[0].plot(xk, sv - hv, color=color_ice_base, linestyle=ls, label=f"Ice base ({name})")
        axes[0].plot(xk, sv, color=color_surface, linestyle=ls, label=f"Surface ({name})")
    axes[0].set_ylabel("Elevation (m)")
    axes[0].set_title("Elevation Transects")
    axes[0].legend(loc="upper left", bbox_to_anchor=(1, 1), fontsize="small")
    axes[0].grid(True, linestyle="--", alpha=0.5)
    for (name, ls), hv in zip((("True", "-"), ("Wrong", ":"), ("EnKF", "--")), h):
        axes[1].plot(xk, hv, color=color_ice_base, linestyle=ls, label=f"Thickness ({name})")
    axes[1].set_xlabel("Distance along centerline (km)")
    axes[1].set_ylabel("Thickness (m)")
    axes[1].set_title("Thickness Transects")
    axes[1].legend()
    axes[1].grid(True, linestyle="--", alpha=0.5)
    fig.tight_layout()
    notes = list(o.warnings) + [
        f"Time step {step} (t = {o.t[step]:g} years). Surface = max(h + b, (1 − ρI/ρW)·h), "
        "as icepack.compute_surface."]
    return ViewResult(figure=fig, notes=notes)


def smb_along_x(ctx: ViewContext, opts: dict) -> ViewResult:
    """cell 17: true, initial (wrong) and assimilated SMB at every node,
    ordered by x."""
    o = ctx.outputs
    hdim = _hdim(ctx)
    xy, _tri = _coords(ctx)
    step = opts["step"]
    order = np.argsort(xy[:, 0])
    _h, _u, smb = _fields(o, hdim, step)
    fig = Figure(figsize=(9, 5))
    ax = fig.subplots()
    ax.plot(xy[order, 0] / 1000, smb[0][order], label="True SMB")
    ax.plot(xy[order, 0] / 1000, smb[1][order], label="Initial/wrong SMB")
    ax.plot(xy[order, 0] / 1000, smb[2][order], label=f"Assimilated SMB, t={o.t[step]:.1f} yr")
    ax.set_xlabel("X (km)")
    ax.set_ylabel("SMB (m/yr)")
    ax.legend()
    ax.grid(True)
    return ViewResult(figure=fig, notes=list(o.warnings))


# ─────────────────────────────────────────────────────────────────────────────
# cells 21 & 22: error tables
# ─────────────────────────────────────────────────────────────────────────────
def error_tables(ctx: ViewContext, opts: dict) -> ViewResult:
    o = ctx.outputs
    hdim = _hdim(ctx)
    t = o.t[:_n(ctx)]
    years = list(opts["years"])
    steps = dx.nearest_steps(t, years)
    rows = []
    for var, b in (("h", 0), ("u", 1), ("v", 2), ("smb", 3)):
        idx = dx.block(b, hdim)
        for year, step in zip(years, steps):
            st = dx.error_stats(o.true_state[idx, step], o.ensemble_mean[idx, step])
            rows.append([var, f"{year:10.1f}".strip(), f"{t[step]:10.1f}".strip(),
                         f"{st['rmse']:.6f}", f"{st['max_error']:.6f}", f"{st['mean_error']:.6f}"])
    table = {
        "caption": f"{o.filter_type} mean minus truth at the requested times (nearest model step)",
        "columns": ["Variable", "Requested", "Actual", "RMSE", "Max error", "Mean error"],
        "rows": rows,
    }
    return ViewResult(table=table, notes=list(o.warnings))


_H_INDX = Option("h_indx", "Grid point (state index)", "int", default=242, minimum=0, maximum=_max_index,
                 help="Index of the node within one field (the notebook uses 242).")
_STEP_LAST = Option("step", "Time step", "int", default=lambda ctx: max(0, _max_step(ctx) - 1),
                    minimum=0, maximum=_max_step, help="Model time-step index (the notebook uses nt − 2).")
_YEARS = Option("years", "Years", "float_list", default=DESIRED_YEARS, minimum=0,
                help="Times (years) shown as rows; the nearest model step is used.")

register(ViewSpec(
    id="icepack.ensemble_profiles", title="Thickness, velocity and SMB at a grid point",
    description="Ensemble members, ensemble mean, truth, the wrong (no-assimilation) run and "
                "the observations over time at one node.",
    source="synthetic_ice_stream/read_results.ipynb, cell 14 (husmb_profile_middle_s.png)",
    build=ensemble_profiles, needs=("ensemble_members",),
    options=(_H_INDX, Option("xlim", "Time axis (years)", "float_list", default=[-1, 125],
                             help="The notebook fixes the axis to −1…125 years.")),
))
register(ViewSpec(
    id="icepack.truth_vs_wrong", title="Truth vs. wrong run at the domain midpoint",
    description="The true and wrong (perturbed) trajectories of thickness (with observations) or "
                "SMB at the notebook's middle profile point.",
    source="synthetic_ice_stream/read_results.ipynb, cells 5–6",
    build=truth_vs_wrong,
    options=(Option("variable", "Variable", "choice", default="thickness", choices=("thickness", "smb")),
             Option("h_indx", "Grid point (state index)", "int", default=_middle_index, minimum=0,
                    maximum=_max_index, help="Defaults to the notebook's middle profile point.")),
))
register(ViewSpec(
    id="icepack.snapshot_maps", title="Thickness and velocity maps",
    description="True, wrong and assimilated thickness and speed at one time step, on common colour scales.",
    source="synthetic_ice_stream/read_results.ipynb, cell 10 (hu_snapshot_middle_s.png)",
    build=snapshot_maps, needs=("dof_coordinates",), options=(_STEP_LAST,),
))
register(ViewSpec(
    id="icepack.difference_maps", title="Errors before and after assimilation",
    description="True fields, then the signed error of the wrong run and of the assimilation "
                "after selected years, per field.",
    source="synthetic_ice_stream/read_results.ipynb, cell 16 (husmb_true_diff.png)",
    build=difference_maps, needs=("dof_coordinates",), options=(_YEARS,),
))
register(ViewSpec(
    id="icepack.transects", title="Centreline transects",
    description="Bed, ice base and surface elevation and thickness along y = 0 for the true, "
                "wrong and assimilated runs.",
    source="synthetic_ice_stream/read_results.ipynb, cell 13",
    build=transects, needs=("dof_coordinates",),
    options=(Option("step", "Time step", "int", default=lambda ctx: min(65, _max_step(ctx)), minimum=0,
                    maximum=_max_step,
                    help="The notebook draws this after its step-65 maps; any step can be chosen."),),
))
register(ViewSpec(
    id="icepack.smb_along_x", title="SMB along the flow direction",
    description="True, initial and assimilated surface mass balance at every node, ordered by x.",
    source="synthetic_ice_stream/read_results.ipynb, cell 17",
    build=smb_along_x, needs=("dof_coordinates",),
    options=(Option("step", "Time step", "int",
                    default=lambda ctx: dx.nearest_steps(ctx.outputs.t[:_n(ctx)], [DESIRED_YEARS[-1]])[0],
                    minimum=0, maximum=_max_step,
                    help="The notebook uses the last difference-map year (85)."),),
))
register(ViewSpec(
    id="icepack.error_tables", title="Error statistics", kind="table",
    description="RMSE, maximum and mean error of the ensemble mean against the truth for "
                "h, u, v and SMB at selected years.",
    source="synthetic_ice_stream/read_results.ipynb, cells 21–22",
    build=error_tables, options=(_YEARS,),
))
