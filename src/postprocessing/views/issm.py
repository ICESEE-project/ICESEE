"""ISSM example views
(applications/issm_model/examples/ISMIP_Choi/read_results.ipynb).

State layout (as in the notebook): six blocks of ``hdim = nd // 6`` rows --
Thickness, Surface, Vx, Vy, bed, friction coefficient.

The notebook rebuilds the time axis and the observation schedule from the
run's ``params.yaml`` (cell 6) rather than using the stored ones; this view
does the same, with ICESEE's own ``UtilsFunctions.generate_observation_schedule``.
"""
from __future__ import annotations

import numpy as np
from matplotlib.figure import Figure

from ..outputs import ResultsError
from . import Option, ViewContext, ViewResult, ViewSpec, register

#: (block, y label, observation handling) per panel, top to bottom (cell 24)
PANELS = (
    (0, "Thickness (m)", None),
    (1, "Surface (m)", "plain"),
    ((2, 3), "Velocity (m/year)", "magnitude"),
    (4, "Bed topography", "mask_zero"),
    (5, "Friction", None),
)


def _schedule(ctx: ViewContext):
    """cell 6: ``t`` and ``ind_m`` from params.yaml."""
    modeling = ctx.section("modeling-parameters")
    enkf = ctx.section("enkf-parameters")
    try:
        tinitial = float(modeling.get("tinitial"))
        num_years = float(modeling.get("num_years"))
        tpy = float(modeling.get("timesteps_per_year"))
    except (TypeError, ValueError):
        raise ResultsError("missing", "params.yaml needs modeling-parameters tinitial, num_years and "
                           "timesteps_per_year to rebuild the time axis") from None
    t = np.linspace(tinitial, num_years, int((num_years - tinitial) / tpy) + 1)
    from ICESEE.src.utils.utils import UtilsFunctions

    _obs_t, ind_m, _n = UtilsFunctions(enkf).generate_observation_schedule(t=t)
    return t, np.asarray(ind_m, dtype=int)


def _style_axis(ax, ylabel=None, xlabel=None):
    """style_axis() of the notebook."""
    if ylabel is not None:
        ax.set_ylabel(ylabel, fontsize=16, fontweight="bold")
    if xlabel is not None:
        ax.set_xlabel(xlabel, fontsize=16, fontweight="bold")
    ax.grid(True, linestyle="--", alpha=0.6)
    ax.tick_params(labelsize=14)
    xticks, yticks = ax.get_xticks(), ax.get_yticks()
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    ax.set_xticklabels([f"{int(x)}" for x in xticks], fontweight="bold")
    ax.set_yticklabels([f"{int(y)}" for y in yticks], fontweight="bold")


def profile_time_series(ctx: ViewContext, opts: dict) -> ViewResult:
    """cell 24: five stacked time series at one node -- ensemble members,
    truth, wrong run, ensemble mean and (where observed) observations."""
    o = ctx.outputs
    t, ind_m = _schedule(ctx)
    w = o.observations[:, :ind_m.shape[0]]
    hdim = o.n_state // 6
    h_indx = opts["h_indx"]
    filter_type = o.filter_type
    _nt = min(499, t.size)
    notes = list(o.warnings)
    avail = min(o.true_state.shape[1], o.ensemble_mean.shape[1])
    if _nt > avail:
        notes.append(f"the outputs cover {avail} of the {_nt} time steps implied by params.yaml")
        _nt = avail

    fig = Figure(figsize=(12, 24))
    axes = fig.subplots(5, 1, sharex=True)
    members = o.ensemble_file is not None
    for ax, (blk, ylabel, obs_mode) in zip(axes, PANELS):
        if obs_mode == "magnitude":
            iu, iv = blk[0] * hdim + h_indx, blk[1] * hdim + h_indx
            y_true = np.sqrt(o.true_state[iu, :] ** 2 + o.true_state[iv, :] ** 2)
            y_nurged = np.sqrt(o.nurged_state[iu, :] ** 2 + o.nurged_state[iv, :] ** 2)
            y_mean = np.sqrt(o.ensemble_mean[iu, :] ** 2 + o.ensemble_mean[iv, :] ** 2)
            y_mem = np.sqrt(o.ensemble_members(iu) ** 2 + o.ensemble_members(iv) ** 2) if members else None
            t_obs, y_obs = t[ind_m], np.sqrt(w[iu, :] ** 2 + w[iv, :] ** 2)
        else:
            idx = blk * hdim + h_indx
            y_true, y_nurged, y_mean = o.true_state[idx, :], o.nurged_state[idx, :], o.ensemble_mean[idx, :]
            y_mem = o.ensemble_members(idx) if members else None
            t_obs = y_obs = None
            if obs_mode is not None:
                vals = w[idx, :]
                if obs_mode == "mask_zero":
                    mask = vals != 0
                    t_obs, y_obs = t[ind_m][mask], vals[mask]
                else:
                    t_obs, y_obs = t[ind_m], vals
        if y_mem is not None:
            for ens in range(y_mem.shape[1]):
                ax.plot(t[:_nt], y_mem[:_nt, ens], lw=0.5, color="gray",
                        label=f"{filter_type}" if ens == 0 else "")
        ax.plot(t[:_nt], y_true[:_nt], "--r", label="True", linewidth=2)
        ax.plot(t[:_nt], y_nurged[:_nt], "--g", label="Wrong", linewidth=2)
        ax.plot(t[:_nt], y_mean[:_nt], "--b", label=f"{filter_type} mean", linewidth=2)
        if t_obs is not None:
            ax.plot(t_obs[1:], y_obs[1:], "kx", label="Observations", markersize=8)
        ax.legend(prop={"size": 12, "weight": "bold"})
        _style_axis(ax, ylabel=ylabel)
    _style_axis(axes[-1], ylabel="Friction", xlabel="Time (years)")
    fig.tight_layout()
    notes.append(f"Node h_indx = {h_indx}; time axis and observation times rebuilt from params.yaml "
                 f"({ind_m.size} observation times).")
    return ViewResult(figure=fig, notes=notes)


register(ViewSpec(
    id="issm.profile_time_series", title="Thickness, surface, velocity, bed and friction at a node",
    description="Ensemble members, ensemble mean, truth, the wrong run and observations over time "
                "for each ISSM state field at one mesh node.",
    source="ISMIP_Choi/read_results.ipynb, cell 24 (profile_middle_time_series.png)",
    build=profile_time_series, needs=("params", "ensemble_members"),
    options=(Option("h_indx", "Mesh node (state index)", "int", default=390, minimum=0,
                    maximum=lambda ctx: ctx.outputs.n_state // 6 - 1,
                    help="Index of the node within one field (the notebook uses 390)."),),
))
