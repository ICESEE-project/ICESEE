"""Flowline example views (applications/flowline_model/examples/flowline_1d).

The example's read_results.ipynb has no working flowline figure (cell 5 is
the Lorenz x/y/z plot, cell 9 uses variables the notebook never defines).
:func:`state_estimation` draws what cell 9 sets out to: grounding-line
position, mid-profile thickness and terminus velocity over time, scaled to
physical units the way run_da_flowline.py derives ``hscale``/``uscale``/
``xscale``.
"""
from __future__ import annotations

import numpy as np
from matplotlib.figure import Figure

from ..outputs import ResultsError
from . import ViewContext, ViewResult, ViewSpec, register


def _scales(ctx: ViewContext) -> dict:
    phys, mod = ctx.section("physical-parameters"), ctx.section("modeling-parameters")
    try:
        hscale, year = float(phys["hscale"]), float(mod["year"])
        uscale = float(phys["rho_ice"]) * float(phys["g"]) * hscale / year / float(phys["C"])
        nx = int(phys["N1"]) + int(phys["N2"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ResultsError("invalid", f"params.yaml lacks the flowline scales ({exc})") from exc
    return {"h": hscale, "u": uscale * year, "xg": uscale * hscale * year / 1e3, "nx": nx}


def state_estimation(ctx: ViewContext, opts: dict) -> ViewResult:
    """Grounding line (km), thickness at the mid-profile node (m) and velocity
    at the terminus (m/yr): ensemble members, truth, background (wrong),
    ensemble-mean analysis and observations."""
    o = ctx.outputs
    sc = _scales(ctx)
    nx = sc["nx"]
    if o.n_state != 2 * nx + 1:
        raise ResultsError("incomplete", f"the state has {o.n_state} rows; a flowline state with "
                           f"N1 + N2 = {nx} nodes has {2 * nx + 1}")
    n = min(o.t.size, o.ensemble_mean.shape[1], o.true_state.shape[1])
    t = o.t[:n]
    panels = [
        (2 * nx, sc["xg"], "Grounding line (km)"),
        (nx // 2 - 1, sc["h"], f"h at node {nx // 2} (m)"),
        (2 * nx - 1, sc["u"], "u at terminus (m/yr)"),
    ]
    notes = list(o.warnings)

    fig = Figure(figsize=(9, 10))
    axes = fig.subplots(nrows=len(panels), ncols=1, sharex=True)
    for ax, (row, scale, label) in zip(axes, panels):
        if o.ensemble_file is not None:
            members = o.ensemble_members(row)[:n] * scale
            for i in range(members.shape[1]):
                ax.plot(t, members[:, i], lw=0.5, color="gray", label="Ensemble" if i == 0 else None)
        ax.plot(t, o.true_state[row, :n] * scale, lw=3, color="black", label="Truth")
        ax.plot(t, o.nurged_state[row, :n] * scale, lw=3, color="red", label="Background")
        ax.plot(t, o.ensemble_mean[row, :n] * scale, lw=3, color="blue", label=f"{o.filter_type} mean")
        ax.scatter(o.t[o.obs_index], o.observations[row, :] * scale, color="green", s=60,
                   edgecolor="black", zorder=3, label="Observations")
        ax.axvspan(0, o.obs_max_time, alpha=0.12, lw=0)
        ax.set_ylabel(label)
    axes[0].legend(loc="center", bbox_to_anchor=(0.5, 1.18), ncol=5, fontsize=10)
    axes[-1].set_xlabel("time (years)")
    fig.tight_layout()
    notes.append(f"Shaded: observation window 0 – {o.obs_max_time:g} years; "
                 f"{len(o.obs_index)} observation times; filter {o.filter_type}.")
    return ViewResult(figure=fig, notes=notes)


register(ViewSpec(
    id="flowline.state_estimation",
    title="Grounding line, thickness and velocity",
    description="Grounding-line position, mid-profile ice thickness and terminus velocity "
                "over time: ensemble members, truth, background, ensemble-mean analysis "
                "and observations, in physical units.",
    source="flowline_1d/read_results.ipynb, cell 9 (its intended figure; the cell itself "
           "does not run)",
    build=state_estimation,
    needs=("params", "ensemble_members"),
    supplementary=True,
))
