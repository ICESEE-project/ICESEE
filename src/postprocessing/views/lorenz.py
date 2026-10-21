"""Lorenz example views (applications/lorenz_model/examples/lorenz96/read_results.ipynb)."""
from __future__ import annotations

from matplotlib.figure import Figure

from . import ViewContext, ViewResult, ViewSpec, register

#: the notebook sets ``mpl.rc("font", weight="bold", size=14)`` and uses
#: mathtext when LaTeX is unavailable; mathtext is used unconditionally here
#: so the figure does not depend on the server's TeX installation.
LORENZ_RC = {
    "text.usetex": False,
    "mathtext.fontset": "dejavusans",
    "font.family": "DejaVu Sans",
    "font.weight": "bold",
    "font.size": 14,
}


def state_estimation(ctx: ViewContext, opts: dict) -> ViewResult:
    """read_results.ipynb cell 5: true, background, observation and analysis
    trajectories of x, y, z, with the observation window shaded."""
    o = ctx.outputs
    t, ind_m, tm_m, w = o.t, o.obs_index, o.obs_max_time, o.observations
    n = min(t.size, o.ensemble_mean.shape[1], o.true_state.shape[1])
    notes = list(o.warnings)

    fig = Figure(figsize=(10, 8))
    ax = fig.subplots(nrows=3, ncols=1).flat
    for k in range(3):
        ax[k].plot(t[:n], o.true_state[k, :n], label="True", linewidth=3)
        ax[k].plot(t[:n], o.nurged_state[k, :n], ":", label="Background", linewidth=3)
        ax[k].plot(t[ind_m], w[k, :], "o", fillstyle="none",
                   label="Observation", markersize=8, markeredgewidth=2)
        ax[k].plot(t[:n], o.ensemble_mean[k, :n], "--", label="Analysis", linewidth=3)
        ax[k].set_xlabel(r"$t$", fontsize=22)
        ax[k].axvspan(0, tm_m, alpha=0.4, lw=0)

    ax[0].legend(loc="center", bbox_to_anchor=(0.5, 1.25), ncol=4, fontsize=15)
    ax[0].set_ylabel(r"$x(t)$")
    ax[1].set_ylabel(r"$y(t)$")
    ax[2].set_ylabel(r"$z(t)$")
    fig.subplots_adjust(hspace=0.5)
    notes.append(f"Shaded: observation window 0 – {tm_m:g}; "
                 f"{len(ind_m)} observation times; filter {o.filter_type}.")
    return ViewResult(figure=fig, notes=notes)


register(ViewSpec(
    id="lorenz.state_estimation",
    title="State estimation: x, y, z",
    description="True and background trajectories, observations and the ensemble-mean "
                "analysis for each state variable; the observation window is shaded.",
    source="lorenz96/read_results.ipynb, cell 5 (figure saved as figures/results.png)",
    build=state_estimation,
    rc=LORENZ_RC,
))
