"""Supplementary diagnostics available for every ICESEE run.

These are NOT part of any example notebook and are labelled as such. They
summarise the assimilation skill (analysis vs. truth, background vs. truth)
and the ensemble spread per state variable, which the notebooks show only
implicitly through the trajectories.
"""
from __future__ import annotations

import numpy as np
from matplotlib.figure import Figure

from .. import diagnostics as dx
from . import ViewContext, ViewResult, ViewSpec, register


def _variables(ctx: ViewContext) -> list[str]:
    names = ctx.section("enkf-parameters").get("vec_inputs")
    if isinstance(names, (list, tuple)) and names:
        return [str(n) for n in names]
    return []


def error_and_spread(ctx: ViewContext, opts: dict) -> ViewResult:
    o = ctx.outputs
    names = _variables(ctx)
    nd = o.n_state
    nblocks = len(names) if names and nd % len(names) == 0 else 1
    size = nd // nblocks
    names = names if nblocks > 1 else ["state"]
    n = min(o.t.size, o.ensemble_mean.shape[1], o.true_state.shape[1])
    t = o.t[:n]

    spread = None
    if o.ensemble_file is not None:
        spread = np.zeros((nblocks, n))
        with o.open_ensemble() as ds:
            chunk = max(1, int(2e7 // max(1, ds.shape[0] * ds.shape[1])))
            for start in range(0, n, chunk):
                stop = min(n, start + chunk)
                std = np.asarray(ds[:, :, start:stop]).std(axis=1, ddof=1) if ds.shape[1] > 1 \
                    else np.zeros((ds.shape[0], stop - start))
                for b in range(nblocks):
                    spread[b, start:stop] = np.sqrt(np.mean(std[dx.block(b, size)] ** 2, axis=0))

    fig = Figure(figsize=(10, 2.6 * nblocks + 0.6))
    axes = np.atleast_1d(fig.subplots(nrows=nblocks, ncols=1, sharex=True))
    rows = []
    for b, ax in enumerate(axes):
        rows_b = dx.block(b, size)
        truth = o.true_state[rows_b, :n]
        an = dx.rmse_over_time(truth, o.ensemble_mean[rows_b, :n])
        bg = dx.rmse_over_time(truth, o.nurged_state[rows_b, :n])
        ax.plot(t, bg, color="tab:green", lw=1.6, label="Background RMSE")
        ax.plot(t, an, color="tab:blue", lw=1.6, label=f"{o.filter_type} analysis RMSE")
        if spread is not None:
            ax.plot(t, spread[b], color="tab:gray", lw=1.2, ls="--", label="Ensemble spread")
        if o.obs_index.size:
            ax.plot(t[o.obs_index[o.obs_index < n]], an[o.obs_index[o.obs_index < n]], "kx",
                    ms=6, label="Observation times")
        ax.set_ylabel(f"{names[b]}")
        ax.grid(True, ls="--", alpha=0.5)
        rows.append([names[b], f"{bg[-1]:.4g}", f"{an[-1]:.4g}",
                     f"{float(np.mean(an)):.4g}", "" if spread is None else f"{spread[b][-1]:.4g}"])
    axes[0].legend(loc="upper right", fontsize=9, ncol=2)
    axes[-1].set_xlabel("Time")
    fig.suptitle("Assimilation error and ensemble spread (supplementary)", fontsize=12)
    table = {
        "caption": "Final-time RMSE against the true state, and time-averaged analysis RMSE",
        "columns": ["Variable", "Background RMSE (final)", "Analysis RMSE (final)",
                    "Analysis RMSE (mean)", "Spread (final)"],
        "rows": rows,
    }
    notes = list(o.warnings) + ["Supplementary diagnostic: not part of the example's "
                                "read_results.ipynb. RMSE is taken over the state rows "
                                "of each variable."]
    return ViewResult(figure=fig, table=table, notes=notes)


register(ViewSpec(
    id="da.error_and_spread",
    title="Assimilation error and spread",
    description="Root-mean-square error of the analysis and of the background against the "
                "true state, with the ensemble spread, per state variable.",
    source="Supplementary (computed by ICESEE post-processing; not in read_results.ipynb)",
    build=error_and_spread,
    supplementary=True,
))
