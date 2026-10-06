# ==============================================================================
# @des: Resolved run schedule (time grid and observation schedule) and CLI
#       override consistency checks. Model and execution-mode agnostic:
#       applications resolve their own time discretization (num_years, dt,
#       nt, t) and then call resolve_run_schedule() before any expensive
#       model work; the dispatcher calls verify_cli_overrides_respected() for
#       every run.
# ==============================================================================
from __future__ import annotations

from typing import Any, Mapping, MutableMapping, Optional

import numpy as np

from ICESEE.src.utils.performance import register_run_metadata


def resolve_run_schedule(
    icesee_kwargs: MutableMapping[str, Any],
    *,
    rank: int = 0,
    details: Optional[Mapping[str, Any]] = None,
) -> int:
    """Rebuild the synthetic observation schedule from the final time grid,
    print the resolved schedule once (on ``rank`` 0), and return the number
    of analysis events the forecast loop will perform.

    Expects ``t`` (physical times, length nt+1), ``nt``, ``dt`` and
    ``num_years`` to be the application's final, resolved values. The
    observation window (``obs_start_time``, ``obs_max_time``, ``freq_obs``)
    is taken as configured -- it is reported, never adjusted. When synthetic
    observations will be generated and the window yields no observation the
    forecast loop can reach (index < nt), a ValueError is raised before any
    model work, instead of running a data-assimilation experiment that
    assimilates nothing. ``details`` are extra lines for the summary.
    """
    from ICESEE.src.utils.utils import UtilsFunctions

    t = np.asarray(icesee_kwargs["t"], dtype=float)
    nt = int(icesee_kwargs["nt"])
    if t.size != nt + 1:
        raise ValueError(f"time grid has {t.size} points but nt={nt} requires {nt + 1}")

    if not icesee_kwargs.get("observations_available", False):
        obs_t, obs_index, count = UtilsFunctions(icesee_kwargs).generate_observation_schedule(
            **icesee_kwargs
        )
        icesee_kwargs.update(
            {"obs_t": obs_t, "obs_index": obs_index, "number_obs_instants": count, "m_obs": count}
        )
    obs_index = np.asarray(icesee_kwargs.get("obs_index", []), dtype=int)
    # An observation at the final state (step nt) is a snapshot of the
    # truth but is never assimilated: the forecast loop ends at step nt-1.
    analysis_steps = obs_index[obs_index < nt]
    analysis_events = int(analysis_steps.size)
    register_run_metadata(observation_snapshots=int(obs_index.size))

    if rank == 0:
        lines = [
            f"num_years       = {icesee_kwargs.get('num_years')}",
            f"dt              = {icesee_kwargs.get('dt')}",
            f"nt              = {nt}",
            f"len(t)          = {t.size}",
            f"time range      = [{t[0]:g}, {t[-1]:g}]",
            f"obs window      = [{icesee_kwargs.get('obs_start_time')}, "
            f"{icesee_kwargs.get('obs_max_time')}] every {icesee_kwargs.get('freq_obs')}",
            f"obs snapshots   = {obs_index.size}"
            + (f" (steps {obs_index.tolist()})" if obs_index.size <= 12 else ""),
            f"analysis events = {analysis_events}"
            + (f" (steps {analysis_steps.tolist()})" if analysis_events <= 12 else "")
            + ("; a snapshot at the final state (step nt) is not assimilated"
               if analysis_events < obs_index.size else ""),
        ]
        lines += [f"{key:<16}= {value}" for key, value in (details or {}).items()]
        print("[ICESEE] Resolved run schedule:\n  " + "\n  ".join(lines), flush=True)

    if analysis_events == 0 and icesee_kwargs.get("generate_synthetic_obs", True):
        raise ValueError(
            "The configured observation window "
            f"[{icesee_kwargs.get('obs_start_time')}, {icesee_kwargs.get('obs_max_time')}] "
            f"(every {icesee_kwargs.get('freq_obs')}) yields no analysis event within the "
            f"simulated interval [{t[0]:g}, {t[-1]:g}) of {nt} steps. Set --obs_start_time, "
            "--obs_max_time and/or --freq_obs inside the run window."
        )
    return analysis_events


def _same_value(requested: Any, actual: Any) -> bool:
    try:
        return bool(np.array_equal(np.asarray(requested), np.asarray(actual)))
    except Exception:
        return requested == actual


def verify_cli_overrides_respected(icesee_kwargs: Mapping[str, Any]) -> None:
    """Raise if any value requested on the command line was replaced by a
    different value before the run started (e.g. an application recomputing
    a derived quantity), so a CLI option is never silently ignored."""
    ignored = {
        key: (requested, icesee_kwargs.get(key))
        for key, requested in dict(icesee_kwargs.get("cli_overrides") or {}).items()
        if not _same_value(requested, icesee_kwargs.get(key))
    }
    if ignored:
        details = "; ".join(
            f"--{key}={requested!r} was replaced by {actual!r}"
            for key, (requested, actual) in ignored.items()
        )
        raise ValueError(
            "Command-line overrides were not honored by this application: "
            f"{details}. Use the option(s) this application derives these values from."
        )
