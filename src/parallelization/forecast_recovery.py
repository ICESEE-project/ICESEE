"""Transactional subcycling for adapters that can restore all forecast state."""
from __future__ import annotations

from dataclasses import dataclass, asdict
import math
import random
import numpy as np


class RecoverableForecastError(RuntimeError):
    """A model adapter explicitly identified a numerical step rejection."""


class ForecastRecoveryExhausted(RuntimeError):
    """The original state was restored after all permitted attempts failed."""


@dataclass(frozen=True)
class ForecastRecoveryPolicy:
    max_halvings: int
    min_dt: float
    max_substeps: int = 256

    def __post_init__(self):
        if (isinstance(self.max_halvings, bool) or not isinstance(self.max_halvings, int)
                or not 0 <= self.max_halvings <= 20):
            raise ValueError("forecast recovery max_halvings must be an integer in [0, 20]")
        if not math.isfinite(self.min_dt) or self.min_dt <= 0:
            raise ValueError("forecast recovery min_dt must be finite and positive")
        if (isinstance(self.max_substeps, bool) or not isinstance(self.max_substeps, int)
                or not 1 <= self.max_substeps <= 2**20):
            raise ValueError("forecast recovery max_substeps must be a positive bounded integer")

    @classmethod
    def from_options(cls, options):
        raw = options.get("forecast_recovery")
        if raw is None or raw is False:
            return None
        if not isinstance(raw, dict) or set(raw) - {"max_halvings", "min_dt", "max_substeps"}:
            raise ValueError("forecast_recovery requires max_halvings/min_dt and optional max_substeps")
        return cls(**raw)

    def describe(self):
        return {"method": "transactional-dyadic-v1", **asdict(self)}


def _errors(comm, error):
    return [e for e in (comm.allgather(error) if comm is not None else [error]) if e]


def advance_with_recovery(*, start_time, interval, policy, capture, restore,
                          advance, validate, comm=None, stable_dt=None, report=None):
    """Advance exactly one external interval, with bounded internal trials.

    ``capture``/``restore`` include auxiliary model state, halos and caches;
    ``advance(time, dt)`` and ``validate`` must not write output or draw noise.
    Each callback must return or raise coherently from its own model collectives.
    Process loss and a model stuck inside MPI require scheduler-level recovery.
    A local CFL bound may be supplied in model time units. All ranks use its
    minimum. Python/NumPy legacy RNG state is restored; adapter-owned generators
    belong in the captured state. A successful trial commits no files itself.
    """
    start_time, interval = float(start_time), float(interval)
    if not math.isfinite(start_time) or not math.isfinite(interval) or interval <= 0:
        raise ValueError("forecast interval must be finite and positive")
    if interval < policy.min_dt and not math.isclose(interval, policy.min_dt, rel_tol=1e-12):
        raise ValueError("forecast interval is smaller than recovery min_dt")
    signature = (start_time, interval, tuple(sorted(policy.describe().items())), stable_dt is not None)
    if comm is not None and any(value != signature for value in comm.allgather(signature)):
        raise ValueError("spatial ranks disagree on forecast recovery interval or policy")
    counts = None
    if stable_dt is not None:
        bounds = comm.allgather(float(stable_dt)) if comm is not None else [float(stable_dt)]
        if any(not math.isfinite(v) or v <= 0 for v in bounds):
            raise ValueError("model stability bounds must be finite and positive")
        counts = max(1, math.ceil(interval / min(bounds)))
    first_level = 0 if counts is None else (counts - 1).bit_length()
    snapshot, problem = None, None
    try:
        snapshot = capture()
    except Exception as exc:
        problem = f"{type(exc).__name__}: {exc}"
    errors = _errors(comm, problem)
    if errors:
        raise RuntimeError("forecast capture failed: " + "; ".join(errors))
    numpy_rng, python_rng = np.random.get_state(), random.getstate()

    def rollback():
        # A restore failure is fatal; no rank starts another model trial.
        problem = None
        try:
            restore(snapshot)
            np.random.set_state(numpy_rng)
            random.setstate(python_rng)
        except Exception as exc:
            problem = f"{type(exc).__name__}: {exc}"
        errors = _errors(comm, problem)
        if errors:
            raise RuntimeError("forecast rollback failed: " + "; ".join(errors))

    last_error = "no subdivision satisfies the configured limits"
    attempts = 0
    for level in range(first_level, policy.max_halvings + 1):
        subdivisions = 2**level
        dt = interval / subdivisions
        if subdivisions > policy.max_substeps or (dt < policy.min_dt and not math.isclose(dt, policy.min_dt, rel_tol=1e-12)):
            break
        rollback()
        attempts += 1
        accepted = True
        for j in range(subdivisions):
            error = None
            try:
                # Recompute from the external grid rather than accumulating dt.
                advance(start_time + interval * j / subdivisions, dt)
                validate()
            except RecoverableForecastError as exc:
                error = ("retry", f"{type(exc).__name__}: {exc}")
            except Exception as exc:
                error = ("fatal", f"{type(exc).__name__}: {exc}")
            errors = _errors(comm, error)
            if errors:
                rollback()
                last_error = "; ".join(e[1] for e in errors)
                if any(e[0] == "fatal" for e in errors):
                    raise RuntimeError("non-retryable forecast failure: " + last_error)
                accepted = False
                if report is not None:
                    report(dict(status="rejected", attempts=attempts, substeps=subdivisions,
                                failed_substep=j, dt=dt, reason=last_error))
                break
        if accepted:
            result = dict(status="accepted", attempts=attempts, substeps=subdivisions,
                          dt=dt, start_time=start_time, end_time=start_time + interval)
            if report is not None:
                report(result)
            return result
    rollback()
    raise ForecastRecoveryExhausted(f"forecast recovery exhausted after {attempts} attempts: {last_error}")
