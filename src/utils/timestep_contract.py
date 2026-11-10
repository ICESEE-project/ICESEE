"""Map forecast cycles to saved state columns and observation columns."""
from __future__ import annotations

import numpy as np

TIME_GRID_CONTRACT = "state-columns-v2"


def observation_cycle_steps(indices, nt=None):
    """Observation columns j>0 are analyzed after forecast cycle j-1.

    Initial-time assimilation needs a separate pre-forecast analysis stage;
    refusing it prevents applying that observation after an unwanted advance.
    """
    raw = np.asarray(indices)
    if raw.ndim != 1 or raw.dtype.kind not in "iu" and raw.size:
        raise ValueError("observation columns must be a one-dimensional integer array")
    columns = raw.astype(np.int64, copy=False)
    if columns.size and columns[0] == 0:
        raise ValueError("initial-time observations require a pre-forecast analysis stage; choose observation times after t[0]")
    if np.any(columns < 1) or np.any(np.diff(columns) <= 0):
        raise ValueError("observation columns must be positive and strictly increasing")
    if nt is not None and np.any(columns > int(nt)):
        raise ValueError("observation column exceeds the final model state")
    return columns - 1


def observation_is_due(cycle, observation_column, indices):
    columns = np.asarray(indices)
    column = int(observation_column)
    if columns.ndim != 1 or columns.size and columns.dtype.kind not in "iu":
        raise ValueError("observation columns must be a one-dimensional integer array")
    if columns.size and columns[0] <= 0:
        raise ValueError("initial-time observations require a pre-forecast analysis stage; choose observation times after t[0]")
    return bool(0 <= column < len(columns) and int(cycle) + 1 == columns[column])


def completed_observations(indices, state_column):
    observation_cycle_steps(indices)
    return int(np.searchsorted(np.asarray(indices, dtype=np.int64), int(state_column), side="right"))


def prepare_timestep_contract(comm, data_path, *, time_grid=None):
    """Record the clock before creating state; reject incompatible old state."""
    import hashlib
    import json
    import os
    from pathlib import Path
    import tempfile

    error = None
    if comm.Get_rank() == 0:
        try:
            root = Path(data_path)
            root.mkdir(parents=True, exist_ok=True)
            marker = root / "time_grid_contract.json"
            signature = None
            if time_grid is not None:
                times = np.asarray(time_grid, dtype="<f8")
                if times.ndim != 1 or times.size == 0 or not np.isfinite(times).all() or np.any(np.diff(times) <= 0):
                    raise ValueError("model time grid must be finite and strictly increasing")
                signature = hashlib.sha256(times.tobytes()).hexdigest()
            if marker.exists():
                recorded = json.loads(marker.read_text())
                if recorded.get("time_grid_contract") != TIME_GRID_CONTRACT or signature is not None and recorded.get("time_grid_signature") != signature:
                    raise ValueError("output uses an incompatible time-grid contract; preserve it and choose a fresh data_path")
            else:
                if ((root / "icesee_ensemble_data.h5").exists()
                        or any(root.glob("icesee_enkf_ens_*.h5"))
                        or (root / "_mode3_state_history").exists()):
                    raise ValueError("existing state has no current time-grid contract; preserve it and choose a fresh data_path")
                fd, name = tempfile.mkstemp(prefix=".time-grid-", suffix=".tmp", dir=root)
                try:
                    with os.fdopen(fd, "w") as stream:
                        json.dump({"time_grid_contract": TIME_GRID_CONTRACT, "time_grid_signature": signature}, stream)
                        stream.write("\n")
                        stream.flush()
                        os.fsync(stream.fileno())
                    os.replace(name, marker)
                    directory = os.open(root, os.O_RDONLY)
                    try:
                        os.fsync(directory)
                    finally:
                        os.close(directory)
                finally:
                    Path(name).unlink(missing_ok=True)
        except Exception as exc:
            error = f"{type(exc).__name__}: {exc}"
    error = comm.bcast(error, root=0)
    if error is not None:
        raise ValueError(f"cannot initialize time-grid contract: {error}")


def format_analysis_event(cycle, observation_column, kwargs):
    """Log the state and observation actually used by a completed analysis."""
    import json

    column = int(cycle) + 1
    slot = int(observation_column)
    observed_state = int(np.asarray(kwargs["obs_index"])[slot])
    if column != observed_state:
        raise ValueError("analysis state and observation time columns disagree")
    return "[ICESEE_ANALYSIS] " + json.dumps({
        "cycle": int(cycle), "state_column": column,
        "model_time": float(np.asarray(kwargs["t"])[column]),
        "observation_column": slot, "observation_state_column": observed_state,
    }, sort_keys=True)
