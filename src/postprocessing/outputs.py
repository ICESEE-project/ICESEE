"""Reading the data-assimilation outputs an ICESEE run writes.

This is the loading half of the examples' ``read_results.ipynb`` notebooks,
moved into a module so the notebooks, the launcher and tests all read a run
the same way. File names and dataset keys are exactly the ones the notebooks
use:

* ``<data_dir>/true-wrong-<model>.h5``    -- ``t``, ``obs_index``,
  ``obs_max_time``, ``run_mode`` (+ ``Lxy``, ``nxy``, ``b_io``)
* ``<data_dir>/true_nurged_states.h5``    -- ``true_state``, ``nurged_state``
* ``<data_dir>/synthetic_obs.h5``         -- ``hu_obs``, or mode 2's compact
  ``hu_obs_compact`` + ``obs_indices`` (scattered back to dense, as
  ``ICESEE.src.utils.tools.load_hu_obs_from_h5`` does)
* ``run_mode == 0``: ``<data_dir>/icesee_ensemble_data.h5`` -- ``ensemble``
  (state, member, time) and ``ensemble_mean`` (state, time)
* ``run_mode != 0``: ``<results_dir>/<filter>-<model>.h5`` --
  ``ensemble_vec_full``, ``ensemble_vec_mean``, ``ensemble_bg``
* execution mode 3: ``<data_dir>/_mode3_state_history/`` -- per-step
  checkpoints (``manifest.json`` + ``shard_*.h5``), read by
  :class:`Mode3History` as the same (state, member, time) ensemble

Nothing is interpreted here: arrays are returned as stored. Problems are
raised as :class:`ResultsError` with a ``kind`` a user interface can act on
(``missing`` / ``corrupt`` / ``incomplete``) instead of a bare ``KeyError``.
"""
from __future__ import annotations

import json
from contextlib import contextmanager
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

TRUE_WRONG = "true-wrong-{model}.h5"
STATES = "true_nurged_states.h5"
OBSERVATIONS = "synthetic_obs.h5"
ENSEMBLE = "icesee_ensemble_data.h5"
FILTER_OUTPUT = "{filter}-{model}.h5"
#: execution mode 3 (ICESEE.src.parallelization.mode3_checkpointing)
MODE3_HISTORY = "_mode3_state_history"
#: written by :func:`ICESEE.src.postprocessing.icepack_mesh.export_dof_coordinates`
DOF_COORDINATES = "dof_coordinates.h5"


class ResultsError(Exception):
    """A run's outputs cannot be read.

    ``kind`` is one of ``missing`` (a required file or dataset is absent),
    ``corrupt`` (a file exists but cannot be read) or ``incomplete`` (the
    files are readable but inconsistent, e.g. a run that stopped early).
    """

    def __init__(self, kind: str, message: str, *, files: tuple = ()):
        super().__init__(message)
        self.kind = kind
        self.message = message
        self.files = tuple(str(f) for f in files)


class Mode3History:
    """The ensemble an execution-mode-3 run keeps in ``_mode3_state_history``,
    indexed like the ``(state, member, time)`` ensemble of modes 0-2.

    ``initial_condition/checkpoint_00000000`` is time index 0 and
    ``steps/checkpoint_<k>`` (the ensemble after cycle ``k``) is time index
    ``k + 1`` (the convention of the mode-3 runners, e.g. Icepack's
    ``idealized_pig/mode3_runner.py``). Only complete checkpoints are used; time slots whose
    checkpoint was pruned (``checkpoint_keep_last``) read as NaN and are
    listed in :attr:`missing`. Each checkpoint is read when it is needed, so
    the history is never loaded whole.
    """

    def __init__(self, root: Path, n_times: int):
        self.root = Path(root)
        self.checkpoints: dict[int, tuple[Path, dict]] = {}
        candidates = [(0, self.root / "initial_condition" / "checkpoint_00000000")]
        for d in (self.root / "steps").glob("checkpoint_*"):
            try:
                candidates.append((int(d.name.split("_")[-1]) + 1, d))
            except ValueError:
                continue
        for index, d in candidates:
            try:
                manifest = json.loads((d / "manifest.json").read_text())
            except (OSError, ValueError):
                continue
            if manifest.get("complete") and 0 <= index < n_times:
                self.checkpoints[index] = (d, manifest)
        if not self.checkpoints:
            raise ResultsError("missing", f"{MODE3_HISTORY} holds no complete checkpoint", files=(self.root,))
        first = next(iter(self.checkpoints.values()))[1]
        self.shape = (int(first["global_size"]), int(first["number_of_members"]), int(n_times))
        self.missing = [i for i in range(n_times) if i not in self.checkpoints]

    def _members(self, index: int) -> np.ndarray:
        """``(state, member)`` at time ``index``; NaN when not kept."""
        out = np.full(self.shape[:2], np.nan)
        if index not in self.checkpoints:
            return out
        import h5py

        d, manifest = self.checkpoints[index]
        try:
            for shard in manifest["shards"]:
                with h5py.File(d / shard["file"], "r") as f:
                    for column, member in enumerate(shard["member_ids"]):
                        for block in shard["blocks"]:
                            start = int(block["global_offset"]) + int(block["owned_start"])
                            stop = int(block["global_offset"]) + int(block["owned_stop"])
                            out[start:stop, int(member)] = f[f"member_{member}/{block['name']}"][()]
        except Exception as exc:
            raise ResultsError("corrupt", f"{d.name} could not be read ({type(exc).__name__}: {exc})",
                               files=(d,)) from exc
        return out

    def __getitem__(self, key) -> np.ndarray:
        rows, members, times = key
        steps = range(self.shape[2])[times] if isinstance(times, slice) else [times]
        stacked = np.stack([self._members(i)[rows][..., members] for i in steps], axis=-1)
        return stacked if isinstance(times, slice) else stacked[..., 0]

    def mean(self) -> np.ndarray:
        """The ensemble mean, ``(state, time)``."""
        return np.stack([self._members(i).mean(axis=1) for i in range(self.shape[2])], axis=1)


def _read_observations(path: Path) -> np.ndarray:
    """The dense ``(state, obs time)`` observation matrix from either layout."""
    try:
        return _read(path, ("hu_obs",))["hu_obs"]
    except ResultsError as exc:
        if exc.kind != "missing" or not path.is_file():
            raise
    import h5py

    compact = _read(path, ("hu_obs_compact", "obs_indices"))
    indices = np.asarray(compact["obs_indices"], dtype=np.int64)
    try:
        with h5py.File(path, "r") as f:
            nd = int(f.attrs.get("state_dimension", (int(indices.max()) + 1) if indices.size else 0))
    except Exception as exc:
        raise ResultsError("corrupt", f"{path.name} could not be read ({type(exc).__name__}: {exc})",
                           files=(path,)) from exc
    dense = np.zeros((nd, compact["hu_obs_compact"].shape[1]), dtype=np.float64)
    dense[indices, :] = compact["hu_obs_compact"]
    return dense


def _read(path: Path, keys: tuple[str, ...], *, optional: tuple[str, ...] = ()) -> dict:
    import h5py

    if not path.is_file():
        raise ResultsError("missing", f"{path.name} was not found", files=(path,))
    try:
        with h5py.File(path, "r") as f:
            absent = [k for k in keys if k not in f]
            if absent:
                raise ResultsError(
                    "missing", f"{path.name} has no dataset {', '.join(absent)}", files=(path,))
            return {k: np.asarray(f[k][()]) for k in (*keys, *optional) if k in f}
    except ResultsError:
        raise
    except Exception as exc:  # truncated / not HDF5 / unreadable
        raise ResultsError("corrupt", f"{path.name} could not be read ({type(exc).__name__}: {exc})",
                           files=(path,)) from exc


@dataclass
class DAOutputs:
    """One run's data-assimilation outputs, as the notebooks load them."""

    data_dir: Path
    model_name: str
    filter_type: str
    t: np.ndarray
    obs_index: np.ndarray
    obs_max_time: float
    run_mode: int
    true_state: np.ndarray
    nurged_state: np.ndarray
    ensemble_mean: np.ndarray
    observations: np.ndarray
    metadata: dict = field(default_factory=dict)
    ensemble_file: Path | None = None
    ensemble_key: str = "ensemble"
    ensemble_bg: np.ndarray | None = None
    warnings: list = field(default_factory=list)

    @property
    def n_state(self) -> int:
        return int(self.true_state.shape[0])

    @property
    def n_times(self) -> int:
        return int(self.true_state.shape[1])

    @property
    def n_ens(self) -> int:
        shape = self._ensemble_shape()
        return int(shape[1]) if shape else 0

    def _ensemble_shape(self):
        if self.ensemble_file is None:
            return None
        with self.open_ensemble() as ds:
            return ds.shape

    @contextmanager
    def open_ensemble(self):
        """The ``(state, member, time)`` ensemble: an HDF5 dataset, or a
        :class:`Mode3History` for an execution-mode-3 run."""
        if self.ensemble_key == MODE3_HISTORY:
            yield Mode3History(self.ensemble_file, self.n_times)
            return
        import h5py

        with h5py.File(self.ensemble_file, "r") as f:
            yield f[self.ensemble_key]

    def ensemble_members(self, index) -> np.ndarray:
        """Members of the ensemble at state row(s) ``index``: an array shaped
        ``(time, member)`` for one row (``ensemble_vec_full[index, :, :].T`` in
        the notebooks). Read straight from the file, so a large ensemble is
        never loaded whole."""
        if self.ensemble_file is None:
            raise ResultsError("missing", "this run has no ensemble members on disk")
        try:
            with self.open_ensemble() as ds:
                return np.asarray(ds[index, :, :]).T
        except ResultsError:
            raise
        except Exception as exc:
            raise ResultsError("corrupt", f"ensemble members could not be read ({exc})",
                               files=(self.ensemble_file,)) from exc


def find_data_dir(run_dir: Path | str, *, model_name: str, data_dir_name: str = "_modelrun_datasets") -> Path | None:
    """The directory holding a run's ``true-wrong-<model>.h5``.

    ICESEE writes its datasets next to the example script, so depending on
    how the run was launched they sit in ``<run>/<data_dir_name>``,
    ``<run>/example/<example>/<data_dir_name>`` or, for outputs fetched from
    a cluster or the cloud, ``<run>/outputs/<data_dir_name>``. The newest
    candidate wins."""
    run_dir = Path(run_dir)
    name = TRUE_WRONG.format(model=model_name)
    candidates = [
        run_dir / data_dir_name,
        run_dir / "outputs" / data_dir_name,
        *sorted(run_dir.glob(f"example/*/{data_dir_name}")),
        *sorted(run_dir.glob(f"_remote_fetch/outputs/{data_dir_name}")),
        run_dir,
    ]
    found = [d for d in candidates if (d / name).is_file()]
    if not found:
        return None
    return max(found, key=lambda d: (d / name).stat().st_mtime)


def load_da_outputs(data_dir: Path | str, *, model_name: str, filter_type: str = "EnKF",
                    results_dir: Path | str | None = None) -> DAOutputs:
    """Load a run the way ``read_results.ipynb`` does (cells "true-wrong",
    "ensemble", "true/nurged states", "observations")."""
    data_dir = Path(data_dir)
    results_dir = Path(results_dir) if results_dir is not None else data_dir

    tw = _read(data_dir / TRUE_WRONG.format(model=model_name),
               ("t", "obs_index", "obs_max_time", "run_mode"),
               optional=("Lxy", "nxy", "b_io"))
    run_mode = int(np.ravel(tw["run_mode"])[0])
    states = _read(data_dir / STATES, ("true_state", "nurged_state"))
    observations = _read_observations(data_dir / OBSERVATIONS)
    ensemble_bg = None
    warnings = []
    history = data_dir / MODE3_HISTORY
    # mode 3 writes run_mode 0; Lorenz still writes the ensemble file, the
    # spatially distributed runners keep only checkpoints
    if run_mode == 0 and history.is_dir() and not (data_dir / ENSEMBLE).is_file():
        mode3 = Mode3History(history, int(states["true_state"].shape[1]))
        ensemble_file, ensemble_key = history, MODE3_HISTORY
        ensemble_mean = mode3.mean()
        if mode3.missing:
            warnings.append(f"{len(mode3.missing)} of {mode3.shape[2]} time steps have no kept checkpoint "
                            "(checkpoint_keep_last: 0 keeps them all) and are left blank")
    elif run_mode != 0:
        path = results_dir / FILTER_OUTPUT.format(filter=filter_type, model=model_name)
        ens = _read(path, ("ensemble_vec_mean",), optional=("ensemble_bg",))
        ensemble_file, ensemble_key = path, "ensemble_vec_full"
        ensemble_mean = ens["ensemble_vec_mean"]
        ensemble_bg = ens.get("ensemble_bg")
    else:
        path = data_dir / ENSEMBLE
        ens = _read(path, ("ensemble_mean",))
        ensemble_file, ensemble_key = path, "ensemble"
        ensemble_mean = ens["ensemble_mean"]

    if ensemble_key != MODE3_HISTORY:
        import h5py

        with h5py.File(ensemble_file, "r") as f:
            if ensemble_key not in f:
                ensemble_file = None

    out = DAOutputs(
        data_dir=data_dir, model_name=model_name, filter_type=filter_type,
        t=np.asarray(tw["t"], dtype=float), obs_index=np.asarray(tw["obs_index"]),
        obs_max_time=float(np.ravel(tw["obs_max_time"])[0]), run_mode=run_mode,
        true_state=states["true_state"], nurged_state=states["nurged_state"],
        ensemble_mean=ensemble_mean, observations=observations,
        metadata={k: tw[k] for k in ("Lxy", "nxy", "b_io") if k in tw},
        ensemble_file=ensemble_file, ensemble_key=ensemble_key, ensemble_bg=ensemble_bg,
        warnings=warnings,
    )
    _check_consistency(out)
    return out


def _check_consistency(out: DAOutputs) -> None:
    nd, nt = out.true_state.shape
    if out.nurged_state.shape != out.true_state.shape:
        raise ResultsError("incomplete", "true and background (nurged) states have different shapes "
                           f"{out.true_state.shape} vs {out.nurged_state.shape}")
    if out.ensemble_mean.shape[0] != nd:
        raise ResultsError("incomplete", f"the analysis has {out.ensemble_mean.shape[0]} state rows, "
                           f"the true state {nd}")
    if out.ensemble_mean.shape[1] < nt:
        out.warnings.append(f"the analysis covers {out.ensemble_mean.shape[1]} of {nt} time steps "
                            "(the run may have stopped early)")
    if out.t.size != nt:
        out.warnings.append(f"the time axis has {out.t.size} entries, the state {nt}")
    if out.observations.shape[0] != nd:
        raise ResultsError("incomplete", f"the observations have {out.observations.shape[0]} state rows, "
                           f"the true state {nd}")
    if out.ensemble_file is None:
        out.warnings.append("no ensemble members were saved; only the ensemble mean is shown")


def load_params(path: Path | str) -> dict:
    """The run's ``params.yaml`` (as ICESEE's ``load_yaml_to_dict`` reads it)."""
    import yaml

    path = Path(path)
    if not path.is_file():
        raise ResultsError("missing", f"{path.name} was not found", files=(path,))
    try:
        data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    except Exception as exc:
        raise ResultsError("corrupt", f"{path.name} is not valid YAML ({exc})", files=(path,)) from exc
    if not isinstance(data, dict):
        raise ResultsError("corrupt", f"{path.name} is not a mapping", files=(path,))
    return data


def load_dof_coordinates(data_dir: Path | str) -> np.ndarray:
    """``(n_dofs, 2)`` physical coordinates of the scalar-space degrees of
    freedom, in the state vector's ordering (Icepack runs)."""
    path = Path(data_dir) / DOF_COORDINATES
    if not path.is_file():
        raise ResultsError(
            "missing",
            "mesh coordinates were not exported for this run (dof_coordinates.h5); "
            "maps need the post-processing step that runs in the Icepack environment",
            files=(path,))
    return _read(path, ("coordinates",))["coordinates"]
