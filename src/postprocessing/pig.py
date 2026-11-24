"""Bounded, MPI-free readers and plots for Idealized PIG output."""
from __future__ import annotations

import hashlib
import json
import os
import tempfile
from collections import OrderedDict
from pathlib import Path

import h5py
import numpy as np

from .outputs import native_history_metadata

FIELDS = ("h", "u", "v", "s", "basal_melt_field")
LABELS = {"h": "Thickness (m)", "u": "x velocity (m/year)",
          "v": "y velocity (m/year)", "s": "Surface elevation (m)",
          "speed": "Speed (m/year)", "basal_melt_field": "Basal melt (m/year)"}


def extract_datasets_from_h5(path, *, max_bytes=32 * 1024**2):
    """Read small metadata files; refuse large arrays before allocating them.

    Use :class:`PigResults` or h5py slices for state and ensemble histories.
    Dataset names include their HDF5 group path.
    """
    with h5py.File(path, "r") as handle:
        datasets = {}
        handle.visititems(lambda name, item: datasets.update({name: item})
                          if isinstance(item, h5py.Dataset) else None)
        if sum(d.size * d.dtype.itemsize for d in datasets.values()) > max_bytes:
            raise MemoryError("File exceeds the metadata read limit; select HDF5 rows/times instead.")
        return {name: dataset[()] for name, dataset in datasets.items()}


class PigResults:
    """Keep histories on disk; materialize only a point or a spatial snapshot.

    ``reference_dir`` can identify shared truth/observations for a cached-input
    run. Otherwise it is read from the run's saved params.yaml when available.
    A missing background is represented by None, never a full NaN trajectory.
    """

    def __init__(self, data_dir, *, reference_dir=None, max_read_bytes=128 * 1024**2):
        self.data_dir = Path(data_dir).expanduser().resolve()
        self.max_read_bytes = int(max_read_bytes)
        self._files = []
        self._cache = OrderedDict()
        try:
            self.params = {}
            for candidate in (self.data_dir.parent / "params.yaml", self.data_dir / "params.yaml"):
                if candidate.is_file():
                    import yaml
                    self.params = yaml.safe_load(candidate.read_text()) or {}
                    break
            flat = {key: value for section in self.params.values() if isinstance(section, dict)
                    for key, value in section.items()}
            if reference_dir is None:
                reference_dir = flat.get("reference_data_path")
            self.reference_dir = self.data_dir
            if reference_dir:
                reference = Path(reference_dir).expanduser()
                if not reference.is_absolute():
                    reference = self.data_dir.parent / reference
                self.reference_dir = reference.resolve()
            self.observed_variables = set(flat.get("observed_vars", [])) | set(flat.get("observed_params", []))
            summary = self.data_dir / "true-wrong-icepack.h5"
            self.summary = extract_datasets_from_h5(summary) if summary.is_file() else {}
            aggregate = self._open(self.data_dir / "icesee_ensemble_data.h5")
            self.ensemble = aggregate["ensemble"]
            self.mean = aggregate["ensemble_mean"]
            if aggregate.attrs.get("format") == "icesee-ensemble-history-v1":
                self.metadata = native_history_metadata(aggregate.filename, summary=self.summary)
            else:
                times = np.asarray(aggregate["t"][:] if "t" in aggregate else self.summary["t"])
                planned = np.asarray(self.summary.get("obs_index", []), dtype=int)
                available = (planned >= 0) & (planned < times.size)
                self.metadata = dict(t=times, obs_index=planned[available],
                                     obs_column_indices=np.flatnonzero(available), warnings=[])
            self.t = self.metadata["t"]
            if self.ensemble.ndim != 3 or self.ensemble.shape[2] != len(self.t):
                raise ValueError("Ensemble and saved time dimensions disagree")
            if self.mean.shape != (self.ensemble.shape[0], len(self.t)):
                raise ValueError("Ensemble mean has an incompatible shape")
            definitions = json.loads(aggregate.attrs["blocks"]) if "blocks" in aggregate.attrs else None
            if definitions is None:
                if self.ensemble.shape[0] % len(FIELDS):
                    raise ValueError("Missing PIG block layout; cannot infer field rows")
                n = self.ensemble.shape[0] // len(FIELDS)
                definitions = [dict(name=name, global_offset=i*n, global_size=n)
                               for i, name in enumerate(FIELDS)]
            self.blocks = {b["name"]: (int(b["global_offset"]), int(b["global_size"])) for b in definitions}
            if not set(FIELDS).issubset(self.blocks):
                raise ValueError("This reader requires PIG's h/u/v/s/basal_melt_field layout")
            self.n_points = self.blocks["h"][1]
            if any(self.blocks[name][1] != self.n_points for name in FIELDS):
                raise ValueError("Point plots require matching scalar-space field sizes")
            if sorted((start, start+size) for start, size in self.blocks.values()) != [
                    (i*self.n_points, (i+1)*self.n_points) for i in range(len(FIELDS))]:
                raise ValueError("PIG field blocks do not cover the packed state")
            self.truth_path = self._reference_file("true_nurged_states.h5")
            states = self._open(self.truth_path)
            self.truth = states["true_state"]
            self.background = states.get("nurged_state")
            if self.truth.shape[0] != self.ensemble.shape[0] or self.truth.shape[1] < len(self.t):
                raise ValueError("Truth shape does not cover this ensemble history")
            if self.background is not None and (self.background.shape[0] != self.truth.shape[0]
                                                or self.background.shape[1] < len(self.t)):
                raise ValueError("Background shape does not cover this ensemble history")
            obs_path = self._reference_file("synthetic_obs.h5", required=False)
            self.observations = self._open(obs_path) if obs_path else None
            self.coordinates = self._coordinates()
        except Exception:
            self.close()
            raise

    def _open(self, path):
        handle = h5py.File(path, "r", rdcc_nbytes=0)
        self._files.append(handle)
        return handle

    def _reference_file(self, name, *, required=True):
        for directory in dict.fromkeys((self.data_dir, self.reference_dir)):
            path = directory / name
            if path.is_file():
                return path
        if required:
            raise FileNotFoundError(f"{name} is missing; set reference_dir to the cached reference data.")
        return None

    def _read(self, dataset, selection, shape):
        if int(np.prod(shape, dtype=np.int64)) * dataset.dtype.itemsize > self.max_read_bytes:
            raise MemoryError("Selected plot exceeds max_read_bytes; increase time_stride or select fewer rows.")
        return np.asarray(dataset[selection])

    def _coordinates(self):
        path = self._reference_file("dof_coordinates.h5", required=False)
        if path is None:
            return None
        with h5py.File(path, "r") as handle:
            if handle.attrs.get("coordinate_order") != "serial_reference":
                raise ValueError("Plot coordinates must declare serial_reference row order")
            reference = handle.attrs.get("reference_truth_file")
            if reference and Path(reference).resolve() != self.truth_path.resolve():
                raise ValueError("Plot coordinates refer to a different truth file")
            if handle["coordinates"].shape != (self.n_points, 2):
                raise ValueError("Plot coordinates do not match the PIG state")
            coords = self._read(handle["coordinates"], np.s_[:, :], (self.n_points, 2))
        if coords.shape != (self.n_points, 2) or not np.isfinite(coords).all():
            raise ValueError("Plot coordinates do not match the PIG state")
        return coords

    def select_point(self, index=None, xy_m=None):
        if xy_m is not None:
            if self.coordinates is None:
                raise ValueError("Physical point selection needs dof_coordinates.h5")
            return int(np.argmin(np.sum((self.coordinates - np.asarray(xy_m))**2, axis=1)))
        if index is None:
            offset, size = self.blocks["basal_melt_field"]
            melt = self._read(self.mean, np.s_[offset:offset+size, -1], (size,))
            active = np.flatnonzero(np.isfinite(melt) & (np.abs(melt) > 0))
            index = active[len(active)//2] if active.size else self.n_points//2
        if int(index) != index or not 0 <= int(index) < self.n_points:
            raise ValueError(f"point_index must be in 0..{self.n_points-1}")
        return int(index)

    def point(self, variable, index, *, stride=1):
        index = self.select_point(index)
        if int(stride) != stride or stride < 1:
            raise ValueError("time_stride must be a positive integer")
        stride = int(stride)
        key = (variable, index, stride)
        if key in self._cache:
            self._cache.move_to_end(key)
            return self._cache[key]
        times = self.t[::stride]
        count = len(times)
        if variable == "speed":
            u, v = (self.point(name, index, stride=stride) for name in ("u", "v"))
            members = np.hypot(u["members"], v["members"])
            result = dict(t=times, members=members, mean=members.mean(axis=1),
                          truth=np.hypot(u["truth"], v["truth"]), background=None, observations=None)
            if u["background"] is not None:
                result["background"] = np.hypot(u["background"], v["background"])
        else:
            row = self.blocks[variable][0] + index
            members = self._read(self.ensemble, np.s_[row, :, ::stride],
                                 (self.ensemble.shape[1], count)).T
            result = dict(t=times, members=members,
                          mean=self._read(self.mean, np.s_[row, ::stride], (count,)),
                          truth=self._read(self.truth, np.s_[row, :len(self.t):stride], (count,)),
                          background=None, observations=self._point_observations(variable, row))
            if self.background is not None:
                result["background"] = self._read(self.background, np.s_[row, :len(self.t):stride], (count,))
        size = sum(x.nbytes for x in result.values() if isinstance(x, np.ndarray))
        while self._cache and self._cache_bytes() + size > self.max_read_bytes:
            self._cache.popitem(last=False)
        if size <= self.max_read_bytes:
            self._cache[key] = result
        return result

    def _cache_bytes(self):
        return sum(x.nbytes for result in self._cache.values() for x in result.values()
                   if isinstance(x, np.ndarray))

    def _point_observations(self, variable, row):
        if self.observations is None or variable not in self.observed_variables:
            return None
        handle = self.observations
        columns = np.asarray(self.metadata["obs_column_indices"], dtype=int)
        if "hu_obs" in handle:
            dataset = handle["hu_obs"]
        elif "hu_obs_compact" in handle and "obs_indices" in handle:
            indices = self._read(handle["obs_indices"], np.s_[:], handle["obs_indices"].shape)
            matches = np.flatnonzero(indices == row)
            if not matches.size:
                return None
            row, dataset = int(matches[0]), handle["hu_obs_compact"]
        else:
            return None
        if not columns.size:
            return None
        if row >= dataset.shape[0] or columns[-1] >= dataset.shape[1]:
            raise ValueError("Observation rows/columns disagree with saved schedule")
        values = self._read(dataset, (row, columns), (len(columns),))
        return self.t[self.metadata["obs_index"]], values

    def plot_point(self, index, *, stride=1, show_members=True):
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(3, 2, figsize=(13, 11), sharex=True, constrained_layout=True)
        for ax, variable in zip(axes.flat, ("h", "s", "u", "v", "speed", "basal_melt_field")):
            series = self.point(variable, index, stride=stride)
            time, members = series["t"], series["members"]
            if show_members:
                ax.plot(time, members, color="0.65", alpha=.3, lw=.6)
            low, high = np.quantile(members, [.1, .9], axis=1)
            ax.fill_between(time, low, high, alpha=.2, color="tab:blue", label="10–90% of members")
            ax.plot(time, series["mean"], color="tab:blue", label="Ensemble mean")
            ax.plot(time, series["truth"], color="tab:red", label="Truth")
            if series["background"] is not None:
                ax.plot(time, series["background"], color="tab:green", ls="--", label="Background")
            if series["observations"] is not None:
                ax.scatter(*series["observations"], color="black", s=14, marker="x", label="Observations")
            ax.set(ylabel=LABELS[variable], xlabel="Model time (years)")
            ax.grid(alpha=.2)
        axes[0, 0].legend(fontsize=8)
        axes[0, 1].legend(fontsize=8)
        fig.suptitle(f"PIG point {index}: member spread (speed mean is the mean of member speeds)")
        return fig

    def plot_point_diagnostics(self, index, *, stride=1):
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(3, 2, figsize=(13, 9), sharex=True, constrained_layout=True)
        for ax, name in zip(axes.flat, ("h", "s", "u", "v", "speed", "basal_melt_field")):
            values = self.point(name, index, stride=stride)
            ax.plot(values["t"], np.abs(values["mean"] - values["truth"]), label="Absolute mean error")
            ax.plot(values["t"], values["members"].std(axis=1, ddof=1), label="Ensemble SD")
            ax.set(ylabel=LABELS[name], xlabel="Model time (years)")
            ax.grid(alpha=.2)
        axes[0, 0].legend()
        fig.suptitle(f"Point {index} diagnostics; these are not domain-wide RMSE")
        return fig

    def _profile(self, path):
        with h5py.File(path, "r") as handle:
            hkey = "h_profiles" if "h_profiles" in handle else "h_nudged_profiles"
            skey = "s_profiles" if "s_profiles" in handle else "s_nudged_profiles"
            h = self._read(handle[hkey], np.s_[:, :], handle[hkey].shape)
            s = self._read(handle[skey], np.s_[:, :], handle[skey].shape)
            if "profile_times" not in handle:
                raise ValueError(f"{path.name} has no profile_times; supply verified profile metadata first")
            times = handle["profile_times"][:]
            steps = handle["written_steps"][:] if "written_steps" in handle else np.full(len(times), -1)
            valid = (steps >= 0) & np.isfinite(times) & (times <= self.t[-1] + 1e-8)
            valid &= np.isfinite(h).all(axis=0) & np.isfinite(s).all(axis=0)
            distances = handle["distances"][:] if "distances" in handle else None
            bed = handle["bed_values"][:] if "bed_values" in handle else None
        return dict(h=h, s=s, t=times, valid=valid, distances=distances, bed=bed)

    def plot_profiles(self, year=None):
        import matplotlib.pyplot as plt
        truth_path = self._reference_file("hs_profiles_true.h5", required=False)
        if truth_path is None:
            truth_path = self._reference_file("hs_profiles_true", required=False)
        paths = [self.data_dir / f"hs_ensemble_profiles{i}.h5" for i in range(self.ensemble.shape[1])]
        if truth_path is None or not all(path.is_file() for path in paths):
            raise FileNotFoundError("Truth or member flowline profiles are missing")
        truth = self._profile(truth_path)
        profiles = [self._profile(path) for path in paths]
        if truth["distances"] is None:
            raise ValueError("Truth profiles have no physical flowline distances")
        for profile in profiles:
            if (profile["distances"] is None or profile["h"].shape[0] != len(truth["distances"])
                    or not np.allclose(profile["distances"], truth["distances"])):
                raise ValueError("Member and truth flowline sampling points disagree")
        common = truth["t"][truth["valid"]]
        for profile in profiles:
            common = np.array([t for t in common if np.any(profile["valid"] & np.isclose(profile["t"], t))])
        if not common.size:
            raise ValueError("No common written flowline time is available")
        chosen = common[-1] if year is None else common[np.argmin(np.abs(common - float(year)))]
        def column(profile, field):
            index = np.flatnonzero(profile["valid"] & np.isclose(profile["t"], chosen))[0]
            return profile[field][:, index]
        x = truth["distances"] / 1000
        fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True, constrained_layout=True)
        for profile in profiles:
            axes[0].plot(x, column(profile, "s"), color="0.6", alpha=.25, lw=.6)
            axes[1].plot(x, column(profile, "h"), color="0.6", alpha=.25, lw=.6)
        surface = np.mean([column(p, "s") for p in profiles], axis=0)
        thickness = np.mean([column(p, "h") for p in profiles], axis=0)
        axes[0].plot(x, column(truth, "s"), color="tab:red", label="True surface")
        axes[0].plot(x, column(truth, "s") - column(truth, "h"), color="tab:red", ls=":", label="True base")
        axes[0].plot(x, surface, color="tab:blue", label="Mean surface")
        axes[0].plot(x, surface-thickness, color="tab:blue", ls=":", label="Mean base")
        if truth["bed"] is not None:
            axes[0].plot(x, truth["bed"], color="saddlebrown", label="Bed")
        axes[1].plot(x, column(truth, "h"), color="tab:red", label="Truth")
        axes[1].plot(x, thickness, color="tab:blue", label="Ensemble mean")
        background_path = self._reference_file("hs_profiles_wrong.h5", required=False)
        if background_path is not None:
            background = self._profile(background_path)
            if np.any(background["valid"] & np.isclose(background["t"], chosen)):
                if (background["distances"] is None
                        or not np.allclose(background["distances"], truth["distances"])):
                    raise ValueError("Background and truth flowline sampling points disagree")
                axes[0].plot(x, column(background, "s"), color="tab:green", ls="--", label="Background surface")
                axes[1].plot(x, column(background, "h"), color="tab:green", ls="--", label="Background")
        for ax in axes:
            ax.legend(fontsize=8)
            ax.grid(alpha=.2)
        axes[0].set_ylabel("Elevation (m)")
        axes[1].set(xlabel="Flowline distance (km)", ylabel="Thickness (m)")
        fig.suptitle(f"Saved flowline profiles at year {chosen:g}; available: {common.tolist()}")
        return fig

    def snapshot(self, variable, year=None):
        column = len(self.t)-1 if year is None else int(np.argmin(np.abs(self.t-float(year))))
        if variable == "speed":
            _, u = self.snapshot("u", self.t[column])
            _, v = self.snapshot("v", self.t[column])
            return self.t[column], {key: None if u[key] is None else np.hypot(u[key], v[key]) for key in u}
        offset, count = self.blocks[variable]
        values = {}
        for name, dataset in (("Truth", self.truth), ("Background", self.background), ("Analysis mean", self.mean)):
            values[name] = None if dataset is None else self._read(dataset, np.s_[offset:offset+count, column], (count,))
        return self.t[column], values

    def plot_snapshot(self, year=None, *, variables=("h", "speed")):
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(len(variables), 3, figsize=(15, 4*len(variables)), squeeze=False, constrained_layout=True)
        for row, name in enumerate(variables):
            chosen, values = self.snapshot(name, year)
            finite = np.concatenate([a[np.isfinite(a)] for a in values.values() if a is not None])
            low, high = (finite.min(), finite.max()) if finite.size else (0, 1)
            for ax, (label, array) in zip(axes[row], values.items()):
                ax.set_title(f"{label}: {LABELS[name]}")
                if array is None:
                    ax.text(.5, .5, "Background not generated", ha="center", transform=ax.transAxes)
                    ax.set_axis_off()
                elif self.coordinates is None:
                    ax.plot(np.arange(len(array)), array, lw=.6)
                    ax.set(xlabel="Scalar-space DOF index (coordinates unavailable)", ylabel=LABELS[name], ylim=(low, high))
                else:
                    artist = ax.scatter(self.coordinates[:, 0]/1000, self.coordinates[:, 1]/1000,
                                        c=array, s=3, vmin=low, vmax=high, rasterized=True)
                    ax.set(xlabel="x (km)", ylabel="y (km)", aspect="equal")
                    fig.colorbar(artist, ax=ax, shrink=.8)
        note = "speed is magnitude of the mean velocity vector" if "speed" in variables else ""
        fig.suptitle(f"PIG snapshot at year {chosen:g}; {note}")
        return fig

    def close(self):
        for handle in getattr(self, "_files", []):
            handle.close()
        self._files = []
        self._cache.clear()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()


def export_reference_coordinates(data_dir, initial_file):
    """Export PIG plotting coordinates only after verifying native row identity.

    This reads raw coordinate vectors from the initial Firedrake checkpoint;
    it does not import Firedrake or rebuild/repartition the mesh.
    """
    data_dir = Path(data_dir).resolve()
    target = data_dir / "dof_coordinates.h5"
    if target.exists():
        raise FileExistsError(target)
    history = data_dir / "_mode3_state_history"
    manifests = sorted((history / "steps").glob("checkpoint_*/manifest.json"), reverse=True)
    manifests += list((history / "initial_condition").glob("checkpoint_*/manifest.json"))
    manifest = next((m for p in manifests if (m := json.loads(p.read_text())).get("complete")), None)
    if manifest is None:
        raise ValueError("A complete native checkpoint is required to verify coordinate order")
    with h5py.File(initial_file, "r") as handle:
        vectors = []
        handle.visititems(lambda name, item: vectors.append(name) if isinstance(item, h5py.Dataset)
                          and name.endswith("/firedrake_default_coordinates/firedrake_default_coordinates") else None)
        if len(vectors) != 1:
            raise ValueError("Initial checkpoint must contain one identifiable coordinate vector")
        coords = np.asarray(handle[vectors[0]][:]).reshape(-1, 2)
    canonical = np.ascontiguousarray(np.round(coords, 3), dtype=np.float64)
    digest = hashlib.sha256(canonical.tobytes())
    digest.update(str(canonical.shape).encode())
    fingerprint = "coords-sha256:" + digest.hexdigest()[:32]
    if any(shard.get("layout_fingerprint") != fingerprint for shard in manifest["shards"]):
        raise ValueError("Initial coordinates do not match the saved native layout fingerprint")
    with h5py.File(history / "reference_layout.h5", "r") as handle:
        if any(handle.attrs.get(key) != manifest[key] for key in ("run_id", "layout_id")):
            raise ValueError("Reference map does not match this run")
        mapping = handle["fields/h"][:]
        if (mapping.dtype.kind not in "iu" or mapping.shape != (len(coords),)
                or not np.array_equal(np.sort(mapping), np.arange(len(coords)))):
            raise ValueError("Invalid reference/native coordinate permutation")
        if any(not np.array_equal(handle[f"fields/{name}"][:], mapping) for name in FIELDS):
            raise ValueError("PIG fields use different scalar-space orderings")
    fd, temporary = tempfile.mkstemp(prefix=".coordinates-", suffix=".h5", dir=data_dir)
    os.close(fd)
    try:
        with h5py.File(temporary, "w") as handle:
            handle["coordinates"] = coords[mapping]
            handle.attrs.update(coordinate_order="serial_reference", source_fingerprint=fingerprint,
                                reference_truth_file=str(data_dir / "true_nurged_states.h5"),
                                run_id=manifest["run_id"], layout_id=manifest["layout_id"])
        os.link(temporary, target)  # fail rather than replace an existing export
    finally:
        os.unlink(temporary)
    return target
