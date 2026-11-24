"""Plot real PIG output layouts without model/MPI imports or full-history reads."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest

from ICESEE.src.postprocessing.pig import PigResults, extract_datasets_from_h5, export_reference_coordinates

ROOT = Path(__file__).resolve().parents[2]
NOTEBOOK = ROOT / "applications/icepack_model/examples/idealized_pig/read_results.ipynb"


def write_pig(data, *, partial=False, background=False, compact=False):
    data.mkdir(parents=True)
    n, members, steps = 6, 3, 5
    t = np.arange(steps) * .25
    truth = np.arange(5*n*steps).reshape(5*n, steps).astype(float)
    ensemble = truth[:, None, :] + np.arange(members)[None, :, None]
    with h5py.File(data / "icesee_ensemble_data.h5", "w") as f:
        f["ensemble"] = ensemble
        f["ensemble_mean"] = ensemble.mean(axis=1)
        f["t"] = t
        f["planned_obs_index"] = [1, 3, 7] if partial else [1, 3]
        f.attrs.update(format="icesee-ensemble-history-v1", complete=not partial, export_complete=True,
                       requested_steps=8 if partial else steps-1,
                       state_order="serial_reference",
                       blocks=json.dumps([dict(name=name, global_offset=i*n, global_size=n)
                                          for i, name in enumerate(("h", "u", "v", "s", "basal_melt_field"))]))
    with h5py.File(data / "true_nurged_states.h5", "w") as f:
        f["true_state"] = truth
        if background:
            f["nurged_state"] = truth + 10
    with h5py.File(data / "synthetic_obs.h5", "w") as f:
        if compact:
            f["hu_obs_compact"] = truth[18:24, [1, 3]][:, :]
            f["obs_indices"] = np.arange(18, 24)
        else:
            f["hu_obs"] = truth[:, [1, 3]]
    (data.parent / "params.yaml").write_text("enkf-parameters:\n  observed_params: ['s']\n")
    for member in [None, *range(members)]:
        name = "hs_profiles_true.h5" if member is None else f"hs_ensemble_profiles{member}.h5"
        with h5py.File(data / name, "w") as f:
            f["h_profiles"] = np.ones((4, 4)) * (30 if member is None else 31+member)
            f["s_profiles"] = np.ones((4, 4)) * 40
            f["bed_values"] = np.arange(4)
            f["distances"] = np.arange(4)*1000
            f["profile_times"] = [0, .5, 1, 2]
            f["written_steps"] = [0, 2, 4, -1]
    return truth, ensemble


@pytest.mark.parametrize("background,compact,partial", [(False, False, False), (True, True, True)])
def test_point_reads_members_times_observed_rows_and_optional_background(tmp_path, background, compact, partial):
    truth, ensemble = write_pig(tmp_path / "data", partial=partial, background=background, compact=compact)
    with PigResults(tmp_path / "data") as r:
        assert isinstance(r.truth, h5py.Dataset)
        assert isinstance(r.ensemble, h5py.Dataset)
        s = r.point("s", 2, stride=2)
        np.testing.assert_array_equal(s["members"], ensemble[20, :, ::2].T)
        np.testing.assert_array_equal(s["truth"], truth[20, ::2])
        np.testing.assert_array_equal(s["observations"][0], [.25, .75])
        np.testing.assert_array_equal(s["observations"][1], truth[20, [1, 3]])
        h = r.point("h", 2)
        assert h["observations"] is None  # dense zero/unobserved rows must not become observations
        if background:
            np.testing.assert_array_equal(h["background"], truth[2] + 10)
        else:
            assert h["background"] is None and r.background is None
        speed = r.point("speed", 2)
        np.testing.assert_allclose(speed["members"], np.hypot(ensemble[8], ensemble[14]).T)
        assert len(r.t) == 5
        assert bool(r.metadata["warnings"]) == partial


def test_plotting_never_materializes_full_state_or_ensemble_history(tmp_path, monkeypatch):
    write_pig(tmp_path / "data")
    getitem = h5py.Dataset.__getitem__
    reads = []
    def bounded(dataset, key):
        if dataset.name in ("/ensemble", "/ensemble_mean", "/true_state", "/nurged_state"):
            reads.append((dataset.name, key))
            assert isinstance(key, tuple)
            assert any(isinstance(item, (int, np.integer)) for item in key), (dataset.name, key)
        return getitem(dataset, key)
    monkeypatch.setattr(h5py.Dataset, "__getitem__", bounded)
    with PigResults(tmp_path / "data") as r:
        point = r.select_point()
        r.plot_point(point)
        before = len(reads)
        r.plot_point_diagnostics(point)
        assert len(reads) == before  # cached local diagnostics do not reread disk
        fig = r.plot_profiles()
        assert "year 1" in fig._suptitle.get_text()
        r.plot_snapshot(.6)
        assert r.coordinates is None
    plt.close("all")


def test_memory_limits_fail_before_read(tmp_path):
    write_pig(tmp_path / "data")
    with PigResults(tmp_path / "data", max_read_bytes=16) as r:
        with pytest.raises(MemoryError, match="exceeds"):
            r.point("h", 0)
    huge = tmp_path / "huge.h5"
    with h5py.File(huge, "w") as f:
        f.create_dataset("ensemble", (100000, 40, 1641), chunks=(100, 1, 1))
    with pytest.raises(MemoryError, match="metadata"):
        extract_datasets_from_h5(huge)


def test_cached_references_and_missing_coordinates(tmp_path):
    write_pig(tmp_path / "reference")
    data = tmp_path / "case/data"
    data.mkdir(parents=True)
    (data / "icesee_ensemble_data.h5").symlink_to(tmp_path / "reference/icesee_ensemble_data.h5")
    (data.parent / "params.yaml").write_text(
        f"enkf-parameters:\n  observed_params: ['s']\n  reference_data_path: {tmp_path / 'reference'}\n")
    with PigResults(data) as r:
        assert r.truth_path.parent == tmp_path / "reference"
        with pytest.raises(ValueError, match="Physical point"):
            r.select_point(xy_m=(0, 0))
        with pytest.raises(ValueError, match="point_index"):
            r.point("h", 6)


def test_legacy_mode12_aggregate_with_summary_metadata(tmp_path):
    truth, ensemble = write_pig(tmp_path / "data", background=True)
    with h5py.File(tmp_path / "data/icesee_ensemble_data.h5", "a") as handle:
        del handle.attrs["format"]
        del handle.attrs["blocks"]
    with h5py.File(tmp_path / "data/true-wrong-icepack.h5", "w") as handle:
        handle["t"] = np.arange(5) * .25
        handle["obs_index"] = [1, 3]
        handle["obs_max_time"] = [.75]
        handle["run_mode"] = [0]
    with PigResults(tmp_path / "data") as results:
        np.testing.assert_array_equal(results.point("h", 1)["members"], ensemble[1].T)
        np.testing.assert_array_equal(results.point("h", 1)["background"], truth[1] + 10)


def test_coordinate_export_verifies_native_identity_and_reference_permutation(tmp_path):
    write_pig(tmp_path / "data")
    coords = np.arange(12).reshape(6, 2).astype(float)
    digest = hashlib.sha256(coords.tobytes()); digest.update(str(coords.shape).encode())
    fingerprint = "coords-sha256:" + digest.hexdigest()[:32]
    history = tmp_path / "data/_mode3_state_history"
    step = history / "steps/checkpoint_00000003"
    step.mkdir(parents=True)
    manifest = dict(complete=True, run_id="test", layout_id="mesh", shards=[dict(layout_fingerprint=fingerprint)])
    (step / "manifest.json").write_text(json.dumps(manifest))
    with h5py.File(history / "reference_layout.h5", "w") as f:
        f.attrs.update(run_id="test", layout_id="mesh")
        for name in ("h", "u", "v", "s", "basal_melt_field"):
            f[f"fields/{name}"] = np.arange(6)[::-1]
    initial = tmp_path / "initial.h5"
    with h5py.File(initial, "w") as f:
        f["topology/vecs/firedrake_default_coordinates/firedrake_default_coordinates"] = coords.ravel()
    export_reference_coordinates(tmp_path / "data", initial)
    with PigResults(tmp_path / "data") as r:
        np.testing.assert_array_equal(r.coordinates, coords[::-1])
        assert r.select_point(xy_m=(0, 1)) == 5
        r.plot_snapshot()
    with pytest.raises(FileExistsError):
        export_reference_coordinates(tmp_path / "data", initial)
    (tmp_path / "data/dof_coordinates.h5").unlink()
    manifest["shards"][0]["layout_fingerprint"] = "wrong"
    (step / "manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="fingerprint"):
        export_reference_coordinates(tmp_path / "data", initial)
    assert not (tmp_path / "data/dof_coordinates.h5").exists()
    plt.close("all")


def test_notebook_run_all_with_real_hdf5_datasets(tmp_path, monkeypatch):
    write_pig(tmp_path / "data")
    monkeypatch.setenv("ICESEE_DATA_PATH", str(tmp_path / "data"))
    namespace = {"__name__": "__main__"}
    try:
        for i, cell in enumerate(json.loads(NOTEBOOK.read_text())["cells"]):
            if cell["cell_type"] == "code":
                exec(compile("".join(cell["source"]), f"{NOTEBOOK}:cell{i}", "exec"), namespace)
        assert len(namespace["figures"]) == 5
        assert namespace["ensemble_nurged_state"] is None
        for name, figure in namespace["figures"].items():
            figure.savefig(tmp_path / f"{name}.png")
    finally:
        if "results" in namespace:
            namespace["results"].close()
        if "_ensemble_file" in namespace:
            namespace["_ensemble_file"].close()
        plt.close("all")


def test_reader_import_does_not_initialize_mpi_or_model():
    code = ("import sys; from ICESEE.src.postprocessing.pig import PigResults, extract_datasets_from_h5; "
            "assert 'mpi4py.MPI' not in sys.modules; assert 'firedrake' not in sys.modules; "
            "assert 'ICESEE.config._utility_imports' not in sys.modules")
    subprocess.run([sys.executable, "-c", code], check=True, env=os.environ.copy())
