"""ICESEE post-processing: loading run outputs, the results contract, and the
result views, on small synthetic runs with the real file layout."""
from __future__ import annotations

import json
from pathlib import Path

import h5py
import numpy as np
import pytest

from ICESEE.src.postprocessing import (
    ContractError, ResultsError, find_data_dir, load_contract, load_da_outputs, parse_contract,
)
from ICESEE.src.postprocessing.render import figure_data, main as render_main, open_context, render_view
from ICESEE.src.postprocessing.views import all_views, get_view

EXAMPLES = Path(__file__).resolve().parents[2] / "applications"


def write_run(d: Path, *, model="lorenz", nd=3, nt=101, nens=4, nobs=5, run_mode=0, seed=0,
              ens_steps=None, extra_tw=None):
    rng = np.random.default_rng(seed)
    d.mkdir(parents=True, exist_ok=True)
    t = np.linspace(0, 1, nt)
    obs_index = np.linspace(nt // 20, nt // 2, nobs).astype(int)
    truth = rng.normal(size=(nd, nt))
    with h5py.File(d / f"true-wrong-{model}.h5", "w") as f:
        f["t"] = t
        f["obs_index"] = obs_index
        f["obs_max_time"] = np.array([0.5])
        f["run_mode"] = np.array([run_mode])
        f["Lxy"] = np.array([1.0, 1.0])
        f["nxy"] = np.array([1, 1])
        f["b_io"] = np.array([0.0, 0.0])
        for k, v in (extra_tw or {}).items():
            del f[k]
            f[k] = v
    with h5py.File(d / "true_nurged_states.h5", "w") as f:
        f["true_state"] = truth
        f["nurged_state"] = truth + 1.0
    with h5py.File(d / "synthetic_obs.h5", "w") as f:
        f["hu_obs"] = truth[:, obs_index] + 0.1
    ns = ens_steps or nt
    members = truth[:, None, :ns] + rng.normal(scale=0.1, size=(nd, nens, ns))
    with h5py.File(d / "icesee_ensemble_data.h5", "w") as f:
        f["ensemble"] = members
        f["ensemble_mean"] = members.mean(axis=1)
    return truth, members


CONTRACT = {"schema": "icesee-results/1", "model_name": "lorenz",
            "views": ["lorenz.state_estimation", "da.error_and_spread"]}


# ── loading ───────────────────────────────────────────────────────────────
def test_load_reads_the_notebook_files(tmp_path):
    truth, members = write_run(tmp_path / "_modelrun_datasets")
    out = load_da_outputs(tmp_path / "_modelrun_datasets", model_name="lorenz")
    assert out.run_mode == 0 and out.n_ens == 4 and out.n_state == 3
    np.testing.assert_array_equal(out.true_state, truth)
    np.testing.assert_allclose(out.ensemble_mean, members.mean(axis=1))
    np.testing.assert_allclose(out.ensemble_members(1), members[1].T)
    assert out.obs_max_time == 0.5 and not out.warnings


def test_missing_file_is_reported_by_name(tmp_path):
    d = tmp_path / "data"
    write_run(d)
    (d / "synthetic_obs.h5").unlink()
    with pytest.raises(ResultsError) as err:
        load_da_outputs(d, model_name="lorenz")
    assert err.value.kind == "missing" and "synthetic_obs.h5" in err.value.message


def test_corrupt_file_is_not_a_crash(tmp_path):
    d = tmp_path / "data"
    write_run(d)
    (d / "true_nurged_states.h5").write_bytes(b"not hdf5 at all")
    with pytest.raises(ResultsError) as err:
        load_da_outputs(d, model_name="lorenz")
    assert err.value.kind == "corrupt"


def test_missing_dataset_key(tmp_path):
    d = tmp_path / "data"
    write_run(d)
    with h5py.File(d / "true_nurged_states.h5", "a") as f:
        del f["nurged_state"]
    with pytest.raises(ResultsError) as err:
        load_da_outputs(d, model_name="lorenz")
    assert err.value.kind == "missing" and "nurged_state" in err.value.message


def test_incomplete_run_is_drawn_with_a_warning(tmp_path):
    d = tmp_path / "data"
    write_run(d, ens_steps=40)
    out = load_da_outputs(d, model_name="lorenz")
    assert any("40 of 101" in w for w in out.warnings)
    ctx = open_context(d, parse_contract(CONTRACT))
    rendered, fig = render_view(ctx, "lorenz.state_estimation", png=tmp_path / "x.png")
    analysis = [ln for ax in figure_data(fig)[:3] for ln in ax["lines"] if ln["label"] == "Analysis"]
    assert all(len(ln["x"]) == 40 for ln in analysis)
    assert any("40 of 101" in n for n in rendered.notes)


def test_inconsistent_shapes_are_incomplete(tmp_path):
    d = tmp_path / "data"
    write_run(d)
    with h5py.File(d / "synthetic_obs.h5", "a") as f:
        del f["hu_obs"]
        f["hu_obs"] = np.zeros((2, 5))
    with pytest.raises(ResultsError) as err:
        load_da_outputs(d, model_name="lorenz")
    assert err.value.kind == "incomplete"


def test_compact_observations_are_scattered_to_dense(tmp_path):
    d = tmp_path / "data"
    write_run(d)
    with h5py.File(d / "synthetic_obs.h5", "a") as f:
        dense = f["hu_obs"][:]
        del f["hu_obs"]
        f["hu_obs_compact"] = dense[[0, 2]]
        f["obs_indices"] = np.array([0, 2])
        f.attrs["state_dimension"] = 3
    out = load_da_outputs(d, model_name="lorenz")
    expected = dense.copy()
    expected[1] = 0.0
    np.testing.assert_array_equal(out.observations, expected)


def write_mode3_history(d: Path, members: np.ndarray, blocks: dict, *, skip=()):
    """Execution-mode-3 checkpoints of ``members`` (state, member, time) in
    the layout of ICESEE.src.parallelization.mode3_checkpointing, two shards
    each holding half of the members and splitting the first block."""
    import json

    nd, nens, nt = members.shape
    (d / ENSEMBLE_FILE).unlink()
    for index in range(nt):
        if index in skip:
            continue
        ckpt = (d / "_mode3_state_history" / "initial_condition" / "checkpoint_00000000" if index == 0
                else d / "_mode3_state_history" / "steps" / f"checkpoint_{index - 1:08d}")
        ckpt.mkdir(parents=True)
        shards = []
        for s, ids in enumerate((list(range(nens // 2)), list(range(nens // 2, nens)))):
            name = f"shard_{s:08d}.h5"
            shard_blocks = [{"name": n, "global_offset": off, "global_size": size,
                             "owned_start": 0, "owned_stop": size} for n, (off, size) in blocks.items()]
            with h5py.File(ckpt / name, "w") as f:
                for m in ids:
                    for b in shard_blocks:
                        f[f"member_{m}/{b['name']}"] = members[b["global_offset"]:b["global_offset"]
                                                               + b["global_size"], m, index]
            shards.append({"file": name, "member_ids": ids, "blocks": shard_blocks})
        (ckpt / "manifest.json").write_text(json.dumps({
            "complete": True, "global_size": nd, "number_of_members": nens, "shards": shards,
            "timestep": index - 1}))


ENSEMBLE_FILE = "icesee_ensemble_data.h5"


def test_mode3_history_reads_like_an_ensemble_file(tmp_path):
    d = tmp_path / "data"
    truth, members = write_run(d, nd=5, nt=12, nens=4)
    write_mode3_history(d, members, {"h": (0, 2), "u": (2, 2), "xg": (4, 1)})
    out = load_da_outputs(d, model_name="lorenz")
    assert out.ensemble_key == "_mode3_state_history" and out.n_ens == 4 and not out.warnings
    np.testing.assert_allclose(out.ensemble_mean, members.mean(axis=1))
    np.testing.assert_allclose(out.ensemble_members(3), members[3].T)
    ctx = open_context(d, parse_contract(CONTRACT))
    rendered, _ = render_view(ctx, "da.error_and_spread", png=tmp_path / "x.png")
    assert (tmp_path / "x.png").is_file()


def test_mode3_pruned_checkpoints_are_blank_and_reported(tmp_path):
    d = tmp_path / "data"
    truth, members = write_run(d, nd=3, nt=10, nens=2)
    write_mode3_history(d, members, {"x": (0, 3)}, skip=(4, 5))
    out = load_da_outputs(d, model_name="lorenz")
    assert np.isnan(out.ensemble_mean[:, 4]).all() and not np.isnan(out.ensemble_mean[:, 6]).any()
    assert any("2 of 10 time steps" in w for w in out.warnings)


def test_flowline_view_scales_to_physical_units(tmp_path):
    d = tmp_path / "data"
    nx = 4
    truth, members = write_run(d, model="flowline", nd=2 * nx + 1, nt=21, nens=4)
    params = {"physical-parameters": {"hscale": 1000.0, "rho_ice": 900.0, "g": 9.81, "C": 3e6,
                                      "N1": 3, "N2": 1},
              "modeling-parameters": {"year": 31536000}}
    contract = parse_contract({"schema": "icesee-results/1", "model_name": "flowline",
                               "views": ["flowline.state_estimation"]})
    ctx = open_context(d, contract, params=params)
    _, fig = render_view(ctx, "flowline.state_estimation", png=tmp_path / "f.png")
    data = figure_data(fig)
    truth_gl = next(ln for ln in data[0]["lines"] if ln["label"] == "Truth")
    xscale = 900.0 * 9.81 * 1000.0 / 3e6 * 1000.0 / 1e3
    np.testing.assert_allclose(truth_gl["y"], truth[2 * nx] * xscale)
    assert [ax["ylabel"] for ax in data] == ["Grounding line (km)", "h at node 2 (m)", "u at terminus (m/yr)"]


def test_filter_outputs_when_run_mode_is_parallel(tmp_path):
    d = tmp_path / "data"
    truth, members = write_run(d, run_mode=1)
    with h5py.File(d / "DEnKF-lorenz.h5", "w") as f:
        f["ensemble_vec_full"] = members * 2
        f["ensemble_vec_mean"] = members.mean(axis=1) * 2
        f["ensemble_bg"] = truth
    out = load_da_outputs(d, model_name="lorenz", filter_type="DEnKF")
    np.testing.assert_allclose(out.ensemble_mean, members.mean(axis=1) * 2)
    assert out.ensemble_key == "ensemble_vec_full"
    with pytest.raises(ResultsError):
        load_da_outputs(d, model_name="lorenz", filter_type="EnKF")


def test_find_data_dir_prefers_the_newest_candidate(tmp_path):
    import os
    import time

    write_run(tmp_path / "example" / "lorenz96" / "_modelrun_datasets")
    write_run(tmp_path / "_modelrun_datasets")
    newer = tmp_path / "example" / "lorenz96" / "_modelrun_datasets" / "true-wrong-lorenz.h5"
    os.utime(newer, (time.time() + 10, time.time() + 10))
    assert find_data_dir(tmp_path, model_name="lorenz") == newer.parent
    assert find_data_dir(tmp_path, model_name="icepack") is None


# ── contract ──────────────────────────────────────────────────────────────
@pytest.mark.parametrize("example", [
    "lorenz_model/examples/lorenz96", "icepack_model/examples/synthetic_ice_stream",
    "issm_model/examples/ISMIP_Choi", "flowline_model/examples/flowline_1d"])
def test_shipped_contracts_are_valid(example):
    contract = load_contract(EXAMPLES / example)
    assert contract is not None
    if contract.status == "supported":
        assert all(v in all_views() for v in contract.views)
        assert (EXAMPLES / example / contract.notebook).is_file()
    else:
        assert contract.reason


@pytest.mark.parametrize("bad, fragment", [
    ({**CONTRACT, "schema": "v0"}, "schema"),
    ({**CONTRACT, "views": ["nope.view"]}, "unknown view"),
    ({**CONTRACT, "outputs": {"data_dir": "../../etc"}}, "relative path"),
    ({**CONTRACT, "outputs": {"data_dir": "/abs"}}, "relative path"),
    ({"schema": "icesee-results/1", "model_name": "x", "status": "unsupported"}, "reason"),
    ({**CONTRACT, "views": []}, "at least one view"),
])
def test_invalid_contracts_are_rejected(bad, fragment):
    with pytest.raises(ContractError, match=fragment):
        parse_contract(bad)


def test_no_contract_means_none(tmp_path):
    assert load_contract(tmp_path) is None


# ── views ─────────────────────────────────────────────────────────────────
def test_lorenz_view_plots_the_notebook_series(tmp_path):
    d = tmp_path / "data"
    truth, members = write_run(d)
    ctx = open_context(d, parse_contract(CONTRACT))
    _r, fig = render_view(ctx, "lorenz.state_estimation", png=tmp_path / "s.png")
    data = figure_data(fig)
    assert [a["ylabel"] for a in data] == ["$x(t)$", "$y(t)$", "$z(t)$"]
    labels = [ln["label"] for ln in data[0]["lines"]]
    assert labels == ["True", "Background", "Observation", "Analysis"]
    np.testing.assert_array_equal(data[1]["lines"][0]["y"], truth[1])
    np.testing.assert_allclose(data[2]["lines"][3]["y"], members.mean(axis=1)[2])
    assert (tmp_path / "s.png").stat().st_size > 10_000


def test_error_and_spread_table(tmp_path):
    d = tmp_path / "data"
    write_run(d)
    ctx = open_context(d, parse_contract(CONTRACT), params={"enkf-parameters": {"vec_inputs": ["x", "y", "z"]}})
    rendered, fig = render_view(ctx, "da.error_and_spread")
    assert [r[0] for r in rendered.table["rows"]] == ["x", "y", "z"]
    assert rendered.supplementary and any("Supplementary" in n for n in rendered.notes)
    # background is truth + 1 everywhere -> RMSE exactly 1
    assert all(r[1] == "1" for r in rendered.table["rows"])


def _icepack_run(d: Path, nx=4, ny=2):
    hdim = (2 * nx + 1) * (2 * ny + 1)
    truth, members = write_run(d, model="icepack", nd=4 * hdim, nt=60, nens=3, nobs=4,
                               extra_tw={"Lxy": np.array([5000, 1200]), "nxy": np.array([nx, ny]),
                                         "b_io": np.array([200.0, -400.0])})
    xs, ys = np.meshgrid(np.linspace(0, 5000, 2 * nx + 1), np.linspace(0, 1200, 2 * ny + 1))
    with h5py.File(d / "dof_coordinates.h5", "w") as f:
        f["coordinates"] = np.column_stack([xs.ravel(), ys.ravel()])
    return truth, members, hdim


ICEPACK = {"schema": "icesee-results/1", "model_name": "icepack",
           "views": ["icepack.ensemble_profiles", "icepack.difference_maps", "icepack.snapshot_maps",
                     "icepack.transects", "icepack.smb_along_x", "icepack.truth_vs_wrong",
                     "icepack.error_tables"]}


@pytest.mark.parametrize("view", ICEPACK["views"])
def test_icepack_views_render(tmp_path, view):
    d = tmp_path / "data"
    _icepack_run(d)
    ctx = open_context(d, parse_contract(ICEPACK))
    options = {"h_indx": 3} if view in ("icepack.ensemble_profiles", "icepack.truth_vs_wrong") else {}
    if view in ("icepack.difference_maps", "icepack.error_tables"):
        options = {"years": [0.1, 0.5]}
    rendered, fig = render_view(ctx, view, options, png=tmp_path / "v.png")
    if get_view(view).kind == "table":
        assert rendered.table["rows"] and fig is None
    else:
        assert (tmp_path / "v.png").is_file()


def test_icepack_profiles_use_the_notebook_blocks(tmp_path):
    d = tmp_path / "data"
    truth, members, hdim = _icepack_run(d)
    ctx = open_context(d, parse_contract(ICEPACK))
    _r, fig = render_view(ctx, "icepack.ensemble_profiles", {"h_indx": 7})
    data = figure_data(fig)
    true_h = next(ln for ln in data[0]["lines"] if ln["label"] == "True")
    np.testing.assert_array_equal(true_h["y"], truth[7])
    true_u = next(ln for ln in data[1]["lines"] if ln["label"] == "True")
    np.testing.assert_allclose(true_u["y"], np.hypot(truth[hdim + 7], truth[2 * hdim + 7]))
    true_smb = next(ln for ln in data[2]["lines"] if ln["label"] == "True")
    np.testing.assert_array_equal(true_smb["y"], truth[3 * hdim + 7])
    assert sum(1 for ln in data[0]["lines"] if ln["color"] == "gray") == 3
    assert [a["ylabel"] for a in data] == ["Thickness (m)", "Velocity (m/year)", "SMB (m/year)"]


def test_maps_need_coordinates(tmp_path):
    d = tmp_path / "data"
    _icepack_run(d)
    (d / "dof_coordinates.h5").unlink()
    ctx = open_context(d, parse_contract(ICEPACK))
    with pytest.raises(ResultsError, match="mesh coordinates"):
        render_view(ctx, "icepack.snapshot_maps")
    render_view(ctx, "icepack.ensemble_profiles", {"h_indx": 1})   # time series still work


@pytest.mark.parametrize("options, fragment", [
    ({"h_indx": 10_000}, "at most"), ({"h_indx": -1}, "at least"), ({"h_indx": "x"}, "Grid point"),
    ({"bogus": 1}, "unknown option")])
def test_bad_options_are_invalid(tmp_path, options, fragment):
    d = tmp_path / "data"
    _icepack_run(d)
    ctx = open_context(d, parse_contract(ICEPACK))
    with pytest.raises(ResultsError, match=fragment) as err:
        render_view(ctx, "icepack.ensemble_profiles", options)
    assert err.value.kind == "invalid"


def test_issm_view_rebuilds_time_from_params(tmp_path):
    d = tmp_path / "data"
    nd = 6 * 20
    write_run(d, model="issm", nd=nd, nt=21, nens=3, nobs=4)
    params = {"modeling-parameters": {"tinitial": 0, "num_years": 10, "timesteps_per_year": 0.5},
              "enkf-parameters": {"freq_obs": 2, "obs_start_time": 2, "obs_max_time": 8}}
    contract = parse_contract({"schema": "icesee-results/1", "model_name": "issm",
                               "outputs": {"ensemble": "icesee_ensemble_data"},
                               "views": ["issm.profile_time_series"]})
    ctx = open_context(d, contract, params=params)
    rendered, fig = render_view(ctx, "issm.profile_time_series", {"h_indx": 5})
    data = figure_data(fig)
    assert [a["ylabel"] for a in data] == ["Thickness (m)", "Surface (m)", "Velocity (m/year)",
                                           "Bed topography", "Friction"]
    assert data[-1]["xlabel"] == "Time (years)"
    np.testing.assert_allclose(next(ln for ln in data[0]["lines"] if ln["label"] == "True")["x"],
                               np.linspace(0, 10, 21))
    ctx_noparams = open_context(d, contract, params=None)
    with pytest.raises(ResultsError, match="params"):
        render_view(ctx_noparams, "issm.profile_time_series")


# ── command line ──────────────────────────────────────────────────────────
def test_cli_writes_json_for_success_and_failure(tmp_path):
    d = tmp_path / "data"
    write_run(d)
    ex = tmp_path / "example"
    ex.mkdir()
    (ex / "results.yaml").write_text(json.dumps(CONTRACT))
    rc = render_main(["--data-dir", str(d), "--contract", str(ex), "--view", "lorenz.state_estimation",
                      "--png", str(tmp_path / "o.png"), "--json", str(tmp_path / "o.json")])
    assert rc == 0 and json.loads((tmp_path / "o.json").read_text())["ok"]
    (d / "icesee_ensemble_data.h5").unlink()
    rc = render_main(["--data-dir", str(d), "--contract", str(ex), "--view", "lorenz.state_estimation",
                      "--json", str(tmp_path / "e.json")])
    err = json.loads((tmp_path / "e.json").read_text())
    assert rc == 2 and err == {"ok": False, "kind": "missing", "view": "lorenz.state_estimation",
                               "message": err["message"]}
    rc = render_main(["--data-dir", str(d), "--contract", str(ex), "--view", "icepack.transects",
                      "--json", str(tmp_path / "u.json")])
    assert rc == 4


def test_user_views_replace_and_add_views(tmp_path):
    d = tmp_path / "data"
    write_run(d)
    mine = tmp_path / "mine"
    mine.mkdir()
    (mine / "results.yaml").write_text(json.dumps({**CONTRACT, "views": ["lorenz.state_estimation", "my.table"]}))
    (mine / "views.py").write_text(
        "from ICESEE.src.postprocessing.views import ViewResult, ViewSpec, get_view, register\n"
        "base = get_view('lorenz.state_estimation')\n"
        "register(ViewSpec(id=base.id, title='Mine', description='', source='me', build=base.build,"
        " options=base.options))\n"
        "register(ViewSpec(id='my.table', title='T', description='', source='me', kind='table',"
        " build=lambda ctx, o: ViewResult(table={'columns': ['n'], 'rows': [[ctx.outputs.n_ens]]})))\n")
    args = ["--data-dir", str(d), "--contract", str(mine), "--user-views", str(mine)]
    assert render_main([*args, "--describe", "--json", str(tmp_path / "d.json")]) == 0
    described = {v["id"]: v for v in json.loads((tmp_path / "d.json").read_text())["views"]}
    assert described["lorenz.state_estimation"]["customized"] and described["my.table"]["customized"]
    assert render_main([*args, "--view", "my.table", "--json", str(tmp_path / "t.json")]) == 0
    assert json.loads((tmp_path / "t.json").read_text())["table"]["rows"] == [[4]]
    (mine / "views.py").write_text("x = (\n")
    assert render_main([*args, "--describe", "--json", str(tmp_path / "e.json")]) == 5
    err = json.loads((tmp_path / "e.json").read_text())
    assert err["kind"] == "user_code" and "views.py line 1" in err["message"]
