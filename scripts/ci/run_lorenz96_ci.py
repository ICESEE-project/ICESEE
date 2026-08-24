from pathlib import Path
import os
import sys
import shutil
import subprocess
import tempfile

import matplotlib
matplotlib.use("Agg")

import matplotlib as mpl
import matplotlib.pyplot as plt
import h5py
import numpy as np
import yaml


REPO_ROOT = Path(__file__).resolve().parents[2]
LORENZ_DIR = REPO_ROOT / "applications" / "lorenz_model" / "examples" / "lorenz96"
CI_DIR = REPO_ROOT / "scripts" / "ci"
DATA_DIR = LORENZ_DIR / "_modelrun_datasets"
FIGURE_DIR = DATA_DIR / "figures"


def _mode_config(mode, temporary_directory):
    source = LORENZ_DIR / "params.yaml"
    with source.open("r", encoding="utf-8") as stream:
        config = yaml.safe_load(stream)
    output = DATA_DIR / f"ci_mode_{mode}"
    enkf = config["enkf-parameters"]
    enkf["execution_mode"] = mode
    enkf["data_path"] = str(output)
    enkf["restart_enabled"] = False
    enkf["force_fresh_start"] = True
    destination = Path(temporary_directory) / f"lorenz_mode_{mode}.yaml"
    with destination.open("w", encoding="utf-8") as stream:
        yaml.safe_dump(config, stream, sort_keys=False)
    return destination, output


def run_lorenz96_modes():
    print(f"Running Lorenz96 modes 0, 1, and 2 in {LORENZ_DIR}")
    outputs = {}
    with tempfile.TemporaryDirectory(prefix="icesee-lorenz-ci-") as temporary:
        for mode in (0, 1, 2):
            config, output = _mode_config(mode, temporary)
            command = [
                sys.executable, "-m",
                "ICESEE.applications.lorenz_model.examples.lorenz96.run_da_lorenz96",
                "-F", str(config),
            ]
            if mode in (1, 2):
                launcher = shutil.which("mpirun") or shutil.which("mpiexec")
                if launcher is None:
                    raise RuntimeError("MPI launcher is required for execution modes 1 and 2")
                command = [launcher, "-np", "2"] + command
            print(f"[ICESEE CI] mode {mode}: {' '.join(command)}")
            subprocess.run(command, cwd=LORENZ_DIR, check=True)
            outputs[mode] = output
    return outputs


def read_h5_dataset(file_path):
    data = {}
    with h5py.File(file_path, "r") as f:
        for key in f.keys():
            data[key] = f[key][:]
    return data


def load_lorenz_outputs(data_dir):

    tw_file = data_dir / "true-wrong-lorenz.h5"
    ensemble_file = data_dir / "icesee_ensemble_data.h5"
    true_nudged_file = data_dir / "true_nurged_states.h5"
    obs_file = data_dir / "synthetic_obs.h5"

    for path in [tw_file, ensemble_file, true_nudged_file, obs_file]:
        if not path.exists():
            raise FileNotFoundError(f"Missing expected Lorenz output file: {path}")

    tw = read_h5_dataset(tw_file)

    with h5py.File(ensemble_file, "r") as f:
        ensemble_vec_mean = f["ensemble_mean"][:]

    with h5py.File(true_nudged_file, "r") as f:
        ensemble_true_state = f["true_state"][:]
        ensemble_nurged_state = f["nurged_state"][:]

    with h5py.File(obs_file, "r") as f:
        if "hu_obs" in f:
            # Execution modes 0 and 1 retain the historical dense layout.
            w = f["hu_obs"][:]
        elif "hu_obs_compact" in f:
            # Execution mode 2 stores only observed state rows.  Reconstruct
            # the dense view expected by this small CI plotting routine; the
            # production analysis continues to consume the compact layout.
            compact = f["hu_obs_compact"][:]
            observed = np.asarray(f["obs_indices"][:], dtype=np.int64)
            active = np.asarray(f["obs_active"][:], dtype=bool)
            w = np.full(
                (ensemble_true_state.shape[0], compact.shape[1]),
                np.nan,
                dtype=compact.dtype,
            )
            compact = np.where(active, compact, np.nan)
            w[observed, :] = compact
        else:
            raise KeyError(
                f"{obs_file} contains neither 'hu_obs' nor 'hu_obs_compact'"
            )

    return {
        "t": tw["t"],
        "ind_m": tw["obs_index"],
        "tm_m": tw["obs_max_time"][0],
        "ensemble_true_state": ensemble_true_state,
        "ensemble_nurged_state": ensemble_nurged_state,
        "ensemble_vec_mean": ensemble_vec_mean,
        "w": w,
    }


def validate_mode_outputs(outputs):
    loaded = {mode: load_lorenz_outputs(path) for mode, path in outputs.items()}
    for mode, data in loaded.items():
        mean = np.asarray(data["ensemble_vec_mean"])
        if mean.size == 0 or not np.all(np.isfinite(mean)):
            raise AssertionError(f"execution mode {mode} produced an invalid ensemble mean")
        print(f"[ICESEE CI] mode {mode}: valid ensemble mean {mean.shape}")

    # Modes 1 and 2 implement the same stochastic filter and are the strict
    # numerical-parity pair. Mode 0 is a separate serial filter path and is a
    # smoke/finite-output check here.
    partial = np.asarray(loaded[1]["ensemble_vec_mean"])
    full = np.asarray(loaded[2]["ensemble_vec_mean"])
    if partial.shape != full.shape:
        raise AssertionError(f"mode 1/2 output shapes differ: {partial.shape} vs {full.shape}")
    max_difference = float(np.max(np.abs(partial - full)))
    if not np.allclose(partial, full, rtol=1e-7, atol=1e-7):
        raise AssertionError(f"mode 1/2 parity failed: max difference={max_difference:.6g}")
    print(f"[ICESEE CI] mode 1/2 parity max difference: {max_difference:.6g}")
    return loaded


def plot_lorenz_outputs(data):
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)

    mpl.rcParams["text.usetex"] = bool(shutil.which("latex"))
    mpl.rcParams["mathtext.fontset"] = "dejavusans"
    mpl.rcParams["font.family"] = "DejaVu Sans"

    font = {"family": "normal", "weight": "bold", "size": 14}
    mpl.rc("font", **font)

    t = data["t"]
    ind_m = data["ind_m"]
    tm_m = data["tm_m"]
    ensemble_true_state = data["ensemble_true_state"]
    ensemble_nurged_state = data["ensemble_nurged_state"]
    ensemble_vec_mean = data["ensemble_vec_mean"]
    w = data["w"]

    fig, ax = plt.subplots(nrows=3, ncols=1, figsize=(10, 8))
    ax = ax.flat

    labels = [r"$x(t)$", r"$y(t)$", r"$z(t)$"]

    for k in range(3):
        ax[k].plot(t, ensemble_true_state[k, :], label="True", linewidth=3)
        ax[k].plot(t, ensemble_nurged_state[k, :], ":", label="Background", linewidth=3)
        ax[k].plot(
            t[ind_m],
            w[k, :],
            "o",
            fillstyle="none",
            label="Observation",
            markersize=8,
            markeredgewidth=2,
        )
        ax[k].plot(t, ensemble_vec_mean[k, :], "--", label="Analysis", linewidth=3)
        ax[k].set_xlabel(r"$t$", fontsize=18)
        ax[k].set_ylabel(labels[k])
        ax[k].axvspan(0, tm_m, alpha=0.25, lw=0)

    ax[0].legend(loc="center", bbox_to_anchor=(0.5, 1.25), ncol=4, fontsize=13)
    fig.subplots_adjust(hspace=0.5)

    outfile = FIGURE_DIR / "lorenz96_ci.png"
    fig.savefig(outfile, dpi=200, bbox_inches="tight")
    plt.close(fig)

    print(f"Wrote {outfile}")


def main():
    parent = REPO_ROOT.parent
    if str(parent) not in sys.path:
        sys.path.insert(0, str(parent))

    outputs = run_lorenz96_modes()
    loaded = validate_mode_outputs(outputs)
    plot_lorenz_outputs(loaded[2])


if __name__ == "__main__":
    main()
