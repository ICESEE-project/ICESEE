from pathlib import Path
import os
import sys
import shutil
import subprocess

import matplotlib
matplotlib.use("Agg")

import matplotlib as mpl
import matplotlib.pyplot as plt
import h5py

# from ICESEE.config._utility_imports import extract_datasets_from_h5
from ICESEE.src.utils.tools import extract_datasets_from_h5


def add_icesee_parent_to_path():
    repo = Path.cwd()
    parent = repo.parent

    if str(parent) not in sys.path:
        sys.path.insert(0, str(parent))


def run_lorenz96_example():
    print("Running Lorenz96 example...")

    subprocess.run(
        [
            sys.executable,
            "-m",
            "ICESEE.applications.lorenz_model.examples.lorenz96.run_da_lorenz96",
        ],
        check=True,
    )


def load_lorenz_outputs():
    data_dir = Path("_modelrun_datasets")
    results_dir = Path("results")

    model_name = "lorenz"
    filter_type = "true-wrong"

    tw_file = results_dir / f"{filter_type}-{model_name}.h5"
    if not tw_file.exists():
        raise FileNotFoundError(f"Missing expected Lorenz truth file: {tw_file}")

    datasets_tw = extract_datasets_from_h5(str(tw_file))

    t = datasets_tw["t"]
    ind_m = datasets_tw["obs_index"]
    tm_m = datasets_tw["obs_max_time"][0]
    run_mode = datasets_tw["run_mode"][0]

    filter_type = "EnTKF"

    if run_mode != 0:
        da_file = results_dir / f"{filter_type}-{model_name}.h5"
        if not da_file.exists():
            raise FileNotFoundError(f"Missing expected DA result file: {da_file}")

        datasets = extract_datasets_from_h5(str(da_file))
        ensemble_vec_mean = datasets["ensemble_vec_mean"]
        ensemble_bg = datasets["ensemble_bg"]
    else:
        ensemble_file = data_dir / "icesee_ensemble_data.h5"
        if not ensemble_file.exists():
            raise FileNotFoundError(f"Missing expected ensemble file: {ensemble_file}")

        with h5py.File(ensemble_file, "r") as f:
            ensemble_vec_mean = f["ensemble_mean"][:]
            ensemble_bg = None

    true_nudged_file = data_dir / "true_nurged_states.h5"
    if not true_nudged_file.exists():
        raise FileNotFoundError(f"Missing expected true/nudged state file: {true_nudged_file}")

    with h5py.File(true_nudged_file, "r") as f:
        ensemble_true_state = f["true_state"][:]
        ensemble_nurged_state = f["nurged_state"][:]

    obs_file = data_dir / "synthetic_obs.h5"
    if not obs_file.exists():
        raise FileNotFoundError(f"Missing expected observation file: {obs_file}")

    with h5py.File(obs_file, "r") as f:
        w = f["hu_obs"][:]

    return {
        "t": t,
        "ind_m": ind_m,
        "tm_m": tm_m,
        "ensemble_true_state": ensemble_true_state,
        "ensemble_nurged_state": ensemble_nurged_state,
        "ensemble_vec_mean": ensemble_vec_mean,
        "ensemble_bg": ensemble_bg,
        "w": w,
    }


def plot_lorenz_outputs(data):
    outdir = Path("docs/figures")
    outdir.mkdir(parents=True, exist_ok=True)

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

    outfile = outdir / "lorenz96_ci.png"
    fig.savefig(outfile, dpi=200, bbox_inches="tight")
    plt.close(fig)

    print(f"Wrote {outfile}")


def main():
    add_icesee_parent_to_path()
    run_lorenz96_example()
    data = load_lorenz_outputs()
    plot_lorenz_outputs(data)


if __name__ == "__main__":
    main()