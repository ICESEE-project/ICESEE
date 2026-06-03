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


REPO_ROOT = Path(__file__).resolve().parents[1]
LORENZ_DIR = REPO_ROOT / "applications" / "lorenz_model" / "examples" / "lorenz96"
FIGURE_DIR = REPO_ROOT / "docs" / "figures"


def run_lorenz96_example():
    print("Running Lorenz96 example...")
    print(f"Working directory: {LORENZ_DIR}")

    subprocess.run(
        [
            sys.executable,
            "-m",
            "ICESEE.applications.lorenz_model.examples.lorenz96.run_da_lorenz96",
        ],
        cwd=LORENZ_DIR,
        check=True,
    )


def read_h5_dataset(file_path):
    data = {}
    with h5py.File(file_path, "r") as f:
        for key in f.keys():
            data[key] = f[key][:]
    return data


def load_lorenz_outputs():
    results_dir = LORENZ_DIR / "results"
    data_dir = LORENZ_DIR / "_modelrun_datasets"

    tw_file = results_dir / "true-wrong-lorenz.h5"
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
        w = f["hu_obs"][:]

    return {
        "t": tw["t"],
        "ind_m": tw["obs_index"],
        "tm_m": tw["obs_max_time"][0],
        "ensemble_true_state": ensemble_true_state,
        "ensemble_nurged_state": ensemble_nurged_state,
        "ensemble_vec_mean": ensemble_vec_mean,
        "w": w,
    }


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

    run_lorenz96_example()
    data = load_lorenz_outputs()
    plot_lorenz_outputs(data)


if __name__ == "__main__":
    main()