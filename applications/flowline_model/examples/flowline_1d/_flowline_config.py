# ==============================================================================
# @des: Static flowline_1d configuration shared by run_da_flowline.py (modes
#       0-2) and mode3_runner.py (mode 3), so both build the same time grid,
#       scales and grid from params.yaml.
# @date: 2026-10-21
# @author: Brian Kyanjo
# ==============================================================================

import numpy as np

from ICESEE.applications.flowline_model.examples.flowline_1d._flowline_model import initialize_model

# Terminus melt (m/yr) over the run: the truth ramps 5 -> 85, the wrong
# ("nurged") model and the forecast ensemble ramp 5 -> 45 (EnKF.ipynb).
FACEMELT_START = 5.0
FACEMELT_TRUE_END = 85.0
FACEMELT_WRONG_END = 45.0


def facemelt_schedule(end, **icesee_kwargs):
    """Terminus melt (m/s) at every model time index 0..nt."""
    return np.linspace(FACEMELT_START, end, int(icesee_kwargs["nt"]) + 1) / float(icesee_kwargs["year"])


def build_flowline_config(icesee_kwargs):
    """Add flowline_1d's derived model values to ``icesee_kwargs`` in place.

    One model step is ``year / timesteps_per_year`` seconds and the run has
    ``num_years * timesteps_per_year`` steps -- the same ``nt`` and ``t``
    (in years) the ICESEE loader uses for the observation schedule.
    Glen's exponent is stored as ``glen_n``: the DA loops reuse ``n``
    for the number of time steps. Calling it again on a built config is a
    no-op.
    """
    if icesee_kwargs.get("flowline_config_built"):
        # run_da_flowline.py builds it before dispatching to mode3_runner;
        # building twice would rescale accum/facemelt by the year again.
        return icesee_kwargs
    num_years = float(icesee_kwargs["num_years"])
    steps_per_year = float(icesee_kwargs.get("timesteps_per_year", 1))
    nt = int(round(num_years * steps_per_year))
    year = float(icesee_kwargs["year"])
    glen_n = int(icesee_kwargs["n"])
    m = 1 / glen_n

    icesee_kwargs.update({
        "nd": int(float(icesee_kwargs["num_state_vars"])),
        "nt": nt,
        "NT": nt,
        "seed": float(icesee_kwargs["seed"]),
        "t": np.linspace(0, num_years, nt + 1),
        "hscale": float(icesee_kwargs["hscale"]),
        "A": float(icesee_kwargs["A"]),
        "n": glen_n,
        "glen_n": glen_n,
        "C": float(icesee_kwargs["C"]),
        "rho_ice": float(icesee_kwargs["rho_ice"]),
        "rho_water": float(icesee_kwargs["rho_water"]),
        "g": float(icesee_kwargs["g"]),
        "accum": float(icesee_kwargs["accum"]) / year,
        "facemelt": float(icesee_kwargs["facemelt"]) / year,
        "m": m,
        "B": float(icesee_kwargs["A"]) ** (-1 / glen_n),
        "ascale": 1.0 / year,
        "N1": int(icesee_kwargs["N1"]),
        "N2": int(icesee_kwargs["N2"]),
        "NX": int(icesee_kwargs["N1"]) + int(icesee_kwargs["N2"]),
        "TF": year / steps_per_year,
        "sigGZ": float(icesee_kwargs["sigGZ"]),
        "sillamp": float(icesee_kwargs["sillamp"]),
        "sillsmooth": float(icesee_kwargs["sillsmooth"]),
        "xsill": float(icesee_kwargs["xsill"]),
        "tcurrent": int(icesee_kwargs["tcurrent"]),
        "transient": int(icesee_kwargs["transient"]),
        "scalar_inputs": icesee_kwargs.get("scalar_inputs", []),
    })
    icesee_kwargs["sigma1"] = np.linspace(icesee_kwargs["sigGZ"] / (icesee_kwargs["N1"] + 0.5),
                                          icesee_kwargs["sigGZ"], icesee_kwargs["N1"])
    icesee_kwargs["sigma2"] = np.linspace(icesee_kwargs["sigGZ"], 1, icesee_kwargs["N2"] + 1)
    # velocity scale balancing basal drag with the driving stress
    icesee_kwargs["uscale"] = (icesee_kwargs["rho_ice"] * icesee_kwargs["g"] * icesee_kwargs["hscale"]
                               * icesee_kwargs["ascale"] / icesee_kwargs["C"]) ** (1 / (m + 1))

    xscale = icesee_kwargs["uscale"] * icesee_kwargs["hscale"] / icesee_kwargs["ascale"]
    sigma = np.concatenate((icesee_kwargs["sigma1"], icesee_kwargs["sigma2"][1:icesee_kwargs["N2"] + 1]))
    icesee_kwargs.update({
        "xscale": xscale,
        "tscale": xscale / icesee_kwargs["uscale"],
        "eps": icesee_kwargs["B"] * ((icesee_kwargs["uscale"] / xscale) ** (1 / glen_n))
               / (2 * icesee_kwargs["rho_ice"] * icesee_kwargs["g"] * icesee_kwargs["hscale"]),
        "lambda": 1 - (icesee_kwargs["rho_ice"] / icesee_kwargs["rho_water"]),
        "dt": icesee_kwargs["TF"],
        "sigma": sigma,
        "grid": {"sigma": sigma,
                 "sigma_elem": np.concatenate(([0], (sigma[:-1] + sigma[1:]) / 2)),
                 "dsigma": np.diff(sigma)},
    })

    huxg_out0 = initialize_model(**icesee_kwargs)
    icesee_kwargs["nd"] = int(huxg_out0.shape[0])
    icesee_kwargs["flowline_config_built"] = True
    icesee_kwargs["var_nd"] = {var: (1 if var in icesee_kwargs["scalar_inputs"] else icesee_kwargs["NX"])
                               for var in icesee_kwargs["vec_inputs"]}
    return icesee_kwargs
