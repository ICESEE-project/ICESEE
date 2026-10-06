# ==============================================================================
# @des: Production execution_mode 3 DA-cycle runner for the flowline_1d
#       application. Registers itself with
#       ``src/parallelization/distributed_mode3_registry.py`` at import time.
# @date: 2026-08-26
# @author: Brian Kyanjo
# ==============================================================================
"""Production ``execution_mode: 3`` DA-cycle runner for flowline_1d.

This module wires ``FLOWLINE_NATIVE_ADAPTER``
(``applications/flowline_model/examples/flowline_1d/_flowline_native.py``)
into ``src/parallelization/distributed_native_cycle.py`` to make mode 3
functional for flowline_1d, and registers itself with
``src/parallelization/distributed_mode3_registry.py`` so
``icesee_da_distributed.py`` can dispatch to it instead of raising
``NotImplementedError``. Structure mirrors the icepack/ISSM mode-3 runners
closely; the differences below are specific to flowline_1d.

Static config is rank-independent -- no root-rebuild-and-broadcast needed
--------------------------------------------------------------------------
Unlike idealized_pig (whose reference ``nd`` depends on a Firedrake mesh
partitioned by ``comm``) and ISSM (whose per-rank config depends on a
per-rank MATLAB server), flowline_1d's static config
(``_build_static_config`` below, mirroring ``run_da_flowline.py``'s config
block) is pure float/int arithmetic over ``params.yaml`` values plus one
deterministic ``scipy.optimize.root``/JAX-autodiff solve with no randomness
of any kind. It therefore produces bit-identical results on every rank
without any MPI dependency, so this runner calls it directly on every rank
-- no ``comm=COMM_SELF`` rebuild, no broadcast.

Setup-phase reuse
------------------
Generating the true/nurged trajectories and synthetic observations reuses
the exact same model-agnostic, single-rank generators execution mode 0
uses (``src/EnKF/_generate_true_wrong_state.py``,
``src/EnKF/_generate_synthetic_observations.py``), called only from the
real MPI world root, exactly like the other three mode-3 runners --
these generators hardcode a single-rank identity internally and would race
or duplicate-write if called from every rank.

Per-member initialization reuses ``initialize_ensemble`` (the same
per-member init function modes 1/2 already use), mirroring ISSM's native
adapter rather than idealized_pig's (which reconstructs its own initial
state independently of the setup-phase HDF5 artifacts).

Observation-error parity
-------------------------
flowline_1d's default ``enkf_observation_error_mode`` (from
``use_ensemble_pertubations: true`` in ``params.yaml``) is
``"legacy_prior_anomalies"``, exactly matching modes 1/2. ``"stochastic_R"``
is also supported, mirroring the other three mode-3 runners.
``"generated_R"`` is not supported here; selecting it raises
``NotImplementedError``.

No inversion
-------------
flowline_1d has no inversion mechanism (confirmed by inspection; see
``_flowline_native.py``'s docstring), so unlike ISSM's runner this one
never sets ``apply_inversion=True`` and calls
``run_native_global_analysis_cycle`` with its default (no inversion).

Output writing
---------------
flowline_1d's per-member state is small (``2*NX+1``), but this runner still
uses the mode-3-native, rank-sharded, no-gather checkpoint primitive
(``src/parallelization/distributed_checkpoint.py``,
``save_distributed_checkpoint``/``load_distributed_checkpoint``) rather than
gathering to root and writing modes 0-2's ``icesee_ensemble_data.h5``
schema, for consistency with idealized_pig's and ISSM's mode-3 output and
to avoid privileging a state-size-dependent code path. This is an honest
schema divergence from modes 0-2, not a compatible drop-in replacement --
downstream plotting/post-processing built against that schema will not read
mode-3 output directly. Directory layout under
``<data_path>/_mode3_state_history/``:

- ``initial_condition/step_00000000/`` -- the native pool's initial
  (pre-forecast) ensemble, written once before the timestep loop.
- ``steps/step_<k:08d>/`` for ``k`` in ``[0, nt)`` -- the analyzed (or, on
  non-observed steps, forecast-only) ensemble *after* native cycle ``k``,
  i.e. the mode-3 equivalent of modes 0-2's ensemble column ``k + 1``.
"""

from __future__ import annotations

import os

import h5py
import numpy as np
from mpi4py import MPI

from ICESEE.applications.flowline_model.examples.flowline_1d._flowline_model import (
    initialize_model,
)
from ICESEE.applications.flowline_model.examples.flowline_1d._flowline_native import (
    FLOWLINE_NATIVE_ADAPTER,
)
from ICESEE.applications.supported_models import SupportedModels
from ICESEE.src.EnKF._generate_synthetic_observations import (
    generate_synthetic_observations,
)
from ICESEE.src.EnKF._generate_true_wrong_state import generate_true_wrong_state
from ICESEE.src.parallelization.distributed_checkpoint import (
    save_distributed_checkpoint,
)
from ICESEE.src.parallelization.distributed_mode3_registry import (
    register_execution_mode_3,
)
from ICESEE.src.parallelization.distributed_native_cycle import (
    NativeObservationBatch,
    run_native_global_analysis_cycle,
)
from ICESEE.src.parallelization.distributed_native_runtime import (
    initialize_native_member_pool,
)
from ICESEE.src.parallelization.distributed_topology import (
    create_distributed_topology,
    register_topology_run_metadata,
)
from ICESEE.src.utils.icesee_context import (
    normalize_execution_mode,
    normalize_icesee_kwargs,
)
from ICESEE.src.utils.localization import active_observation_std
from ICESEE.src.utils.tools import save_all_data
from ICESEE.src.utils.performance import emit_performance_report, register_run_metadata
from ICESEE.src.utils.utils import UtilsFunctions

_SUPPORTED_ERROR_MODES = {"legacy_prior_anomalies", "stochastic_r"}
_DEFAULT_RUN_ID = "flowline-1d-mode3"


def _build_static_config(icesee_kwargs: dict) -> dict:
    """Rebuild flowline_1d's static, rank-independent model config in place.

    Mirrors ``run_da_flowline.py``'s config block exactly (see module
    docstring for why this has no MPI/comm dependency and can safely run
    identically on every rank).
    """

    icesee_kwargs.update({"nd": int(float(icesee_kwargs["num_state_vars"]))})

    nt = int(float(icesee_kwargs["num_years"]))
    icesee_kwargs.update(
        {
            "nt": nt,
            "NT": nt,
            "nd": icesee_kwargs["nd"],
            "seed": float(icesee_kwargs["seed"]),
            "t": np.linspace(0, nt, nt + 1),
            "hscale": float(icesee_kwargs["hscale"]),
            "A": float(icesee_kwargs["A"]),
            "n": int(icesee_kwargs["n"]),
            "C": float(icesee_kwargs["C"]),
            "rho_ice": float(icesee_kwargs["rho_ice"]),
            "rho_water": float(icesee_kwargs["rho_water"]),
            "g": float(icesee_kwargs["g"]),
            "accum": float(icesee_kwargs["accum"]) / float(icesee_kwargs["year"]),
            "facemelt": float(icesee_kwargs["facemelt"]) / float(icesee_kwargs["year"]),
            "m": 1 / int(icesee_kwargs["n"]),
            "B": float(icesee_kwargs["A"]) ** (-1 / int(icesee_kwargs["n"])),
            "ascale": 1.0 / float(icesee_kwargs["year"]),
            "N1": int(icesee_kwargs["N1"]),
            "N2": int(icesee_kwargs["N2"]),
            "NX": int(icesee_kwargs["N1"]) + int(icesee_kwargs["N2"]),
            "TF": float(icesee_kwargs["year"]),
            "sigGZ": float(icesee_kwargs["sigGZ"]),
            "sigma1": np.linspace(
                float(icesee_kwargs["sigGZ"]) / (int(icesee_kwargs["N1"]) + 0.5),
                float(icesee_kwargs["sigGZ"]),
                int(icesee_kwargs["N1"]),
            ),
            "sigma2": np.linspace(
                float(icesee_kwargs["sigGZ"]), 1, int(icesee_kwargs["N2"]) + 1
            ),
            "sillamp": float(icesee_kwargs["sillamp"]),
            "sillsmooth": float(icesee_kwargs["sillsmooth"]),
            "xsill": float(icesee_kwargs["xsill"]),
            "tcurrent": int(icesee_kwargs["tcurrent"]),
            "transient": int(icesee_kwargs["transient"]),
            "uscale": float(icesee_kwargs["rho_ice"])
            * float(icesee_kwargs["g"])
            * float(icesee_kwargs["hscale"])
            * (1.0 / float(icesee_kwargs["year"]))
            / float(icesee_kwargs["C"]),
            "scalar_inputs": icesee_kwargs.get("scalar_inputs", []),
        }
    )

    xscale = icesee_kwargs["uscale"] * icesee_kwargs["hscale"] / icesee_kwargs["ascale"]
    sigma = np.concatenate(
        (icesee_kwargs["sigma1"], icesee_kwargs["sigma2"][1 : icesee_kwargs["N2"] + 1])
    )
    icesee_kwargs.update(
        {
            "xscale": xscale,
            "tscale": xscale / icesee_kwargs["uscale"],
            "eps": icesee_kwargs["B"]
            * ((icesee_kwargs["uscale"] / xscale) ** (1 / icesee_kwargs["n"]))
            / (2 * icesee_kwargs["rho_ice"] * icesee_kwargs["g"] * icesee_kwargs["hscale"]),
            "lambda": 1 - (icesee_kwargs["rho_ice"] / icesee_kwargs["rho_water"]),
            "dt": icesee_kwargs["TF"] / nt,
            "sigma": sigma,
            "grid": {
                "sigma": sigma,
                "sigma_elem": np.concatenate(([0], (sigma[:-1] + sigma[1:]) / 2)),
                "dsigma": np.diff(sigma),
            },
        }
    )

    huxg_out0 = initialize_model(**icesee_kwargs)
    icesee_kwargs["nd"] = int(huxg_out0.shape[0])
    var_nd = {
        var: (1 if var in icesee_kwargs["scalar_inputs"] else icesee_kwargs["NX"])
        for var in icesee_kwargs["vec_inputs"]
    }
    icesee_kwargs["var_nd"] = var_nd
    return icesee_kwargs


def run_flowline_execution_mode_3(**icesee_kwargs):
    """Drive flowline_1d's full mode-3 DA cycle. Registered under ``"flowline"``."""

    icesee_kwargs = normalize_icesee_kwargs(icesee_kwargs)
    normalize_execution_mode(icesee_kwargs, expected=3)

    world = MPI.COMM_WORLD
    world_rank = int(world.Get_rank())

    global_start_time = MPI.Wtime()
    time_forecast_step = 0.0
    time_forecast_file_writing = 0.0

    model = icesee_kwargs.get("model_name")
    Nens = int(icesee_kwargs["Nens"])

    model_module = SupportedModels(
        model=model, verbose=icesee_kwargs.get("verbose")
    ).call_model()
    icesee_kwargs.update(
        {
            "model_module": model_module,
            "vec_inputs_old": icesee_kwargs.get("vec_inputs"),
            "model_nprocs": icesee_kwargs.get("model_nprocs") or 1,
        }
    )
    icesee_kwargs["observed_vars_params"] = list(
        icesee_kwargs.get("observed_vars", [])
    ) + list(icesee_kwargs.get("observed_params", []))
    icesee_kwargs["all_observed"] = icesee_kwargs["observed_vars_params"]

    # --- static config: rank-independent, no broadcast needed (see module
    # docstring). ---
    icesee_kwargs = _build_static_config(icesee_kwargs)
    nt = int(icesee_kwargs["nt"])
    nd = int(icesee_kwargs["nd"])

    _modelrun_datasets = icesee_kwargs.get("data_path") or "_modelrun_datasets"
    os.environ["ICESEE_RESULTS_DIR"] = str(_modelrun_datasets)
    if world_rank == 0 and not os.path.exists(_modelrun_datasets):
        os.makedirs(_modelrun_datasets, exist_ok=True)
    world.Barrier()

    _true_nurged = f"{_modelrun_datasets}/true_nurged_states.h5"
    _synthetic_obs = f"{_modelrun_datasets}/synthetic_obs.h5"
    icesee_kwargs.update(
        {"true_nurged_file": _true_nurged, "synthetic_obs_file": _synthetic_obs}
    )

    # --- setup phase: real world root only (see module docstring -- the
    # two generic generators hardcode a single-rank identity and would race
    # or duplicate-write if run from every rank). Non-root ranks keep this
    # timing at 0.0; the MPI.MAX reduction below picks up root's real
    # value, matching modes 0/1/2's timing-aggregation pattern. ---
    true_wrong_time = 0.0
    observation_time = 0.0
    if world_rank == 0:
        _t = MPI.Wtime()
        icesee_kwargs = generate_true_wrong_state(**icesee_kwargs)
        true_wrong_time = MPI.Wtime() - _t
        _t = MPI.Wtime()
        icesee_kwargs = generate_synthetic_observations(**icesee_kwargs)
        observation_time = MPI.Wtime() - _t
    world.Barrier()

    with h5py.File(_synthetic_obs, "r") as f:
        hu_obs = f["hu_obs"][:]
        error_R = f["R"][:]
    icesee_kwargs["error_R"] = error_R

    topology = create_distributed_topology(world, spatial_ranks=1)
    adapter = FLOWLINE_NATIVE_ADAPTER
    pool = initialize_native_member_pool(adapter, topology, icesee_kwargs)

    obs_indices = np.asarray(
        UtilsFunctions(icesee_kwargs).JObs_indices(nd), dtype=np.int64
    )
    error_mode = str(
        icesee_kwargs.get("enkf_observation_error_mode", "legacy_prior_anomalies")
    ).lower()
    if error_mode not in _SUPPORTED_ERROR_MODES:
        raise NotImplementedError(
            "execution_mode 3's flowline runner supports "
            "'legacy_prior_anomalies' and 'stochastic_R' observation-error "
            f"modes; got {error_mode!r}. See the module docstring for why "
            "'generated_R' is not (yet) supported here."
        )
    base_seed = int(icesee_kwargs.get("base_seed", 42))

    obs_t, ind_m, m_obs = UtilsFunctions(icesee_kwargs).generate_observation_schedule(
        **icesee_kwargs
    )
    icesee_kwargs.update(
        {"obs_t": obs_t, "obs_index": ind_m, "number_obs_instants": m_obs, "m_obs": m_obs}
    )
    km = 0

    run_id = str(icesee_kwargs.get("run_id", _DEFAULT_RUN_ID))
    history_root = os.path.join(_modelrun_datasets, "_mode3_state_history")
    initial_root = os.path.join(history_root, "initial_condition")
    steps_root = os.path.join(history_root, "steps")

    # --- initial (pre-forecast) native ensemble, written once before the
    # timestep loop -- the mode-3 equivalent of modes 0-2's ensemble column
    # 0 (see module docstring for the checkpoint directory layout). ---
    init_file_time = 0.0
    _t = MPI.Wtime()
    save_distributed_checkpoint(
        initial_root,
        0,
        pool.snapshot_owned(),
        topology,
        run_id=run_id,
        metadata={"Nens": Nens, "description": "initial (pre-forecast) ensemble"},
    )
    init_file_time = MPI.Wtime() - _t

    for k in range(nt):
        do_analysis = bool(km < m_obs and k == ind_m[km])
        batches = []
        if do_analysis:
            y_full = np.asarray(hu_obs[:, km], dtype=np.float64).copy()
            y_full[np.isnan(y_full)] = 0.0
            d = y_full[obs_indices]
            member_errors = None
            if error_mode == "stochastic_r":
                sigma = active_observation_std(icesee_kwargs, km, obs_indices)
                obs_seed = base_seed + 1000003 * (km + 1)
                rng = np.random.default_rng(obs_seed)
                eta = rng.standard_normal((obs_indices.size, Nens)) * sigma[:, None]
                eta -= np.mean(eta, axis=1, keepdims=True)
                member_errors = {
                    member_id: eta[:, member_id] for member_id in pool.member_ids
                }
            batches.append(
                NativeObservationBatch(
                    observation_ids=obs_indices,
                    values=d,
                    member_errors=member_errors,
                )
            )

        # run_native_global_analysis_cycle fuses the forecast and (when
        # scheduled) analysis update into one call at this API boundary, so
        # unlike modes 0-2 the two cannot be timed separately here; the
        # combined per-step cost is reported as forecast_step_time.
        _t = MPI.Wtime()
        result = run_native_global_analysis_cycle(
            pool,
            adapter,
            k,
            batches,
            number_of_batches=(1 if do_analysis else 0),
            topology=topology,
            icesee_kwargs=icesee_kwargs,
            error_mode=error_mode,
        )
        time_forecast_step += MPI.Wtime() - _t
        if do_analysis:
            km += 1

        _t = MPI.Wtime()
        save_distributed_checkpoint(
            steps_root,
            k,
            result.local_analysis,
            topology,
            run_id=run_id,
            metadata={
                "Nens": Nens,
                "observation_rows": result.observation_rows,
                "error_mode": error_mode,
                "did_analysis": do_analysis,
            },
        )
        time_forecast_file_writing += MPI.Wtime() - _t

    world.Barrier()

    analysis_file_time = 0.0
    if world_rank == 0:
        Lx = icesee_kwargs.get("Lx", 1.0)
        Ly = icesee_kwargs.get("Ly", 1.0)
        nx = icesee_kwargs.get("nx", 1)
        ny = icesee_kwargs.get("ny", 1)
        b_in = icesee_kwargs.get("b_in", 0.0)
        b_out = icesee_kwargs.get("b_out", 0.0)
        _t = MPI.Wtime()
        save_all_data(
            icesee_kwargs,
            data_path=_modelrun_datasets,
            nofilter=True,
            t=icesee_kwargs["t"],
            b_io=np.array([b_in, b_out]),
            Lxy=np.array([Lx, Ly]),
            nxy=np.array([nx, ny]),
            obs_max_time=np.array([icesee_kwargs["obs_max_time"]]),
            obs_index=icesee_kwargs["obs_index"],
            run_mode=np.array([icesee_kwargs["execution_flag"]]),
        )
        analysis_file_time = MPI.Wtime() - _t

    # ─────────────────────────────────────────────────────────────
    #  End Timer and Aggregate Elapsed Time Across Ranks
    # ─────────────────────────────────────────────────────────────
    global_elapsed_time = MPI.Wtime() - global_start_time

    register_run_metadata(
        execution_mode=icesee_kwargs.get("execution_mode"),
        model=icesee_kwargs.get("model_name"),
        forecast_steps=nt,
        analysis_events=km,
    )
    register_topology_run_metadata(topology, Nens)
    # Analyses are fused into the forecast step, so there is no separate
    # analysis timer: not measured when analyses ran, no events otherwise.
    emit_performance_report(
        world,
        elapsed_s=global_elapsed_time,
        phases={"true_wrong_state": true_wrong_time, "observation_generation": observation_time, "forecast_step": time_forecast_step, "init_file_io": init_file_time, "forecast_file_io": time_forecast_file_writing, "analysis_file_io": analysis_file_time, "analysis_step": None if km else 0.0},
        counts={"forecast_step": nt, "analysis_step": km},
        output_dir=_modelrun_datasets,
    )

    return icesee_kwargs


register_execution_mode_3(
    "flowline",
    run_flowline_execution_mode_3,
    notes=(
        "Production mode-3 DA-cycle runner for flowline_1d using "
        "FLOWLINE_NATIVE_ADAPTER (_flowline_native.py) and "
        "src/parallelization/distributed_native_cycle.py. Per-timestep "
        "output uses rank-sharded distributed checkpoints (no gather), not "
        "modes 0-2's icesee_ensemble_data.h5 schema -- see this module's "
        "docstring for the checkpoint directory layout and other "
        "divergences. No inversion mechanism exists for this application."
    ),
)
