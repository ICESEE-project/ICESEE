# ==============================================================================
# @des: Production execution_mode 3 DA-cycle runner for the Lorenz96
#       application. Registers itself with
#       ``src/parallelization/distributed_mode3_registry.py`` at import time.
# @date: 2026-08-24
# @author: Brian Kyanjo
# ==============================================================================
"""Production ``execution_mode: 3`` DA-cycle runner for Lorenz96.

This module wires ``LorenzNativeAdapter``
(``applications/lorenz_model/lorenz_utils/distributed_native_adapter.py``)
into ``src/parallelization/distributed_native_cycle.py`` to make mode 3
actually functional for the Lorenz96 application, and registers itself with
``src/parallelization/distributed_mode3_registry.py`` so
``icesee_da_distributed.py`` can dispatch to it instead of raising
``NotImplementedError``.

Setup-phase reuse
------------------
Generating the true/nurged states, synthetic observations, and initial
ensemble is *not* reimplemented here.  The exact same single-rank generators
used by execution mode 0 (``src/EnKF/_generate_true_wrong_state.py``,
``src/EnKF/_generate_synthetic_observations.py``,
``src/EnKF/_ensemble_initialization.py``) are called, but only from the real
MPI world root -- these three functions hardcode ``rank_world = 0`` and
``size_world = 1`` internally and write shared HDF5 artifacts
(``true_nurged_states.h5``, ``synthetic_obs.h5``, ``icesee_ensemble_data.h5``)
under ``data_path``, so calling them concurrently from every mode-3 rank
would race on the same files.  After an ``MPI.COMM_WORLD.Barrier()``, every
rank reads the two artifacts it needs back from disk (cheap for Lorenz's
``nd=3`` state and ``Nens`` on the order of tens of members).

Observation-error parity
-------------------------
The Lorenz96 example's default ``enkf_observation_error_mode`` (set by
``config/_utility_imports.py`` from ``use_ensemble_pertubations: true`` in
``params.yaml``) is ``"legacy_prior_anomalies"``, exactly matching modes
1/2's ``EnKF_X5`` branch in
``src/parallelization/_mpi_analysis_functions.py``: the observation-error
term is the forecast ensemble's own anomalies
(``eta = HA - mean(HA, axis=1)``), not independently drawn noise. This is
already implemented identically inside
``StochasticAnalysisProducts.add_chunk`` (``src/parallelization/
distributed_analysis.py``), so no ``member_errors`` need to be supplied for
this mode -- it is deterministic given the forecast, and therefore
automatically rank/order independent.

``"stochastic_R"`` is also supported for completeness (should a user select
it): ``active_observation_std``/``stochastic_observation_terms``
(``src/utils/localization.py``) show that its perturbed-observation
``eta`` depends only on ``(shape, sigma, seed)`` -- not on which rank
computes it -- using
``obs_seed = base_seed + 1000003 * (k_obs + 1)`` (the exact formula
``EnKF_X5`` uses). Every rank can therefore redundantly compute the full
``(len(obs_indices), Nens)`` eta matrix (trivially cheap at Lorenz's scale)
and slice out only the columns for its locally scheduled members, giving
bit-identical parity with modes 1/2 -- unlike the AR(1) process-noise case
in ``distributed_native_adapter.py``, which cannot be made bit-identical
under genuine per-member parallelism.

``"generated_R"`` is not supported by this runner (Lorenz96 does not use it
by default); selecting it raises ``NotImplementedError``.

Output writing
--------------
Lorenz96's ``nd=3, Nens`` on the order of tens is trivially small regardless
of mode 3's bounded-*history* invariant (which concerns per-member *time
history*, not one single-timestep full-ensemble snapshot). At the end of
every timestep this runner gathers each rank's locally owned analyzed member
state to the world root and writes one column into the same
``icesee_ensemble_data.h5`` schema modes 0-2 already produce, preserving
downstream plotting/post-processing compatibility.
"""

from __future__ import annotations

import os

import h5py
import numpy as np
from mpi4py import MPI

from ICESEE.applications.lorenz_model.lorenz_utils.distributed_native_adapter import (
    LorenzNativeAdapter,
)
from ICESEE.applications.supported_models import SupportedModels
from ICESEE.src.EnKF._ensemble_initialization import ensemble_initialization
from ICESEE.src.EnKF._generate_synthetic_observations import (
    generate_synthetic_observations,
)
from ICESEE.src.EnKF._generate_true_wrong_state import generate_true_wrong_state
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


def _write_ensemble_column(
    topology,
    local_members: dict,
    column: int,
    Nens: int,
    nd: int,
    ensemble_data_path: str,
) -> None:
    """Gather one timestep's analyzed members to the world root and write it.

    Only ever holds one ``(nd, Nens)`` column in memory (on the root rank),
    never the full ``(nd, Nens, nt+1)`` history -- consistent with mode 3's
    bounded-per-rank-memory invariant.
    """

    gathered = topology.world.gather(local_members, root=0)
    if topology.world_rank != 0:
        return
    full = np.zeros((nd, Nens), dtype=np.float64)
    for rank_members in gathered:
        for member_id, values in rank_members.items():
            full[:, int(member_id)] = np.asarray(values, dtype=np.float64)
    with h5py.File(ensemble_data_path, "r+") as f:
        f["ensemble"][:, :, column] = full
        f["ensemble_mean"][:, column] = full.mean(axis=1)


def run_lorenz96_execution_mode_3(**icesee_kwargs):
    """Drive Lorenz96's full mode-3 DA cycle. Registered under ``"lorenz"``."""

    icesee_kwargs = normalize_icesee_kwargs(icesee_kwargs)
    normalize_execution_mode(icesee_kwargs, expected=3)

    world = MPI.COMM_WORLD
    world_rank = int(world.Get_rank())

    # start the timer (mirrors icesee_da_serial.py / icesee_da_partial_parallel.py
    # so all four execution modes share one performance-metrics display)
    global_start_time = MPI.Wtime()
    time_forecast_step = 0.0
    time_forecast_file_writing = 0.0

    model = icesee_kwargs.get("model_name")
    Nens = int(icesee_kwargs["Nens"])
    nd = int(icesee_kwargs.get("nd", icesee_kwargs.get("num_state_vars", 3)))
    nt = int(icesee_kwargs["nt"])

    # --- setup phase: mirrors icesee_da_serial.py's model-agnostic setup
    # block (dim_list/observed_vars_params/all_observed/model_module),
    # cheap and deterministic so every rank computes it identically. ---
    model_module = SupportedModels(
        model=model, verbose=icesee_kwargs.get("verbose")
    ).call_model()
    icesee_kwargs.update(
        {
            "dim_list": [nd],
            "model_module": model_module,
            "vec_inputs_old": icesee_kwargs.get("vec_inputs"),
            "model_nprocs": icesee_kwargs.get("model_nprocs") or 1,
        }
    )
    icesee_kwargs["observed_vars_params"] = list(
        icesee_kwargs.get("observed_vars", [])
    ) + list(icesee_kwargs.get("observed_params", []))
    icesee_kwargs["all_observed"] = icesee_kwargs["observed_vars_params"]

    _modelrun_datasets = icesee_kwargs.get("data_path") or "_modelrun_datasets"
    os.environ["ICESEE_RESULTS_DIR"] = str(_modelrun_datasets)
    if world_rank == 0 and not os.path.exists(_modelrun_datasets):
        os.makedirs(_modelrun_datasets, exist_ok=True)
    world.Barrier()

    _true_nurged = f"{_modelrun_datasets}/true_nurged_states.h5"
    _synthetic_obs = f"{_modelrun_datasets}/synthetic_obs.h5"
    _ensemble_data = f"{_modelrun_datasets}/icesee_ensemble_data.h5"
    icesee_kwargs.update(
        {"true_nurged_file": _true_nurged, "synthetic_obs_file": _synthetic_obs}
    )

    # --- file-writing setup calls: real world root only (see module
    # docstring -- these three functions hardcode a single-rank identity
    # internally and would race if called from every rank). Non-root ranks
    # keep these timings at 0.0; the MPI.MAX reduction below picks up
    # root's real value, matching modes 0/1/2's timing-aggregation pattern.
    true_wrong_time = 0.0
    observation_time = 0.0
    time_init_ensemble_mean = 0.0
    init_file_time = 0.0
    if world_rank == 0:
        _t = MPI.Wtime()
        icesee_kwargs = generate_true_wrong_state(**icesee_kwargs)
        true_wrong_time = MPI.Wtime() - _t
        _t = MPI.Wtime()
        icesee_kwargs = generate_synthetic_observations(**icesee_kwargs)
        observation_time = MPI.Wtime() - _t

        _t = MPI.Wtime()
        (
            icesee_kwargs,
            _ens_vec,
            _time_init_noise_generation,
            time_init_ensemble_mean,
            init_file_time,
            _shape_ens,
            _ensemble_bg,
            _ensemble_vec_mean,
            _ensemble_vec_full,
        ) = ensemble_initialization(**icesee_kwargs)
        ensemble_init_time = MPI.Wtime() - _t
    else:
        ensemble_init_time = 0.0
    world.Barrier()

    with h5py.File(_ensemble_data, "r") as f:
        ensemble_vec = f["ensemble"][:, :, 0]
    with h5py.File(_synthetic_obs, "r") as f:
        hu_obs = f["hu_obs"][:]
        error_R = f["R"][:]

    icesee_kwargs["error_R"] = error_R
    icesee_kwargs["native_initial_ensemble"] = {
        i: np.asarray(ensemble_vec[:, i], dtype=np.float64) for i in range(Nens)
    }

    topology = create_distributed_topology(world, spatial_ranks=1)
    adapter = LorenzNativeAdapter()
    pool = initialize_native_member_pool(adapter, topology, icesee_kwargs)

    obs_indices = np.asarray(
        UtilsFunctions(icesee_kwargs).JObs_indices(nd), dtype=np.int64
    )
    error_mode = str(
        icesee_kwargs.get("enkf_observation_error_mode", "legacy_prior_anomalies")
    ).lower()
    if error_mode not in _SUPPORTED_ERROR_MODES:
        raise NotImplementedError(
            "execution_mode 3's Lorenz96 runner supports "
            "'legacy_prior_anomalies' and 'stochastic_R' observation-error "
            f"modes; got {error_mode!r}. See the module docstring for why "
            "'generated_R' is not (yet) supported here."
        )
    base_seed = int(icesee_kwargs.get("base_seed", 42))

    obs_t, ind_m, m_obs = UtilsFunctions(icesee_kwargs).generate_observation_schedule(
        **icesee_kwargs
    )
    # Modes 0-2 all re-publish this freshly computed schedule back onto
    # icesee_kwargs['obs_index']/'number_obs_instants' before their DA loop
    # and final file write (icesee_da_serial.py, synchronize_observation_schedule
    # for modes 1/2) because the value config-loading computes at import time
    # can go stale (it is derived from an early, transient 't'/'nt' pair that
    # gets corrected later by the application entry point). Mirror that here
    # so the 'obs_index' written to the true-wrong-<model>.h5 metadata file
    # (used by read_results.ipynb) matches hu_obs's actual column count
    # instead of the stale config-time value.
    icesee_kwargs.update(
        {"obs_t": obs_t, "obs_index": ind_m, "number_obs_instants": m_obs, "m_obs": m_obs}
    )
    km = 0

    # --- initial condition (column 0) is already written by
    # ensemble_initialization; write it again from the native pool's initial
    # snapshot is unnecessary since it is bit-identical by construction. ---

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
        _write_ensemble_column(
            topology, result.local_analysis.members, k + 1, Nens, nd, _ensemble_data
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
        phases={"true_wrong_state": true_wrong_time, "observation_generation": observation_time, "ensemble_init": ensemble_init_time, "forecast_step": time_forecast_step, "init_file_io": init_file_time, "forecast_file_io": time_forecast_file_writing, "analysis_file_io": analysis_file_time, "init_ensemble_mean": time_init_ensemble_mean, "analysis_step": None if km else 0.0},
        counts={"forecast_step": nt, "analysis_step": km},
        output_dir=_modelrun_datasets,
    )

    return icesee_kwargs


register_execution_mode_3(
    "lorenz",
    run_lorenz96_execution_mode_3,
    notes=(
        "Production mode-3 DA-cycle runner using LorenzNativeAdapter "
        "(distributed_native_adapter.py) and "
        "src/parallelization/distributed_native_cycle.py. Every native "
        "member owns its full 3-element state locally (spatial_ranks=1); "
        "mode 3's actual contract here is bounded per-rank ensemble memory, "
        "not spatial decomposition -- see distributed_native_adapter.py's "
        "module docstring."
    ),
)
