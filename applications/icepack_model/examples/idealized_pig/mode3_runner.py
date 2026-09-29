# ==============================================================================
# @des: Production execution_mode 3 DA-cycle runner for the idealized_pig
#       Icepack application. Registers itself with
#       ``src/parallelization/distributed_mode3_registry.py`` at import time.
# @date: 2026-08-25
# @author: Brian Kyanjo
# ==============================================================================
"""Production ``execution_mode: 3`` DA-cycle runner for idealized_pig.

This module wires ``IDEALIZED_PIG_NATIVE_ADAPTER``
(``applications/icepack_model/examples/idealized_pig/_icepack_native.py``)
into ``src/parallelization/distributed_native_cycle.py`` to make mode 3
actually functional for idealized_pig, and registers itself with
``src/parallelization/distributed_mode3_registry.py`` so
``icesee_da_distributed.py`` can dispatch to it instead of raising
``NotImplementedError``. Structure mirrors
``applications/lorenz_model/lorenz_utils/mode3_runner.py`` closely; the
differences below are specific to idealized_pig's Firedrake/PETSc state.

Setup-phase reuse -- and one deliberate divergence from Lorenz96
------------------------------------------------------------------
Generating the true/nurged trajectories and synthetic observations reuses
the exact same model-agnostic, single-rank generators execution mode 0 uses
(``src/EnKF/_generate_true_wrong_state.py``,
``src/EnKF/_generate_synthetic_observations.py``), called only from the real
MPI world root, exactly like Lorenz96's mode-3 runner.

Unlike Lorenz96, this runner does **not** call
``src/EnKF/_ensemble_initialization.py``: idealized_pig's native adapter
(``_icepack_native.initialize_member``) does not depend on that generic
setup-phase HDF5 artifact at all -- it reconstructs each member's initial
perturbed thickness directly and deterministically via
``_member_initialization_context`` (the same member-keyed RNG stream modes
1/2 use), so there is nothing for ``ensemble_initialization`` to produce
that the native adapter would consume.

``run_da_icepack.py``'s unconditional pre-dispatch call to
``initialize_model(comm=MPI.COMM_WORLD)`` genuinely MPI-partitions the
Firedrake mesh across the full world communicator -- correct for modes 0-2
(which forecast on that same communicator) but wrong for mode 3 whenever
``world size > 1``: ``icesee_kwargs['nd']`` as computed there
(``h0.dat.data.size * total_state_param_vars``) would be each rank's
*local* mesh-node count, not the true global count the true/nurged/
observation HDF5 artifacts must be shaped with. Rather than touch
``run_da_icepack.py`` (and risk regressing modes 0-2), this runner
independently rebuilds a correct, single-rank reference context by calling
the exact same ``initialize_model`` function modes 0-2 use, but with
``comm`` overridden to ``MPI.COMM_SELF`` and only from the real world root
-- giving the true global node count regardless of ``world size``, exactly
as a serial (mode 0) run would compute it. The resulting ``nd`` is broadcast
to every rank, overwriting whatever (possibly wrong, for world size > 1)
value ``run_da_icepack.py``'s own pre-dispatch call left in
``icesee_kwargs``. This reference mesh/solver/fields are local variables of
the setup phase only; they are never handed to the native adapter, which
builds its own per-ensemble-group Firedrake objects on
``topology.spatial_comm`` (see ``_icepack_native._build_shared_context``).

Observation-error parity
-------------------------
idealized_pig's default ``enkf_observation_error_mode`` (from
``use_ensemble_pertubations: true`` in ``params.yaml``) is
``"legacy_prior_anomalies"``, exactly matching modes 1/2. ``"stochastic_R"``
is also supported, exactly mirroring Lorenz96's runner (see that module's
docstring for why both are rank/order-independent). ``"generated_R"`` is not
supported here; selecting it raises ``NotImplementedError``.

Output writing -- the second deliberate divergence from Lorenz96
-------------------------------------------------------------------
Lorenz96's ``nd=3`` state is trivially small, so its runner gathers each
timestep's analyzed ensemble to the world root and writes one column into
the same ``icesee_ensemble_data.h5`` schema modes 0-2 produce. idealized_pig's
state (mesh nodes x 5 variables) does not fit that "trivially small
regardless" justification, and gathering a complete member to one rank on
every step would reintroduce the unbounded-per-rank-memory pattern mode 3
exists to avoid. This runner instead uses the mode-3-native, rank-sharded,
no-gather checkpoint primitive (``src/parallelization/distributed_checkpoint.py``,
``save_distributed_checkpoint``/``load_distributed_checkpoint``): every rank
writes only its own locally owned interval of every scheduled member to its
own HDF5 shard, and a manifest is published only once every shard for that
step is durable. This is an honest schema divergence from modes 0-2's single
``icesee_ensemble_data.h5`` file, not a compatible drop-in replacement --
downstream plotting/post-processing built against that schema will not read
mode-3 output directly. Directory layout under
``<data_path>/_mode3_state_history/``:

- ``initial_condition/step_00000000/`` -- the native pool's initial
  (pre-forecast) ensemble, written once before the timestep loop.
- ``steps/step_<k:08d>/`` for ``k`` in ``[0, nt)`` -- the analyzed (or, on
  non-observed steps, forecast-only) ensemble *after* native cycle ``k``,
  i.e. the mode-3 equivalent of modes 0-2's ensemble column ``k + 1``.

Scope caveat: every rank still reads back the *full* ``hu_obs``/``R``
synthetic-observation arrays in memory after the setup phase (matching
Lorenz96's pattern and modes 0-2's existing behavior -- not a new
regression). This is a conscious scope-limiting choice for this first
working version, justified because idealized_pig's default
``model_nprocs``/``spatial_ranks`` is 1 (no actual spatial mesh partitioning
is exercised yet), making "this rank's spatially local observations" equal
to the whole observation set regardless. Revisit if ``model_nprocs > 1`` is
exercised for this application.
"""

from __future__ import annotations

import os

import h5py
import numpy as np
from mpi4py import MPI

from ICESEE.applications.icepack_model.examples.idealized_pig._icepack_model import (
    initialize_model,
)
from ICESEE.applications.icepack_model.examples.idealized_pig._icepack_native import (
    IDEALIZED_PIG_NATIVE_ADAPTER,
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
from ICESEE.src.parallelization.distributed_streaming_runtime import (
    run_native_store_streaming_analysis_cycle,
)
from ICESEE.src.parallelization.distributed_topology import (
    create_distributed_topology,
)
from ICESEE.src.utils.icesee_context import (
    normalize_execution_mode,
    normalize_icesee_kwargs,
)
from ICESEE.src.utils.localization import active_observation_std
from ICESEE.src.utils.tools import save_all_data, display_timing_verbose
from ICESEE.src.utils.utils import UtilsFunctions

_SUPPORTED_ERROR_MODES = {"legacy_prior_anomalies", "stochastic_r"}
_DEFAULT_RUN_ID = "icepack-idealized-pig-mode3"


def _build_reference_setup_context(icesee_kwargs: dict) -> int:
    """Rebuild a correct, single-rank reference context and return ``nd``.

    Only ever called from the real MPI world root (see module docstring).
    Overrides every Firedrake-object kwarg ``generate_true_wrong_state``/
    ``generate_synthetic_observations`` need, using ``comm=MPI.COMM_SELF``
    so the resulting node counts are global regardless of world size.
    """

    ref_kwargs = dict(icesee_kwargs)
    ref_kwargs["comm"] = MPI.COMM_SELF
    (
        h,
        h0,
        s,
        s0,
        u,
        bed,
        zF,
        grounded,
        floating,
        A0,
        beta0,
        smb,
        basal_melt_field,
        Q,
        V,
        forward_solver,
    ) = initialize_model(**ref_kwargs)

    nd = int(h0.dat.data.size * icesee_kwargs["total_state_param_vars"])
    icesee_kwargs.update(
        {
            "h": h,
            "h0": h0,
            "s": s,
            "s0": s0,
            "u": u,
            "bed": bed,
            "zF": zF,
            "grounded": grounded,
            "floating": floating,
            "A0": A0,
            "beta0": beta0,
            "smb": smb,
            "basal_melt_field": basal_melt_field,
            "Q": Q,
            "V": V,
            "solver": forward_solver,
            "nd": nd,
        }
    )
    return nd


def run_icepack_execution_mode_3(**icesee_kwargs):
    """Drive idealized_pig's full mode-3 DA cycle. Registered under ``"icepack"``."""

    icesee_kwargs = normalize_icesee_kwargs(icesee_kwargs)
    normalize_execution_mode(icesee_kwargs, expected=3)

    world = MPI.COMM_WORLD
    world_rank = int(world.Get_rank())

    global_start_time = MPI.Wtime()
    time_forecast_step = 0.0
    time_forecast_file_writing = 0.0

    model = icesee_kwargs.get("model_name")
    Nens = int(icesee_kwargs["Nens"])
    nt = int(icesee_kwargs["nt"])

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

    # --- setup phase: real world root only (see module docstring -- both
    # the reference-context rebuild and the two file-writing generators
    # hardcode/require a single-rank identity and would race or compute a
    # wrong global nd if run from every rank). Non-root ranks keep these
    # timings at 0.0; the MPI.MAX reduction below picks up root's real
    # value, matching modes 0/1/2's timing-aggregation pattern. ---
    # generate_true_wrong_state/generate_synthetic_observations are generic
    # (modes-0-2-shared) functions that rebuild their returned kwargs dict
    # and do not preserve every key they don't recognize -- confirmed
    # directly (2026-09-28 PACE calibration debugging) that they silently
    # drop `model_nprocs` on world_rank==0's own copy of icesee_kwargs
    # only (every other rank's copy is untouched, since this block is
    # rank-0-only). Left uncorrected, the topology built below from
    # icesee_kwargs["model_nprocs"] would then be built from an
    # inconsistent value on rank 0 vs every other rank. Rather than
    # changing the shared generic function (out of scope -- it is used by
    # modes 0-2 too, not Mode-3-specific), re-assert the pre-block value
    # uniformly on every rank once the rank-0-only block has run.
    #
    # GOTCHA while verifying this fix: a console line "[ICESEE] ICESEE
    # topology: world=..., ranks/model=..." appears early in every run's
    # output and looks like it should reflect this -- it does NOT. That
    # line is src/parallelization/parallel_mpi/resource_plan.py's own
    # ResourcePlan display (generic modes-0-2 `ranks_per_model`), a
    # completely decoupled concept from Mode-3's own `model_nprocs`/
    # `topology.spatial_ranks` below (see this module's own docstring).
    # The authoritative value for Mode-3 is `topology.spatial_ranks`
    # itself, or equivalently a checkpoint's own
    # manifest.json["source_topology"]["spatial_ranks"].
    _model_nprocs_before_true_wrong_state = icesee_kwargs.get("model_nprocs")

    true_wrong_time = 0.0
    nd = None
    if world_rank == 0:
        _t = MPI.Wtime()
        nd = _build_reference_setup_context(icesee_kwargs)
        icesee_kwargs = generate_true_wrong_state(**icesee_kwargs)
        icesee_kwargs = generate_synthetic_observations(**icesee_kwargs)
        true_wrong_time = MPI.Wtime() - _t
    world.Barrier()

    icesee_kwargs["model_nprocs"] = _model_nprocs_before_true_wrong_state

    nd = int(world.bcast(nd, root=0))
    icesee_kwargs["nd"] = nd

    with h5py.File(_synthetic_obs, "r") as f:
        hu_obs = f["hu_obs"][:]
        error_R = f["R"][:]
    icesee_kwargs["error_R"] = error_R

    topology = create_distributed_topology(
        world, spatial_ranks=int(icesee_kwargs["model_nprocs"])
    )
    adapter = IDEALIZED_PIG_NATIVE_ADAPTER
    pool = initialize_native_member_pool(adapter, topology, icesee_kwargs)

    obs_indices = np.asarray(
        UtilsFunctions(icesee_kwargs).JObs_indices(nd), dtype=np.int64
    )
    # Route the canonical (world-identical) observation-index array to this
    # spatial rank's owned mesh-node subset ONCE, since Firedrake's mesh
    # partition is fixed for the whole run (see _icepack_native._shared_
    # context's per-ensemble-group cache). `obs_positions` are indices back
    # into `obs_indices` (and therefore into any Nobs-sized array keyed by
    # the same canonical ordering, e.g. `d`/`member_errors` below); every
    # analysis step subsets by `obs_positions` before building its
    # NativeObservationBatch, whose own `observation_ids` must already be
    # spatially local -- see run_native_global_analysis_cycle's docstring.
    # This is the fix for the P>1 spatial-decomposition
    # "observation rows are not owned by this spatial rank" ValueError:
    # previously every rank was hand the full, unpartitioned `obs_indices`.
    obs_positions, obs_local_indices = pool.route_observation_ids(
        obs_indices, adapter, topology=topology, icesee_kwargs=icesee_kwargs
    )
    error_mode = str(
        icesee_kwargs.get("enkf_observation_error_mode", "legacy_prior_anomalies")
    ).lower()
    if error_mode not in _SUPPORTED_ERROR_MODES:
        raise NotImplementedError(
            "execution_mode 3's idealized_pig runner supports "
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
            # `d`/`sigma`/`eta` are computed over the FULL, world-identical
            # `obs_indices` (as before the spatial-decomposition fix) so the
            # deterministic RNG stream and scientific values are unchanged
            # for any spatial-rank count; only the final batch handed to
            # run_native_global_analysis_cycle is subset to this rank's
            # owned observations (`obs_positions`/`obs_local_indices`,
            # routed once above).
            y_full = np.asarray(hu_obs[:, km], dtype=np.float64).copy()
            y_full[np.isnan(y_full)] = 0.0
            d = y_full[obs_indices][obs_positions]
            member_errors = None
            if error_mode == "stochastic_r":
                sigma = active_observation_std(icesee_kwargs, km, obs_indices)
                obs_seed = base_seed + 1000003 * (km + 1)
                rng = np.random.default_rng(obs_seed)
                eta = rng.standard_normal((obs_indices.size, Nens)) * sigma[:, None]
                eta -= np.mean(eta, axis=1, keepdims=True)
                member_errors = {
                    member_id: eta[obs_positions, member_id]
                    for member_id in pool.member_ids
                }
            batches.append(
                NativeObservationBatch(
                    observation_ids=obs_local_indices,
                    values=d,
                    member_errors=member_errors,
                )
            )

        # run_native_global_analysis_cycle fuses the forecast and (when
        # scheduled) analysis update into one call at this API boundary, so
        # unlike modes 0-2 the two cannot be timed separately here; the
        # combined per-step cost is reported as forecast_step_time.
        #
        # use_store_streaming_analysis (opt-in, requires use_member_streaming
        # too): selects run_native_store_streaming_analysis_cycle, which
        # pulls/pushes state-row blocks directly through the pool's
        # InactiveMemberStore (transform_members_via_store) instead of
        # requiring this rank's round-assigned members' full owned arrays
        # simultaneously resident (transform_local_members). Verified
        # bit-for-bit equivalent to the original path -- see
        # test_distributed_streaming_runtime.py and this session's real-
        # Icepack equivalence run.
        _t = MPI.Wtime()
        cycle_fn = (
            run_native_store_streaming_analysis_cycle
            if icesee_kwargs.get("use_store_streaming_analysis")
            else run_native_global_analysis_cycle
        )
        result = cycle_fn(
            pool,
            adapter,
            k,
            batches,
            number_of_batches=(1 if do_analysis else 0),
            topology=topology,
            icesee_kwargs=icesee_kwargs,
            error_mode=error_mode,
        )
        _step_wall = MPI.Wtime() - _t
        time_forecast_step += _step_wall
        # PACE calibration instrumentation (2026-09-28, additive only --
        # does not affect control flow or any existing metric). JIT/warm-up
        # cost is concentrated in the first forecast step; printing each
        # step's own wall time lets a post-run parser separate "first step"
        # from "steady-state per-step" cost instead of only ever seeing the
        # aggregated forecast_step_time average.
        if icesee_kwargs.get("verbose") and world_rank == 0:
            print(f"[ICESEE] mode3 forecast step {k}: wall={_step_wall:.4f}s did_analysis={do_analysis}", flush=True)
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

    # Inactive-member store I/O stats (Gate 11) and cleanup (Gate 3):
    # every rank owns and cleans up only its OWN store file -- no
    # additional cross-rank coordination needed beyond the barrier above,
    # since ownership never changes and no rank ever opens another rank's
    # file. Durable checkpoints (steps_root/initial_root above) are a
    # completely separate, untouched mechanism -- this only ever removes
    # the TEMPORARY working store, never a checkpoint.
    _store = getattr(pool, "_store", None)
    if _store is not None and hasattr(_store, "stats"):
        if icesee_kwargs.get("member_store_backend", "memory") != "memory":
            print(
                f"[ICESEE] mode3 member-store stats (world_rank={world_rank}): "
                f"{_store.stats}",
                flush=True,
            )
        if hasattr(_store, "close"):
            _store.close()
        if hasattr(_store, "file_paths"):
            for _fp in _store.file_paths():
                try:
                    os.remove(_fp)
                except OSError:
                    pass
                # Cosmetic only (Phase 25): try to remove the now-possibly-
                # empty run/parent directories. No new synchronization --
                # whichever rank happens to remove the last file in a
                # shared directory succeeds; every other rank's rmdir on
                # a still-non-empty directory just raises OSError, caught
                # and ignored. Never removes data_path itself.
                _parent = os.path.dirname(_fp)
                for _ in range(2):
                    try:
                        os.rmdir(_parent)
                    except OSError:
                        break
                    _parent = os.path.dirname(_parent)

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
    total_elapsed_time = world.allreduce(global_elapsed_time, op=MPI.SUM)
    total_wall_time = world.allreduce(global_elapsed_time, op=MPI.MAX)
    true_wrong_time = world.allreduce(true_wrong_time, op=MPI.MAX)
    forecast_step_time = world.allreduce(time_forecast_step, op=MPI.MAX)
    forecast_file_time = world.allreduce(time_forecast_file_writing, op=MPI.MAX)
    analysis_file_time = world.allreduce(analysis_file_time, op=MPI.MAX)
    init_file_time = world.allreduce(init_file_time, op=MPI.MAX)
    # No ensemble_initialization() call in this runner (see module
    # docstring), so there is no separate ensemble-mean-computation cost to
    # report; the initial-checkpoint write time is attributed to
    # init_file_time above instead.
    ensemble_init_time = 0.0
    analysis_step_time = 0.0  # fused into forecast_step_time; see comment above
    assimilation_time = ensemble_init_time + forecast_step_time + analysis_step_time
    total_file_time = init_file_time + forecast_file_time + analysis_file_time

    world.Barrier()
    if world_rank == 0:
        display_timing_verbose(
            computational_time=total_elapsed_time,
            wallclock_time=total_wall_time,
            true_wrong_time=true_wrong_time,
            assimilation_time=assimilation_time,
            forecast_step_time=forecast_step_time,
            analysis_step_time=analysis_step_time,
            ensemble_init_time=ensemble_init_time,
            init_file_time=init_file_time,
            forecast_file_time=forecast_file_time,
            analysis_file_time=analysis_file_time,
            total_file_time=total_file_time,
            forecast_noise_time=0.0,
            time_init_ensemble_mean_computation=0.0,
            time_forecast_ensemble_mean_computation=0.0,
            time_analysis_ensemble_mean_computation=0.0,
            comm=world,
            # Mode 3's `model_nprocs` (topology.spatial_ranks) is already
            # part of `world`'s own size (no separately spawned children,
            # unlike modes 1/2's coordinator+worker convention this shared
            # display function's "(model_nprocs+1)" multiplier was written
            # for) -- model_nprocs=0 makes that multiplier a no-op so the
            # header reports the true world size instead of double-counting.
            model_nprocs=0,
        )

    return icesee_kwargs


register_execution_mode_3(
    "icepack",
    run_icepack_execution_mode_3,
    notes=(
        "Production mode-3 DA-cycle runner for idealized_pig using "
        "IDEALIZED_PIG_NATIVE_ADAPTER (_icepack_native.py) and "
        "src/parallelization/distributed_native_cycle.py. Per-timestep "
        "output uses rank-sharded distributed checkpoints (no gather), not "
        "modes 0-2's icesee_ensemble_data.h5 schema -- see this module's "
        "docstring for the checkpoint directory layout and other "
        "divergences."
    ),
)
