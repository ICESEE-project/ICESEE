# ==============================================================================
# @des: Production execution_mode 3 DA-cycle runner for the idealized_pig
#       Icepack application. Registers itself with
#       ``src/parallelization/distributed_mode3_registry.py`` at import time.
# @date: 2026-08-25
# @author: Brian Kyanjo
# ==============================================================================
"""Native distributed assimilation for Idealized PIG.

World rank zero generates truth and observations on COMM_SELF. Each model
group initializes persistent Firedrake fields on its spatial communicator;
initial thickness is unperturbed. Observation rows are mapped from serial
geometry to distributed ownership before analysis.

The supported observation-error modes are legacy_prior_anomalies and
stochastic_R. Full hu_obs and R arrays are currently loaded on every rank.

Restart state is written as owned-row HDF5 shards under
<data_path>/_mode3_state_history/, with a manifest committed after all shards.
Initial state is in initial_condition/checkpoint_00000000; step k checkpoints
contain state after cycle k (history column k+1). Retention and resume use
Mode3CheckpointManager. Full-history retention can additionally materialize
icesee_ensemble_data.h5 in serial reference numbering after the run. Flowline
diagnostics are separate HDF5 files."""

from __future__ import annotations

import os
import gc

import h5py
import numpy as np
from mpi4py import MPI
from firedrake import COMM_SELF as _REFERENCE_COMM

from ICESEE.applications.icepack_model.examples.idealized_pig._icepack_model import (
    initialize_model,
)
from ICESEE.applications.icepack_model.examples.idealized_pig._icepack_native import (
    IDEALIZED_PIG_NATIVE_ADAPTER,
    node_coordinates,
    owned_node_layout,
    release_model_context,
)
from ICESEE.applications.icepack_model.icepack_utils._reference_numbering import (
    match_owned_nodes_to_reference,
    reference_to_layout_permutation,
    translate_observation_ids,
)
from ICESEE.applications.supported_models import SupportedModels
from ICESEE.src.EnKF._generate_synthetic_observations import (
    generate_synthetic_observations,
)
from ICESEE.src.EnKF._generate_true_wrong_state import generate_true_wrong_state
from ICESEE.src.parallelization._mpi_analysis_functions import apply_analysis_controls_local
from ICESEE.src.parallelization.native_process_noise import (
    NativeProcessNoiseAdapter, canonical_noise_identity,
)
from ICESEE.src.parallelization.mode3_checkpointing import (
    HISTORY_DIRNAME,
    CheckpointPolicy,
    Mode3CheckpointManager,
    fresh_start_plan,
)
from ICESEE.src.parallelization.distributed_history_export import (
    finalize_ensemble_history,
    prepare_reference_layout,
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
    native_layout_fingerprint,
    restore_native_member_pool,
)
from ICESEE.src.parallelization.distributed_streaming_runtime import (
    run_native_store_streaming_analysis_cycle,
)
from ICESEE.src.parallelization.parallel_mpi.model_capabilities import (
    validate_node_local_model_groups,
)
from ICESEE.src.parallelization.distributed_topology import (
    create_distributed_topology,
    register_topology_run_metadata,
)
from ICESEE.src.utils.icesee_context import (
    normalize_execution_mode,
    normalize_icesee_kwargs,
)
from ICESEE.src.utils.init_trace import init_trace_enabled, trace_init, traced_init
from ICESEE.src.utils.localization import active_observation_std
from ICESEE.src.utils.performance import (
    emit_performance_report,
    record_phase,
    register_io_provider,
    register_package_versions,
    register_run_metadata,
)
from ICESEE.src.utils.tools import save_all_data
from ICESEE.src.utils.utils import UtilsFunctions

_SUPPORTED_ERROR_MODES = {"legacy_prior_anomalies", "stochastic_r"}
_DEFAULT_RUN_ID = "icepack-idealized-pig-mode3"


def _release_model_topology(topology):
    """Collectively release owned splits before PyOP2's interpreter teardown.

    PyOP2 attaches an inner communicator to the user spatial communicator.
    Free the owning outer communicator while those attributes and inner
    communicators are still valid; never free COMM_WORLD here.
    """
    release_model_context(topology)
    gc.collect()
    seen = set()
    for name in ("ensemble_comm", "spatial_comm"):
        comm = getattr(topology, name, None)
        if comm is None or comm is getattr(topology, "world", None) or id(comm) in seen:
            continue
        seen.add(id(comm))
        comm.Free()


def _build_reference_setup_context(icesee_kwargs: dict) -> int:
    """Rebuild a correct, single-rank reference context and return ``nd``.

    Only ever called from the real MPI world root (see module docstring).
    Overrides every Firedrake-object kwarg ``generate_true_wrong_state``/
    ``generate_synthetic_observations`` need, using ``comm=MPI.COMM_SELF``
    so the resulting node counts are global regardless of world size.
    """

    ref_kwargs = dict(icesee_kwargs)
    ref_kwargs["comm"] = _REFERENCE_COMM
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


def _setup_artifacts_reusable(icesee_kwargs: dict, nd: int) -> bool:
    """Whether a resumed run can reuse the original truth/observation files.

    Reusing them (instead of regenerating the full truth trajectory) keeps
    a resume cheap and assimilates exactly the observations the
    interrupted run used. Shapes must match this run's configuration.
    """

    try:
        with h5py.File(icesee_kwargs["synthetic_obs_file"], "r") as f:
            hu_obs = f["hu_obs"]
            if "R" not in f or hu_obs.shape[0] != nd:
                return False
            if hu_obs.shape[1] != int(icesee_kwargs.get("number_obs_instants", hu_obs.shape[1])):
                return False
        if icesee_kwargs.get("generate_true_state", True):
            with h5py.File(icesee_kwargs["true_nurged_file"], "r") as f:
                if f["true_state"].shape != (nd, int(icesee_kwargs["nt"]) + 1):
                    return False
    except (OSError, KeyError):
        return False
    return True


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
    export_history = bool(icesee_kwargs.get("mode3_export_ensemble_history", False))
    export_rows = int(icesee_kwargs.get("mode3_export_row_chunk_size", 8192))
    checkpoint_policy = CheckpointPolicy.from_kwargs(icesee_kwargs)
    if export_history and (checkpoint_policy.every != 1 or checkpoint_policy.keep_last != 0):
        raise ValueError("Mode-3 plotting history requires checkpoint_every=1 and checkpoint_keep_last=0")
    if export_history and export_rows < 1:
        raise ValueError("mode3_export_row_chunk_size must be positive")

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
    trace_init("mode3 runner entered", Nens=Nens, nt=nt,
               model_nprocs=icesee_kwargs.get("model_nprocs"))
    with traced_init("mode3 barrier (data_path ready)"):
        world.Barrier()

    _true_nurged = f"{_modelrun_datasets}/true_nurged_states.h5"
    _synthetic_obs = f"{_modelrun_datasets}/synthetic_obs.h5"
    icesee_kwargs.update(
        {"true_nurged_file": _true_nurged, "synthetic_obs_file": _synthetic_obs}
    )

    # Preserve model_nprocs across truth and observation setup.
    _model_nprocs_before_true_wrong_state = icesee_kwargs.get("model_nprocs")

    resume = bool(icesee_kwargs.get("resume_from_checkpoint", False))
    true_wrong_time = 0.0
    observation_time = 0.0
    nd = None
    reference_coords = None
    if world_rank == 0:
        # The single-rank reference model the truth is generated with is set
        # up here; it is reported as its own phase so that True/Wrong State
        # Time covers trajectory generation only, as in modes 0-2.
        _t = MPI.Wtime()
        with traced_init("truth model setup (COMM_SELF reference)"):
            nd = _build_reference_setup_context(icesee_kwargs)
            # Truth and observation rows use this serial numbering.
            reference_coords = node_coordinates(icesee_kwargs["Q"])
        record_phase("truth_model_setup", MPI.Wtime() - _t)

        if resume and _setup_artifacts_reusable(icesee_kwargs, nd):
            print(
                "[ICESEE] resume: reusing existing "
                f"{_true_nurged} and {_synthetic_obs}",
                flush=True,
            )
            record_phase("truth_generation", 0.0, operations=0)
            record_phase("wrong_reference_generation", 0.0, operations=0)
        else:
            _t = MPI.Wtime()
            with traced_init("truth generation", nd=nd):
                icesee_kwargs = generate_true_wrong_state(**icesee_kwargs)
            true_wrong_time = MPI.Wtime() - _t

            _t = MPI.Wtime()
            with traced_init("observation generation"):
                icesee_kwargs = generate_synthetic_observations(**icesee_kwargs)
            observation_time = MPI.Wtime() - _t

    with traced_init("mode3 barrier (truth/observations ready)"):
        world.Barrier()

    icesee_kwargs["model_nprocs"] = _model_nprocs_before_true_wrong_state

    nd = int(world.bcast(nd, root=0))
    icesee_kwargs["nd"] = nd
    reference_coords = world.bcast(reference_coords, root=0)

    with h5py.File(_synthetic_obs, "r") as f:
        hu_obs = f["hu_obs"][:]
        error_R = f["R"][:]
    icesee_kwargs["error_R"] = error_R

    with traced_init("topology communicator creation",
                     model_nprocs=icesee_kwargs["model_nprocs"]):
        topology = create_distributed_topology(
            world, spatial_ranks=int(icesee_kwargs["model_nprocs"])
        )
    # The per-group model setup below opens the initialization checkpoint
    # collectively on topology.spatial_comm; refuse node-spanning groups
    # before reaching it (see model_capabilities.py).
    _spatial_ranks = int(icesee_kwargs["model_nprocs"])
    validate_node_local_model_groups(
        model, world, lambda r: r // _spatial_ranks, _spatial_ranks
    )
    _group = getattr(topology, "ensemble_slot", None)
    if init_trace_enabled():
        trace_init(
            "topology",
            group=_group,
            groups=getattr(topology, "ensemble_groups", None),
            model_rank=getattr(topology, "spatial_rank", None),
            model_comm=getattr(topology, "spatial_comm", None),
            ensemble_comm=getattr(topology, "ensemble_comm", None),
        )
    adapter = IDEALIZED_PIG_NATIVE_ADAPTER

    obs_t, ind_m, m_obs = UtilsFunctions(icesee_kwargs).generate_observation_schedule(
        **icesee_kwargs
    )
    icesee_kwargs.update(
        {"obs_t": obs_t, "obs_index": ind_m, "number_obs_instants": m_obs, "m_obs": m_obs}
    )

    run_id = str(icesee_kwargs.get("run_id", _DEFAULT_RUN_ID)) + ":native-controls-v1:" + canonical_noise_identity(icesee_kwargs)
    checkpoints = Mode3CheckpointManager(
        os.path.join(_modelrun_datasets, HISTORY_DIRNAME),
        run_id=run_id,
        topology=topology,
        number_of_members=Nens,
        nt=nt,
        analysis_steps=[int(step) for step in np.asarray(ind_m)[:m_obs]],
        policy=checkpoint_policy,
        time_grid=icesee_kwargs.get("t"),
        layout_fingerprint=native_layout_fingerprint(adapter, topology, icesee_kwargs),
    )
    plan = checkpoints.prepare_resume() if resume else fresh_start_plan()
    checkpoints.report_plan(plan)

    ensemble_init_time = 0.0
    restart_restore_time = 0.0
    _t = MPI.Wtime()

    if plan.checkpoint_path is not None:
        # Resume from the last valid committed checkpoint.
        # Do not regenerate the initial ensemble.
        with traced_init("ensemble restoration (native member pool)", group=_group):
            pool, _ = restore_native_member_pool(
                adapter,
                topology,
                icesee_kwargs,
                plan.checkpoint_path,
                expected_run_id=run_id,
            )

        restart_restore_time = MPI.Wtime() - _t

    else:
        with traced_init("ensemble initialization (native member pool)", group=_group):
            pool = initialize_native_member_pool(
                adapter, topology, icesee_kwargs
            )

        if init_trace_enabled():
            trace_init(
                "members assigned",
                group=_group,
                members=list(getattr(pool, "member_ids", [])),
            )

        ensemble_init_time = MPI.Wtime() - _t

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
    # Observation ids are rows of the serial reference numbering (truth and
    # synthetic observations are generated on COMM_SELF). A model group with
    # several spatial ranks renumbers the mesh, so translate each id to the
    # distributed row of the same physical node before routing; order is
    # kept, so `obs_positions` still index the canonical arrays.
    owned_start, global_nodes, owned_coords = owned_node_layout(
        topology=topology, icesee_kwargs=icesee_kwargs
    )
    if global_nodes * int(icesee_kwargs["total_state_param_vars"]) != nd:
        raise RuntimeError(
            f"distributed layout has {global_nodes} nodes per variable but the "
            f"reference state has nd={nd}"
        )
    owned_reference_nodes = match_owned_nodes_to_reference(reference_coords, owned_coords)
    reference_to_layout = reference_to_layout_permutation(
        owned_start,
        owned_reference_nodes,
        topology.spatial_comm,
        global_nodes,
    )
    if export_history:
        prepare_reference_layout(
            world, checkpoints.history_root, run_id=run_id,
            layout_id=pool.layout.layout_id,
            row_maps={block.name: reference_to_layout for block in pool.layout.blocks},
        )
    # Apply the same configured forecast noise as modes 1/2, at the native
    # model boundary. Initial thickness remains unperturbed. Noise is keyed
    # in serial reference numbering, then scattered to physical owned rows.
    canonical_rows = np.concatenate([
        block.global_offset + owned_reference_nodes for block in pool.layout.blocks
    ])
    adapter = NativeProcessNoiseAdapter(
        adapter, topology=topology, canonical_rows=canonical_rows,
        reference_coords=reference_coords, icesee_kwargs=icesee_kwargs,
        workspace_budget_mb=icesee_kwargs.get("native_noise_workspace_mb", 256),
    )
    if world_rank == 0:
        print("[ICESEE] Native forecast noise: canonical configured schedule; "
              "initial fields unchanged; bounded leader compatibility workspace",flush=True)
    del reference_coords
    obs_positions, obs_local_indices = pool.route_observation_ids(
        translate_observation_ids(obs_indices, reference_to_layout),
        adapter,
        topology=topology,
        icesee_kwargs=icesee_kwargs,
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
    # Analyses already applied before the restored checkpoint (0 when
    # starting fresh); the loop continues the same schedule from there.
    km = plan.analyses_completed

    # --- initial (pre-forecast) native ensemble, written once before the
    # timestep loop -- the mode-3 equivalent of modes 0-2's ensemble column
    # 0 (see module docstring for the checkpoint directory layout). A
    # resumed run keeps the one its original run wrote. ---
    # Checkpoint write counters for the performance summary; bytes are the
    # owned state arrays handed to the checkpoint writer.
    checkpoint_io = {"bytes_written": 0.0, "write_time_s": 0.0, "writes": 0}
    register_io_provider("checkpoint", lambda: checkpoint_io)

    def _record_checkpoint(snapshot, seconds):
        members = getattr(snapshot, "members", None)
        if isinstance(members, dict):
            checkpoint_io["bytes_written"] += sum(
                getattr(array, "nbytes", 0) for array in members.values()
            )
        checkpoint_io["write_time_s"] += seconds
        checkpoint_io["writes"] += 1

    checkpoints.storage_preflight(
        start_step=plan.start_step,
        bytes_per_checkpoint=int(pool.layout.global_size) * Nens * np.dtype(np.float64).itemsize,
        files_per_checkpoint=int(world.Get_size()) + 1,
        include_initial=True,
    )

    init_file_time = 0.0
    if not plan.resumed:
        _initial_snapshot = pool.snapshot_owned()
        _t = MPI.Wtime()
        checkpoints.write_initial(
            _initial_snapshot, extra={
                "description": "initial (pre-forecast) ensemble",
                "observation_columns": np.asarray(ind_m, dtype=int).tolist(),
                "observation_max_time": float(icesee_kwargs["obs_max_time"]),
                "execution_flag": int(icesee_kwargs["execution_flag"]),
            }
        )
        init_file_time = MPI.Wtime() - _t
        _record_checkpoint(_initial_snapshot, init_file_time)
        del _initial_snapshot
    time_analysis_cycle_steps = 0.0
    analyses_this_run = 0

    for k in range(plan.start_step, nt):
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

        # Report individual step wall times when verbose output is enabled.
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
            analysis_row_controls=apply_analysis_controls_local,
        )
        _step_wall = MPI.Wtime() - _t
        time_forecast_step += _step_wall
        if do_analysis:
            time_analysis_cycle_steps += _step_wall
        # Report individual step wall times when verbose output is enabled.
        if icesee_kwargs.get("verbose") and world_rank == 0:
            print(f"[ICESEE] mode3 forecast step {k}: wall={_step_wall:.4f}s did_analysis={do_analysis}", flush=True)
        if do_analysis:
            km += 1
            analyses_this_run += 1

        # Restart checkpoint of the state after cycle k. idealized_pig has no
        # state-modifying analysis finalizer, so result.local_analysis is
        # exactly the pool's live state. Retention prunes older rolling
        # checkpoints only after this one has committed.
        _t = MPI.Wtime()
        _written = checkpoints.commit_step(
            k,
            result.local_analysis,
            did_analysis=do_analysis,
            analyses_completed=km,
            extra={"observation_rows": result.observation_rows, "error_mode": error_mode},
        )
        if _written is not None:
            _checkpoint_s = MPI.Wtime() - _t
            time_forecast_file_writing += _checkpoint_s
            _record_checkpoint(result.local_analysis, _checkpoint_s)

    world.Barrier()

    history_export_time = 0.0
    if export_history:
        _t = MPI.Wtime()
        with traced_init("native ensemble history finalization"):
            finalize_ensemble_history(world, _modelrun_datasets, row_chunk_size=export_rows)
        if world_rank == 0:
            history_export_time = MPI.Wtime() - _t

    # Inactive-member store I/O stats (Gate 11) and cleanup (Gate 3):
    # every rank owns and cleans up only its OWN store file -- no
    # additional cross-rank coordination needed beyond the barrier above,
    # since ownership never changes and no rank ever opens another rank's
    # file. Durable checkpoints (the checkpoint manager's roots) are a
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

    analysis_file_time = history_export_time
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
        analysis_file_time += MPI.Wtime() - _t

    # ─────────────────────────────────────────────────────────────
    #  End Timer and Aggregate Elapsed Time Across Ranks
    # ─────────────────────────────────────────────────────────────
    global_elapsed_time = MPI.Wtime() - global_start_time

    steps_run = nt - plan.start_step
    register_run_metadata(
        execution_mode=icesee_kwargs.get("execution_mode"),
        model=icesee_kwargs.get("model_name"),
        forecast_steps=steps_run,
        analysis_events=analyses_this_run,
        checkpoints=checkpoint_io["writes"],
        checkpoint_policy=checkpoints.policy.describe(),
        checkpoints_pruned=len(checkpoints.pruned),
        ensemble_history_export=export_history,
        resumed_from_step=plan.completed_step if plan.resumed else None,
        member_store_backend=icesee_kwargs.get("member_store_backend", "memory"),
        store_streaming_analysis=bool(icesee_kwargs.get("use_store_streaming_analysis", False)),
    )
    adapter.close()
    register_topology_run_metadata(topology, Nens)
    register_package_versions("petsc4py", "firedrake", "icepack")
    # The native cycle fuses each analysis into its forecast step, so there
    # is no separate analysis timer: "analysis_step" is not measured when
    # analyses ran (and has no events when none did), and the steps that
    # carried an analysis are timed as "forecast_step/with_analysis", a
    # subset of "forecast_step". There
    # is no separate ensemble-mean phase. Forecast timing includes process
    # noise and deterministic AR replay when restoring a checkpoint.
    emit_performance_report(
        world,
        elapsed_s=global_elapsed_time,
        phases={
            "true_wrong_state": true_wrong_time,
            "observation_generation": observation_time,
            "ensemble_init": ensemble_init_time,
            "restart_restore": restart_restore_time,
            "forecast_step": time_forecast_step,
            "analysis_step": None if analyses_this_run else 0.0,
            "forecast_step/with_analysis": time_analysis_cycle_steps,
            "init_file_io": init_file_time,
            "forecast_file_io": time_forecast_file_writing,
            "analysis_file_io": analysis_file_time,
            "analysis_file_io/ensemble_history_export": history_export_time,
        },
        counts={
            "forecast_step": steps_run,
            "analysis_step": analyses_this_run,
            "forecast_step/with_analysis": analyses_this_run,
            # A resumed run restores the ensemble instead of initializing it.
            "ensemble_init": 0 if plan.resumed else 1,
            "restart_restore": 1 if plan.resumed else 0,
            "analysis_file_io/ensemble_history_export": int(export_history and world_rank == 0),
        },
        output_dir=_modelrun_datasets,
    )

    pool = None
    _release_model_topology(topology)
    return icesee_kwargs


register_execution_mode_3(
    "icepack",
    run_icepack_execution_mode_3,
    notes=(
        "Production mode-3 DA-cycle runner for idealized_pig using "
        "IDEALIZED_PIG_NATIVE_ADAPTER (_icepack_native.py) and "
        "src/parallelization/distributed_native_cycle.py. Per-timestep "
        "output uses rank-sharded distributed checkpoints. With full retention, "
        "bounded finalization publishes icesee_ensemble_data.h5 in serial "
        "reference row order for plotting."
    ),
)
