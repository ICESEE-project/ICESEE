# ==============================================================================
# @des: Production execution_mode 3 DA-cycle runner for the ISSM ISMIP_Choi
#       hybrid application. Registers itself with
#       ``src/parallelization/distributed_mode3_registry.py`` at import time.
# @date: 2026-08-26
# @author: Brian Kyanjo
# ==============================================================================
"""Whole-member Mode-3 execution for the ISMIP_Choi application.

Python ranks reuse MATLAB servers across scheduled members. ISSM controls its
internal FEM parallelism; Python state remains complete-member storage.
Managed continuation commits full MATLAB backend snapshots before their vector
checkpoints and restores them without regenerating ensemble initial conditions.
"""

from __future__ import annotations

from ICESEE.src.utils.timestep_contract import (
    TIME_GRID_CONTRACT, observation_is_due, observation_cycle_steps, completed_observations,
    prepare_timestep_contract,
)
import os
import shutil
import socket

import h5py
import numpy as np
import scipy.io as sio
from mpi4py import MPI

from ICESEE.applications.issm_model.examples.ISMIP_Choi._issm_model import (
    initialize_model,
)
from ICESEE.applications.issm_model.examples.ISMIP_Choi._issm_native import (
    ISMIP_CHOI_NATIVE_ADAPTER,
)
from ICESEE.applications.issm_model.issm_utils.matlab2python.mat2py_utils import (
    MatlabServer,
    add_issm_dir_to_sys_path,
    setup_example_directory,
)
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
    matlab_icesee_kwargs,
    normalize_execution_mode,
    normalize_icesee_kwargs,
)
from ICESEE.src.utils.inference_plugin import resolve_analysis_cycle_time
from ICESEE.src.utils.localization import active_observation_std
from ICESEE.src.utils.tools import save_all_data
from ICESEE.src.utils.performance import emit_performance_report, register_run_metadata
from ICESEE.src.utils.utils import UtilsFunctions

_SUPPORTED_ERROR_MODES = {"legacy_prior_anomalies", "stochastic_r"}
_DEFAULT_RUN_ID = "issm-ismip-choi-mode3"


def _build_static_config(icesee_kwargs: dict, *, icesee_cwd: str,
                          issm_dir: str, issm_examples_dir: str) -> None:
    """Rank-independent config, copied verbatim from ``run_da_issm.py``.

    Every value here is a deterministic transform of already-configured
    ``icesee_kwargs`` entries (no per-member or per-rank identity), so it is
    safe -- and required for byte-identical behavior with modes 0-2 -- for
    every rank to compute it identically.
    """

    icesee_kwargs.update({
        'Lx': int(float(icesee_kwargs.get('Lx'))), 'Ly': int(float(icesee_kwargs.get('Ly'))),
        'nx': int(float(icesee_kwargs.get('nx'))), 'ny': int(float(icesee_kwargs.get('ny'))),
        'ParamFile': icesee_kwargs.get('ParamFile'),
        'cluster_name': socket.gethostname().replace('-', ''),
        'steps': int(float(icesee_kwargs.get('steps'))),
        'dt': float(icesee_kwargs.get('timesteps_per_year')),
        'tinitial': float(icesee_kwargs.get('tinitial')),
        'tfinal': float(icesee_kwargs.get('num_years')),
        't': np.linspace(icesee_kwargs.get('tinitial'), icesee_kwargs.get('num_years'), int((icesee_kwargs.get('num_years') - icesee_kwargs.get('tinitial'))/icesee_kwargs.get('timesteps_per_year'))+1),
        'nt': int((icesee_kwargs.get('num_years') - icesee_kwargs.get('tinitial'))/icesee_kwargs.get('timesteps_per_year')),
        'icesee_path': icesee_cwd,
        'data_path': icesee_kwargs.get('data_path'),
        'issm_dir': issm_dir,
        'issm_examples_dir': issm_examples_dir,
        'hpcmode': icesee_kwargs.get('hpcmode', False),
        'devmode': icesee_kwargs.get('devmode', False),
        'use_reference_data': icesee_kwargs.get('use_reference_data', False),
        'reference_data_dir': icesee_kwargs.get('reference_data_dir', 'data'),
        'reference_data': icesee_kwargs.get('reference_data'),
        'sill_friction': icesee_kwargs.get('sill_friction', 90000),
        'range_friction': icesee_kwargs.get('range_friction', 5000),
        'mean_friction': icesee_kwargs.get('mean_friction', 2500),
        'nugget_friction': icesee_kwargs.get('nugget_friction', 0),
        'sill_bed': icesee_kwargs.get('sill_bed', 4000),
        'range_bed': icesee_kwargs.get('range_bed', 50000),
        'nugget_bed': icesee_kwargs.get('nugget_bed', 200),
        'deepwater_melting_rate': float(icesee_kwargs.get('deepwater_melting_rate', 200)),
        'smb': float(icesee_kwargs.get('smb', 0.0)),
        'vel_idx': int(float(icesee_kwargs.get('vel_idx', 2))),
        'inversion_flag': icesee_kwargs.get('inversion_flag', False),
        'inversion_enabled': icesee_kwargs.get('inversion_flag', False),
        'inversion_start_time': float(
            icesee_kwargs.get('inversion_start_time', 0.0)
        ),
        'ensemble_spinup_dt': float(
            icesee_kwargs.get(
                'ensemble_spinup_dt',
                icesee_kwargs.get('timesteps_per_year'),
            )
        ),
        'ensemble_spinup_years': float(
            icesee_kwargs.get(
                'ensemble_spinup_years',
                icesee_kwargs.get('timesteps_per_year'),
            )
        ),
        'friction_idx': int(float(icesee_kwargs.get('friction_idx', 5))),
        'min_friction': float(icesee_kwargs.get('min_friction', 2000)),
        'max_friction': float(icesee_kwargs.get('max_friction', 4000)),
        'Nens': int(float(icesee_kwargs.get('Nens'))),
        'bed_relaxation_factor': float(icesee_kwargs.get('bed_relaxation_factor', 0.05)),
        'initial_bed_bias': float(icesee_kwargs.get('initial_bed_bias', 0.0015)),
        'initial_thickness_scale': float(icesee_kwargs.get('initial_thickness_scale', 1.0)),
        'initial_bed_offset_m': float(icesee_kwargs.get('initial_bed_offset_m', 0.0)),
        'initial_bed_background_domain': str(
            icesee_kwargs.get('initial_bed_background_domain', 'all')
        ),
        'initial_bed_gl_buffer_m': float(
            icesee_kwargs.get('initial_bed_gl_buffer_m', 0.0)
        ),
        'initial_floating_bed_anomaly_factor': float(
            icesee_kwargs.get('initial_floating_bed_anomaly_factor', 0.0)
        ),
        'initial_floating_bed_max_error_m': float(
            icesee_kwargs.get('initial_floating_bed_max_error_m', 100.0)
        ),
        'initial_floating_bed_transition_m': float(
            icesee_kwargs.get('initial_floating_bed_transition_m', 25000.0)
        ),
        'initial_floating_bed_flotation_margin_m': float(
            icesee_kwargs.get(
                'initial_floating_bed_flotation_margin_m', 5.0
            )
        ),
        'initial_bed_smoothing_iterations': int(
            icesee_kwargs.get('initial_bed_smoothing_iterations', 35)
        ),
        'initial_bed_smoothing_strength': float(
            icesee_kwargs.get('initial_bed_smoothing_strength', 0.65)
        ),
        'initial_bed_seed_max_x_m': float(
            icesee_kwargs.get('initial_bed_seed_max_x_m', 300000.0)
        ),
        'initial_bed_downstream_anomaly_factor': float(
            icesee_kwargs.get(
                'initial_bed_downstream_anomaly_factor', 0.60
            )
        ),
        'initial_thickness_anomaly_fraction': float(
            icesee_kwargs.get('initial_thickness_anomaly_fraction', 0.0)
        ),
        'initial_thickness_anomaly_m': float(
            icesee_kwargs.get('initial_thickness_anomaly_m', 0.0)
        ),
        'initial_thickness_delta_min_m': float(
            icesee_kwargs.get('initial_thickness_delta_min_m', -500.0)
        ),
        'initial_thickness_delta_max_m': float(
            icesee_kwargs.get('initial_thickness_delta_max_m', 500.0)
        ),
        'initial_floating_thickness_anomaly_factor': float(
            icesee_kwargs.get('initial_floating_thickness_anomaly_factor', 1.0)
        ),
        'initial_gl_seaward_thickness_m': float(
            icesee_kwargs.get('initial_gl_seaward_thickness_m', 0.0)
        ),
        'initial_gl_seaward_width_m': float(
            icesee_kwargs.get('initial_gl_seaward_width_m', 50000.0)
        ),
        'initial_bed_anomaly_m': float(
            icesee_kwargs.get('initial_bed_anomaly_m', 0.0)
        ),
        'initial_bed_delta_min_m': float(
            icesee_kwargs.get('initial_bed_delta_min_m', -500.0)
        ),
        'initial_bed_delta_max_m': float(
            icesee_kwargs.get('initial_bed_delta_max_m', 500.0)
        ),
        'initial_prior_length_x_m': float(
            icesee_kwargs.get('initial_prior_length_x_m', 120000.0)
        ),
        'initial_prior_length_y_m': float(
            icesee_kwargs.get('initial_prior_length_y_m', 40000.0)
        ),
        'initial_prior_pattern_phase': float(
            icesee_kwargs.get('initial_prior_pattern_phase', 0.0)
        ),
        'initial_thickness_factor_min': float(
            icesee_kwargs.get('initial_thickness_factor_min', 0.60)
        ),
        'initial_thickness_factor_max': float(
            icesee_kwargs.get('initial_thickness_factor_max', 1.25)
        ),
        'abs_vel_weight': float(icesee_kwargs.get('abs_vel_weight', 1.0)),
        'rel_vel_weight': float(icesee_kwargs.get('rel_vel_weight', 1.0)),
        'tikhonov_regularization_weight': float(icesee_kwargs.get('tikhonov_regularization_weight', 1e-13)),
        'b_nurge': float(icesee_kwargs.get('b_nurge', 0)),
        's_nurge': float(icesee_kwargs.get('s_nurge', 0)),
    })


def _initialize_rank_server_and_mesh(icesee_kwargs: dict, *, icesee_cwd: str,
                                      issm_examples_dir: str,
                                      ensemble_slot: int) -> tuple[object, int]:
    """Build this rank's server and initialize the shared mesh setup once."""

    ens_id = int(ensemble_slot)
    icesee_kwargs_file = f'icesee_kwargs_{ens_id}.mat'
    sio.savemat(os.path.join(issm_examples_dir, icesee_kwargs_file),
                matlab_icesee_kwargs(icesee_kwargs))

    shutil.copy(os.path.join(icesee_cwd, '..', '..', 'issm_utils', 'matlab2python', 'issm_env.m'), issm_examples_dir)
    shutil.copy(os.path.join(icesee_cwd, '..', '..', 'issm_utils', 'matlab2python', 'matlab_server.m'), issm_examples_dir)
    shutil.copy(os.path.join(icesee_cwd, 'Domain.exp'), issm_examples_dir)
    shutil.copy(os.path.join(icesee_cwd, icesee_kwargs.get('ParamFile')), issm_examples_dir)

    os.chdir(issm_examples_dir)

    server = MatlabServer(
        color=ens_id,
        Nens=icesee_kwargs['Nens'],
        comm=icesee_kwargs['icesee_comm'],
        verbose=icesee_kwargs.get('verbose'),
    )

    icesee_kwargs.update({
        'server': server, 'ens_id': ens_id,
        'rank': 0, 'nprocs': 1,
    })

    variable_size = initialize_model(**icesee_kwargs)
    nd = int(variable_size) * int(icesee_kwargs.get('total_state_param_vars'))

    os.chdir(icesee_cwd)
    return server, nd


def run_issm_execution_mode_3(**icesee_kwargs):
    """Drive ISMIP_Choi's round-robin mode-3 DA cycle. Registered under ``"issm"``."""

    icesee_kwargs = normalize_icesee_kwargs(icesee_kwargs)
    normalize_execution_mode(icesee_kwargs, expected=3)
    if icesee_kwargs.get("model_module") is None:
        from ICESEE.applications.issm_model.examples.ISMIP_Choi import _issm_enkf
        icesee_kwargs["model_module"] = _issm_enkf
    icesee_kwargs["vec_inputs_old"] = list(icesee_kwargs.get("vec_inputs", []))
    icesee_kwargs["observed_vars_params"] = list(icesee_kwargs.get("observed_vars", [])) + list(icesee_kwargs.get("observed_params", []))
    icesee_kwargs["all_observed"] = icesee_kwargs["observed_vars_params"]

    world = MPI.COMM_WORLD
    world_rank = int(world.Get_rank())

    global_start_time = MPI.Wtime()
    time_forecast_step = 0.0
    time_forecast_file_writing = 0.0

    topology = create_distributed_topology(world, spatial_ranks=1)
    ensemble_slot = int(topology.ensemble_slot)

    icesee_cwd = os.getcwd()
    if os.path.isabs(icesee_kwargs.get("data_path", "")):
        # The MATLAB bridge joins this value to the scientific source directory.
        icesee_kwargs["data_path"] = os.path.relpath(icesee_kwargs["data_path"],icesee_cwd)
    issm_dir = os.environ.get('ISSM_DIR')
    add_issm_dir_to_sys_path(issm_dir)
    issm_examples_dir = setup_example_directory(issm_dir, icesee_kwargs.get('example_name'))

    _build_static_config(
        icesee_kwargs,
        icesee_cwd=icesee_cwd,
        issm_dir=issm_dir,
        issm_examples_dir=issm_examples_dir,
    )
    icesee_kwargs.update({
        'icesee_comm': MPI.COMM_SELF,
        'comm': MPI.COMM_SELF,
        'model_nprocs': icesee_kwargs.get('model_nprocs'),
    })

    obs_t, obs_idx, num_observations = UtilsFunctions(icesee_kwargs).generate_observation_schedule(**icesee_kwargs)
    icesee_kwargs["obs_index"] = obs_idx
    icesee_kwargs["number_obs_instants"] = num_observations

    # --- every rank builds its own persistent server and bootstraps the
    # shared mesh setup exactly once (see module docstring). ---
    server, nd = _initialize_rank_server_and_mesh(
        icesee_kwargs,
        icesee_cwd=icesee_cwd,
        issm_examples_dir=issm_examples_dir,
        ensemble_slot=ensemble_slot,
    )
    icesee_kwargs['nd'] = nd

    Nens = int(icesee_kwargs["Nens"])
    nt = int(icesee_kwargs["nt"])

    from ICESEE.src.parallelization.forecast_recovery import ForecastRecoveryPolicy
    from ICESEE.src.parallelization.run_control import CycleStopController, atomic_json
    from ICESEE.src.parallelization.mode3_checkpointing import Mode3CheckpointManager, CheckpointPolicy
    from ICESEE.src.parallelization.distributed_native_runtime import restore_native_member_pool
    from ICESEE.applications.issm_model.issm_utils.recovery import commit_backend, recovery_identity, prune_backend
    recovery_policy = ForecastRecoveryPolicy.from_options(icesee_kwargs)
    managed = (recovery_policy is not None or bool(icesee_kwargs.get("issm_native_continuation", False))
               or bool(icesee_kwargs.get("resume_from_checkpoint", False))
               or icesee_kwargs.get("run_segment_cycles") is not None
               or icesee_kwargs.get("run_walltime_seconds") is not None)
    resume = bool(icesee_kwargs.get("resume_from_checkpoint", False))
    configured_recovery_id = recovery_identity(icesee_kwargs) if managed else None
    if managed and (icesee_kwargs.get("use_member_streaming") or icesee_kwargs.get("use_store_streaming_analysis")):
        raise NotImplementedError("ISSM backend recovery currently requires resident whole-member state")

    _modelrun_datasets = icesee_kwargs.get("data_path") or "_modelrun_datasets"
    os.environ["ICESEE_RESULTS_DIR"] = str(_modelrun_datasets)
    prepare_timestep_contract(world, _modelrun_datasets, time_grid=icesee_kwargs.get("t"))
    if world_rank == 0 and not os.path.exists(_modelrun_datasets):
        os.makedirs(_modelrun_datasets, exist_ok=True)
    world.Barrier()

    _true_nurged = f"{_modelrun_datasets}/true_nurged_states.h5"
    _synthetic_obs = f"{_modelrun_datasets}/synthetic_obs.h5"
    icesee_kwargs.update(
        {"true_nurged_file": _true_nurged, "synthetic_obs_file": _synthetic_obs}
    )

    # --- setup phase: real world root only. Its own ensemble slot is
    # always 0 (spatial_ranks=1), so its own server already has ens_id=0,
    # exactly the identity generate_true_wrong_state forces internally. ---
    true_wrong_time = 0.0
    observation_time = 0.0
    if world_rank == 0 and not resume:
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

    adapter = ISMIP_CHOI_NATIVE_ADAPTER
    pool = None if managed else initialize_native_member_pool(adapter, topology, icesee_kwargs)

    obs_indices = np.asarray(
        UtilsFunctions(icesee_kwargs).JObs_indices(nd), dtype=np.int64
    )
    error_mode = str(
        icesee_kwargs.get("enkf_observation_error_mode", "legacy_prior_anomalies")
    ).lower()
    if error_mode not in _SUPPORTED_ERROR_MODES:
        raise NotImplementedError(
            "execution_mode 3's ISSM runner supports 'legacy_prior_anomalies' "
            f"and 'stochastic_R' observation-error modes; got {error_mode!r}."
        )
    base_seed = int(icesee_kwargs.get("base_seed", 42))
    inversion_enabled = bool(icesee_kwargs.get(
        "inversion_enabled", icesee_kwargs.get("inversion_flag", False)
    ))
    inversion_start_time = float(icesee_kwargs.get("inversion_start_time", 0.0))
    km = 0

    run_id = str(icesee_kwargs.get("run_id", _DEFAULT_RUN_ID))
    history_root = os.path.join(_modelrun_datasets, "_mode3_state_history")
    initial_root = os.path.join(history_root, "initial_condition")
    steps_root = os.path.join(history_root, "steps")

    init_file_time = 0.0
    start_step = 0
    manager = None
    if managed:
        run_id += ":" + configured_recovery_id
        manager = Mode3CheckpointManager(history_root,run_id=run_id,topology=topology,
            number_of_members=Nens,nt=nt,
            analysis_steps=(np.asarray(obs_idx[:num_observations],dtype=int)-1).tolist()
                if icesee_kwargs.get("analysis_enabled", True) else [],
            policy=CheckpointPolicy.from_kwargs(icesee_kwargs),time_grid=icesee_kwargs["t"])
        if int(manager.policy.every) != 1:
            raise ValueError("ISSM managed continuation requires checkpoint_every=1")
        if resume:
            plan = manager.prepare_resume()
            from ICESEE.src.parallelization.distributed_native_cycle import _collective_native_call
            pool, restored = _collective_native_call(topology,
                lambda:restore_native_member_pool(adapter,topology,icesee_kwargs,
                    plan.checkpoint_path,expected_run_id=run_id),phase="ISSM backend restore")
            start_step, km = plan.start_step, plan.analyses_completed
        else:
            pool = initialize_native_member_pool(adapter,topology,icesee_kwargs)
            backend = commit_backend(pool,topology,icesee_kwargs,-1)
            manager.write_initial(pool.snapshot_owned(),extra={"issm_backend":backend})
    else:
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

    obs_index = icesee_kwargs["obs_index"]
    m_obs = icesee_kwargs["number_obs_instants"]

    completed_step = start_step
    with CycleStopController(icesee_kwargs, world, elapsed=MPI.Wtime()-global_start_time) as control:
        for k in range(start_step, nt):
            icesee_kwargs.update({"k": k, "km": km})
            do_analysis = bool(icesee_kwargs.get("analysis_enabled", True)) and observation_is_due(k, km, obs_index[:m_obs])
            batches = []
            apply_inversion = False
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
                        observation_variances=(sigma**2 if error_mode == "stochastic_r" else None),
                    )
                )

                cycle_time, model_cycle_time = resolve_analysis_cycle_time(
                    icesee_kwargs, k, km
                )
                apply_inversion = (
                    inversion_enabled
                    and cycle_time + 1.0e-12 >= inversion_start_time
                )
                icesee_kwargs["inversion_flag"] = apply_inversion
                if world_rank == 0 and inversion_enabled and not apply_inversion:
                    print(
                        "[ICESEE] Deferring friction inversion at "
                        f"observation t={cycle_time:g} yr "
                        f"(model t={model_cycle_time:g} yr); configured start is "
                        f"{inversion_start_time:g} yr."
                    )
                elif world_rank == 0 and apply_inversion:
                    print(
                        "[ICESEE] Friction inversion enabled at "
                        f"observation t={cycle_time:g} yr "
                        f"(model t={model_cycle_time:g} yr)."
                    )

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
                apply_inversion=apply_inversion,
            )
            time_forecast_step += MPI.Wtime() - _t
            if do_analysis:
                km += 1

            _t = MPI.Wtime()
            if manager is not None:
                backend = commit_backend(pool,topology,icesee_kwargs,k)
                manager.commit_step(k,result.local_analysis,did_analysis=do_analysis,
                    analyses_completed=km,extra={"issm_backend":backend,"did_inversion":apply_inversion})
                prune_backend(manager,topology,icesee_kwargs)
            else:
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
                        "did_inversion": apply_inversion,
                    },
                )
            time_forecast_file_writing += MPI.Wtime() - _t

            completed_step = k+1
            if managed and control.enabled and control.after_commit(completed_step-start_step):
                break
    world.Barrier()

    if managed:
        report_error = None
        if world_rank == 0:
            try:
                atomic_json(os.path.join(_modelrun_datasets,"_continuation","segment.json"),
                    dict(run_id=run_id,next_step=completed_step,nt=nt,
                         model_time=float(icesee_kwargs["t"][completed_step]),
                         analyses_completed=km,complete=completed_step==nt,
                         reason=control.reason or "complete"))
            except Exception as exc:
                report_error = f"{type(exc).__name__}: {exc}"
        report_error = world.bcast(report_error,root=0)
        if report_error:
            raise RuntimeError("ISSM continuation report failed: " + report_error)

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

    try:
        server.shutdown()
    except Exception:
        pass

    # ─────────────────────────────────────────────────────────────
    #  End Timer and Aggregate Elapsed Time Across Ranks
    # ─────────────────────────────────────────────────────────────
    global_elapsed_time = MPI.Wtime() - global_start_time

    register_run_metadata(
        execution_mode=icesee_kwargs.get("execution_mode"),
        model=icesee_kwargs.get("model_name"),
        forecast_steps=completed_step-start_step,
        analysis_events=km,
    )
    register_topology_run_metadata(topology, Nens)
    # Analyses are fused into the forecast step, so there is no separate
    # analysis timer: not measured when analyses ran, no events otherwise.
    emit_performance_report(
        world,
        elapsed_s=global_elapsed_time,
        phases={"true_wrong_state": true_wrong_time, "observation_generation": observation_time, "forecast_step": time_forecast_step, "init_file_io": init_file_time, "forecast_file_io": time_forecast_file_writing, "analysis_file_io": analysis_file_time, "analysis_step": None if km else 0.0},
        counts={"forecast_step": completed_step-start_step, "analysis_step": km},
        output_dir=_modelrun_datasets,
    )

    return icesee_kwargs


register_execution_mode_3(
    "issm",
    run_issm_execution_mode_3,
    notes=(
        "Round-robin, P_e-bounded mode-3 DA-cycle runner for ISMIP_Choi using "
        "ISMIP_CHOI_NATIVE_ADAPTER (a WholeMemberNativeAdapter, "
        "_issm_native.py) and members_for_ensemble_slot round-robin "
        "scheduling. Does NOT satisfy mode-3's strict no-rank-holds-a-"
        "complete-member invariant (ISSM's MATLAB/PETSc boundary only ever "
        "exposes a complete member); concurrently resident complete members "
        "are instead bounded by the configured ensemble-group count P_e "
        "(world size, since spatial_ranks=1), not Nens. See this module's "
        "docstring for the full rationale and the collective-barrier-safety "
        "argument."
    ),
)
