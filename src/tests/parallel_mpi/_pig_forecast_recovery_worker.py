"""PACE physical recovery validation; requires real PIG inputs and Firedrake."""
from pathlib import Path
import json
import os
import sys

SOURCE = Path(__file__).resolve().parents[3]
PIG = SOURCE / "applications/icepack_model/examples/idealized_pig"
sys.path[:0] = [str(SOURCE.parent), str(SOURCE), str(PIG)]
os.chdir(PIG)

from ICESEE.config._utility_imports import icesee_kwargs
import numpy as np
from mpi4py import MPI
from firedrake.exceptions import ConvergenceError
from ICESEE.src.utils.pyop2_compile_cache import install_consistent_compile_cache
from ICESEE.src.parallelization.distributed_topology import create_distributed_topology
from ICESEE.src.parallelization.distributed_runtime import LocalMemberEnsemble
from ICESEE.src.parallelization.distributed_checkpoint import load_distributed_checkpoint
from ICESEE.src.parallelization.mode3_checkpointing import Mode3CheckpointManager, CheckpointPolicy
from ICESEE.src.parallelization.run_control import CycleStopController, atomic_json
from ICESEE.applications.icepack_model.examples.idealized_pig import _icepack_native as native
from ICESEE.applications.icepack_model.examples.idealized_pig._melt_forcing import configure_melt_experiment


def main():
    install_consistent_compile_cache()
    world = MPI.COMM_WORLD
    options = dict(icesee_kwargs)
    options.update(dt=float(options["timesteps_per_year"]), nt=2,
                   run_stop_time=None, save_steps=[])
    configure_melt_experiment(options)
    topology = create_distributed_topology(world, spatial_ranks=2)
    member_id = topology.ensemble_slot
    state = native.initialize_member(member_id, topology=topology, icesee_kwargs=options)
    ctx = native._shared_context(topology, options)
    initial = state.registry.pack_owned().copy()
    layout = state.registry.layout
    reports = []

    def restore(values):
        state.registry.unpack_owned(values, synchronize=True)
        model = native.icepack.models.IceStream(friction=native.schoofFriction,
                                               viscosity=native.regViscosity)
        ctx.forward_solver = native.icepack.solvers.FlowSolver(model, dirichlet_ids=[1],
            diagnostic_solver_parameters={"max_iterations":150,"tolerance":1e-6})

    def step(k, cfg, local_dt=None):
        result = native.forecast_member(state,k,topology=topology,
            icesee_kwargs=cfg,_local_dt=local_dt)
        state.update_from_forecast(result)
        return state.registry.pack_owned().copy()

    def compare(label, expected, actual):
        difference = float(world.allreduce(float(np.max(np.abs(expected-actual),initial=0)),op=MPI.MAX))
        exact = bool(world.allreduce(bool(np.array_equal(expected,actual)),op=MPI.LAND))
        reports.append(dict(check=label,bitwise_equal=exact,max_absolute_difference=difference))
        if not exact:
            raise AssertionError(f"{label}: expected identical operations; max difference={difference}")

    recovery = dict(options,forecast_recovery=dict(max_halvings=3,min_dt=options["dt"]/8,max_substeps=8))
    restore(initial)
    baseline = step(0,options)
    restore(initial)
    accepted = step(0,recovery)
    compare("no_rejection_vs_original",baseline,accepted)

    restore(initial)
    step(0,options,options["dt"]/2)
    half_step_reference = step(.5,options,options["dt"]/2)
    restore(initial)
    original = native.Icepack
    attempts = []
    def inject(*args,**kwargs):
        attempts.append(float(args[6]))
        if len(attempts) == 1:
            # Every rank in this model group rejects before model collectives.
            state.thickness.dat.data[:] += 700
            raise ConvergenceError("physical recovery validation: injected rejection")
        return original(*args,**kwargs)
    native.Icepack = inject
    try:
        retried = step(0,recovery)
    finally:
        native.Icepack = original
    compare("rejected_trial_vs_explicit_half_steps",half_step_reference,retried)
    assert attempts == [options["dt"],options["dt"]/2,options["dt"]/2]

    restore(initial)
    step(0,recovery)
    uninterrupted = step(1,recovery)
    restore(initial)
    manager = Mode3CheckpointManager(Path(options["data_path"])/"physical_history",
        run_id="pig-recovery-validation-v1",topology=topology,
        number_of_members=topology.ensemble_groups,nt=2,analysis_steps=[],
        policy=CheckpointPolicy(every=1,keep_last=0),time_grid=[0,options["dt"],2*options["dt"]])
    def snapshot():
        return LocalMemberEnsemble(layout,{member_id:state.registry.pack_owned().copy()})
    manager.write_initial(snapshot())
    with CycleStopController(dict(run_segment_cycles=1,checkpoint_every=1),world) as control:
        step(0,recovery)
        manager.commit_step(0,snapshot(),did_analysis=False,analyses_completed=0)
        assert control.after_commit(1)
    plan = manager.prepare_resume()
    assert plan.start_step == 1 and plan.analyses_completed == 0
    state.thickness.dat.data[:] = -888
    restored,_ = load_distributed_checkpoint(plan.checkpoint_path,layout,topology,
        expected_run_id=manager.run_id,number_of_members=topology.ensemble_groups)
    restore(restored.members[member_id])
    resumed = step(plan.start_step,recovery)
    manager.commit_step(1,snapshot(),did_analysis=False,analyses_completed=0)
    compare("physical_checkpoint_resume",uninterrupted,resumed)
    minimum_h = world.allreduce(float(np.min(state.thickness.dat.data_ro,initial=np.inf)),op=MPI.MIN)
    maximum_h = world.allreduce(float(np.max(state.thickness.dat.data_ro,initial=-np.inf)),op=MPI.MAX)
    assert minimum_h >= float(options["hThresh"]) - 1e-10
    assert world.allreduce(bool(np.isfinite(resumed).all()),op=MPI.LAND)
    if world.rank == 0:
        report = dict(passed=True,ranks=world.size,spatial_ranks=2,
            members=topology.ensemble_groups,dt=options["dt"],checks=reports,
            h_min=minimum_h,h_max=maximum_h,solver_settings_changed=False,
            process_noise_tested=False,analysis_events_tested=False)
        atomic_json(Path(options["data_path"]).parent/"physical-verification.json",report)
        print(json.dumps(report,indent=2),flush=True)


if __name__ == "__main__":
    try:
        main()
    except Exception:
        import traceback
        traceback.print_exc()
        MPI.COMM_WORLD.Abort(1)
