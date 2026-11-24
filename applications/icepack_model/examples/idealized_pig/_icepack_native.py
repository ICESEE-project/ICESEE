# @author: Brian Kyanjo
# Native Icepack adapter with persistent owned fields and spatial communicators.
from __future__ import annotations

from dataclasses import dataclass
import hashlib
from typing import Any, Mapping

import numpy as np
from firedrake import Function, SpatialCoordinate

import icepack
from ._melt_forcing import bind_melt_reference

from ICESEE.applications.icepack_model.examples.idealized_pig._icepack_model import (
    EXPERIMENT_WRONG,
    BasalMeltRate,
    estimated_melt_forecast, initial_estimated_melt, melt_forecast_mode,
    depth_law_melt, project_depth_law_melt,
    Icepack,
    initializeMesh,
    initializeRun,
    regViscosity,
    schoofFriction,
)
import ICESEE.applications.icepack_model.examples.idealized_pig.modelfunc as mf
from ICESEE.src.run_model_da._error_generation import coordinate_keyed_white_noise
# ICESEE.-prefixed (not the bare "applications...."/"src...." forms): see
# _distributed_fields.py's own import-path comment -- these two modules
# define classes (IcepackNativeState, DistributedStateLayout, ...) whose
# isinstance() identity depends on being imported via one consistent
# top-level package path everywhere in this stack. A real mode-3 run
# reproduced a "native member returned an unsupported distributed layout"
# TypeError traced to exactly this kind of split import.
from ICESEE.applications.icepack_model.icepack_utils._distributed_fields import (
    IcepackForecastFields,
    IcepackNativeAdapter,
    IcepackNativeState,
    scalar_ownership,
)
from ICESEE.src.parallelization._mpi_ensemble_intialization import (
    _member_initialization_context,
)


@dataclass
class _SharedPigContext:
    """Mesh, solver, and forcing fields shared by every member of one
    ensemble group (i.e. by every member scheduled on this rank)."""

    Q: Any
    V: Any
    forward_solver: Any
    bed: Any
    h0: Any
    s0: Any
    u0: Any
    A0: Any
    beta0: Any
    floating0: Any
    smb: Any
    coords: Any


_SHARED_CONTEXT_CACHE: dict[int, _SharedPigContext] = {}


def _build_shared_context(
    topology: Any, icesee_kwargs: Mapping[str, Any]
) -> _SharedPigContext:
    """Build the mesh/solver/forcing fields once for this ensemble group.

    Mirrors ``_icepack_model.initialize_model`` exactly, except the mesh
    communicator is ``topology.spatial_comm`` -- only the ranks inside this
    member's ensemble group -- instead of the world communicator, and the
    outputs are kept as persistent Firedrake objects instead of being
    packed into a flat per-member state vector.
    """

    from ICESEE.src.utils.init_trace import init_trace_enabled, trace_init

    mesh_kwargs = dict(icesee_kwargs)
    mesh_kwargs["comm"] = topology.spatial_comm
    if init_trace_enabled():
        trace_init("ENTER per-group model context",
                   group=getattr(topology, "ensemble_slot", None),
                   model_rank=getattr(topology, "spatial_rank", None),
                   model_comm=topology.spatial_comm)
    mesh, _meshOpts, Q, V = initializeMesh(**mesh_kwargs)

    forward_model = icepack.models.IceStream(
        friction=schoofFriction, viscosity=regViscosity
    )
    opts = {
        "dirichlet_ids": [1],
        "diagnostic_solver_parameters": {"max_iterations": 150, "tolerance": 1e-6},
    }
    forward_solver = icepack.solvers.FlowSolver(forward_model, **opts)

    _h, h0, _s, s0, u0, bed, _zF, _grounded0, floating0, A0, beta0, smb = (
        initializeRun(mesh_kwargs, forward_solver, mesh, Q, V)
    )

    # This rank's OWNED physical mesh-node coordinates only -- O(Nx/P_model)
    # memory, never a global-sized array. Computed once per ensemble group
    # (cached alongside the rest of the shared context) and used for the
    # restart layout fingerprint and the reference-numbering map.
    coords = node_coordinates(Q, mesh)

    if init_trace_enabled():
        trace_init("EXIT  per-group model context",
                   group=getattr(topology, "ensemble_slot", None),
                   model_rank=getattr(topology, "spatial_rank", None),
                   owned_nodes=coords.shape[0])
    return _SharedPigContext(
        Q=Q,
        V=V,
        forward_solver=forward_solver,
        bed=bed,
        h0=h0,
        s0=s0,
        u0=u0,
        A0=A0,
        beta0=beta0,
        floating0=floating0,
        smb=smb,
        coords=coords,
    )


def node_coordinates(Q: Any, mesh: Any = None) -> np.ndarray:
    """Owned nodal coordinates of scalar space ``Q``, in owned-row order."""

    x_expr, y_expr = SpatialCoordinate(Q.mesh() if mesh is None else mesh)
    return np.stack(
        (
            np.asarray(Function(Q).interpolate(x_expr).dat.data_ro),
            np.asarray(Function(Q).interpolate(y_expr).dat.data_ro),
        ),
        axis=1,
    )


def owned_node_layout(*, topology: Any, icesee_kwargs: Mapping[str, Any]):
    """``(owned_start, global_nodes, owned_coords)`` of this rank's scalar
    rows in the distributed member layout (see build_icepack_field_registry)."""

    ctx = _shared_context(topology, icesee_kwargs)
    ownership = scalar_ownership(ctx.h0)
    return ownership.owned_start, ownership.global_size, ctx.coords


def _shared_context(
    topology: Any, icesee_kwargs: Mapping[str, Any]
) -> _SharedPigContext:
    key = id(topology)
    context = _SHARED_CONTEXT_CACHE.get(key)
    if context is None:
        context = _build_shared_context(topology, icesee_kwargs)
        _SHARED_CONTEXT_CACHE[key] = context
    return context


def release_model_context(topology):
    """Release the group's cached Firedrake objects before its MPI communicator."""
    _SHARED_CONTEXT_CACHE.pop(id(topology), None)


def _basal_melt_for_step(
    icesee_kwargs: Mapping[str, Any],
    step: int,
    floating: Any,
    Q: Any,
    s: Any,
    h: Any,
    *,
    experiment: str = EXPERIMENT_WRONG,
) -> Any:
    """Evaluate the configured true/wrong basal-melt scenario at this step.

    Ensemble forecasts default to EXPERIMENT_WRONG. Truth generation uses the
    single-rank reference generator with EXPERIMENT_TRUE."""

    dt = float(icesee_kwargs["dt"])
    if step < (6 / dt):
        scenario = "control"
    elif (6 / dt) <= step < (15 / dt):
        scenario = "warm"
    elif (15 / dt) <= step < (18 / dt):
        scenario = "control"
    elif (18 / dt) <= step < (20 / dt):
        scenario = "warm"
    elif (20 / dt) <= step < (25 / dt):
        scenario = "control"
    elif (25 / dt) <= step < (27 / dt):
        scenario = "warm"
    elif (27 / dt) <= step < (31 / dt):
        scenario = "control"
    elif (31 / dt) <= step < (40 / dt):
        scenario = "warm"
    elif (40 / dt) <= step < (48 / dt):
        scenario = "control"
    elif (48 / dt) <= step < (50 / dt):
        scenario = "warm"
    elif (50 / dt) <= step < (59 / dt):
        scenario = "control"
    elif (59 / dt) <= step < (64 / dt):
        scenario = "warm"
    elif (64 / dt) <= step < (69 / dt):
        scenario = "control"
    elif (69 / dt) <= step < (76 / dt):
        scenario = "warm"
    else:
        scenario = "control"
    basal_melt_field, _melt_max = BasalMeltRate(
        icesee_kwargs, step, floating, Q, s, h, scenario=scenario, experiment=experiment
    )
    return basal_melt_field


def _basal_melt_in_state(icesee_kwargs: Mapping[str, Any]) -> bool:
    """Return whether vec_inputs includes basal_melt_field, independent of
    joint_estimation bookkeeping."""

    return "basal_melt_field" in icesee_kwargs.get("vec_inputs", ())


def _recompute_floating(ctx: "_SharedPigContext", surface: Any) -> Any:
    """Recompute the floating mask from the CURRENT surface, mirroring
    ``_icepack_model.Icepack``'s own internal recomputation. floating/
    grounded are pure diagnostics of (surface, bed) -- never persistent
    DA state -- so they are always recomputed fresh here rather than
    cached from a frozen initial value (see ``forecast_member``)."""

    zF = mf.flotationHeight(ctx.bed, ctx.Q)
    floating, _grounded = mf.flotationMask(surface, zF, ctx.Q)
    return floating


def _icepack_run_kwargs(
    ctx: _SharedPigContext, icesee_kwargs: Mapping[str, Any]
) -> dict[str, Any]:
    """Small, bounded scalar/reference-field kwargs ``Icepack()`` needs.

    Built fresh per call (never proportional to mesh size) so the shared
    application-level ``icesee_kwargs`` mapping is never mutated by
    ``Icepack()``'s internal ``icesee_kwargs["h"] = h`` side effect.
    """

    from firedrake import Constant
    threshold = icesee_kwargs["uThresh"]
    if not hasattr(threshold, "ufl_shape"):
        threshold = Constant(threshold)
    return {
        "Q": ctx.Q,
        "s0": ctx.s0,
        "beta0": ctx.beta0,
        "A0": ctx.A0,
        "bed": ctx.bed,
        "uThresh": threshold,
        "GLThresh": icesee_kwargs["GLThresh"],
        "hThresh": icesee_kwargs["hThresh"],
        "water_to_ice": icesee_kwargs.get("water_to_ice", 1.0),
    }


# Apply parameter nudging only when joint estimation is enabled.
def _nudge_basal_melt(basal_melt_field: Any, icesee_kwargs: Mapping[str, Any], Q: Any) -> Any:
    """Joint-estimation-specific nudge: initializes the ensemble's WRONG
    prior estimate of a jointly-estimated parameter. Same operation as
    ``_icepack_enkf.initialize_ensemble``."""

    nudged = basal_melt_field.dat.data_ro + icesee_kwargs["wrong_basal_melt_field"]
    result = Function(Q)
    result.dat.data[:] = nudged
    return result


def initialize_member(
    member_id: int, *, topology: Any, icesee_kwargs: Mapping[str, Any]
) -> IcepackNativeState:
    """Copy the unperturbed physical fields at the initial model time."""

    ctx = _shared_context(topology, icesee_kwargs)
    icesee_kwargs = bind_melt_reference(icesee_kwargs, ctx)

    h_perturbed = ctx.h0.dat.data_ro.copy()

    h_p = Function(ctx.Q)
    h_p.dat.data[:] = h_perturbed

    # Field existence: does THIS state carry a basal_melt_field slot at
    # all? Driven by the application's configured state definition
    # (vec_inputs), never by joint_estimation.
    basal_melt_in_state = _basal_melt_in_state(icesee_kwargs)
    # Bootstrap floating mask (ctx.floating0): there is no evolved member
    # state yet, so this mirrors _icepack_enkf.initialize_ensemble's own
    # use of the bootstrap floating/h0/s0 fields -- not the frozen-forever
    # ctx.floating0 that forecast_member used to (incorrectly) keep reusing.
    basal_melt_field = _basal_melt_for_step(
        icesee_kwargs, 0, ctx.floating0, ctx.Q, ctx.s0, ctx.h0,
        experiment=EXPERIMENT_WRONG,
    )
    # Joint-estimation algorithm: initialize the ensemble's WRONG prior
    # belief about the jointly-estimated parameter (see _nudge_basal_melt's
    # own docstring). Independent of whether basal_melt_field exists in
    # this state -- if joint_estimation is on, the field always exists too
    # (it IS the jointly-estimated parameter), so no ordering hazard here.
    if melt_forecast_mode(icesee_kwargs) in {"estimated_field", "estimated_depth_law"}:
        basal_melt_field = initial_estimated_melt(basal_melt_field, member_id, icesee_kwargs)
    elif icesee_kwargs.get("joint_estimation"):
        basal_melt_field = _nudge_basal_melt(basal_melt_field, icesee_kwargs, ctx.Q)

    h = h_p
    u = ctx.u0.copy(deepcopy=True)
    s = ctx.s0.copy(deepcopy=True)

    return IcepackNativeState(
        thickness=h,
        velocity=u,
        surface=s,
        basal_melt=basal_melt_field if basal_melt_in_state else None,
        layout_id="icepack-idealized-pig-v1",
    )


def forecast_member(
    state: IcepackNativeState,
    timestep: int,
    *,
    topology: Any,
    icesee_kwargs: Mapping[str, Any],
    _local_dt: float | None = None,
) -> IcepackForecastFields:
    """Advance persistent native fields by one step.

    Use the shared reference thickness for inflow boundaries. Recompute the
    floating mask from the member's current surface and bed each step; this
    diagnostic mask is not independently checkpointed assimilation state."""

    from ICESEE.src.parallelization.forecast_recovery import ForecastRecoveryPolicy
    policy = ForecastRecoveryPolicy.from_options(icesee_kwargs)
    if policy is not None and _local_dt is None:
        return _recover_forecast_member(state, timestep, topology=topology,
                                        icesee_kwargs=icesee_kwargs, policy=policy)
    ctx = _shared_context(topology, icesee_kwargs)
    icesee_kwargs = bind_melt_reference(icesee_kwargs, ctx)
    step = int(timestep) if _local_dt is None else float(timestep)
    local_dt = float(icesee_kwargs["dt"]) if _local_dt is None else float(_local_dt)

    floating = _recompute_floating(ctx, state.surface)
    mode = melt_forecast_mode(icesee_kwargs)
    if mode == "estimated_depth_law":
        basal_melt_field, amplitude = project_depth_law_melt(
            state.basal_melt, floating, ctx.Q, state.surface, state.thickness, step, icesee_kwargs)
    elif mode == "estimated_field":
        basal_melt_field = estimated_melt_forecast(
            state.basal_melt, floating, ctx.Q, icesee_kwargs)
    else:
        basal_melt_field = _basal_melt_for_step(
            icesee_kwargs, step, floating, ctx.Q, state.surface, state.thickness,
            experiment=EXPERIMENT_WRONG,
        )
    run_kwargs = _icepack_run_kwargs(ctx, icesee_kwargs)
    h, u, s, _floating, _grounded = Icepack(
        ctx.forward_solver,
        state.thickness,
        state.velocity,
        ctx.smb,
        basal_melt_field,
        ctx.bed,
        local_dt,
        ctx.h0,
        run_kwargs,
    )

    if mode == "estimated_depth_law":
        from ._melt_forcing import advance_melt_amplitude
        amplitude = advance_melt_amplitude(amplitude, 1.0, local_dt, icesee_kwargs)
        next_floating = _recompute_floating(ctx, s)
        basal_melt_field = depth_law_melt(amplitude, next_floating, ctx.Q, s, h,
            step + local_dt / float(icesee_kwargs["dt"]), icesee_kwargs)
    return IcepackForecastFields(
        thickness=h,
        velocity=u,
        surface=s,
        basal_melt=basal_melt_field if state.basal_melt is not None else None,
    )


def _recover_forecast_member(state, timestep, *, topology, icesee_kwargs, policy):
    """Retry physical solves before process noise, profiles or analysis run."""
    import json
    from firedrake.exceptions import ConvergenceError
    from ICESEE.src.parallelization.forecast_recovery import (
        advance_with_recovery, RecoverableForecastError,
    )
    ctx = _shared_context(topology, icesee_kwargs)
    result = [None]
    macro_dt = float(icesee_kwargs["dt"])

    def restore(snapshot):
        state.registry.unpack_owned(snapshot, synchronize=True)
        # Solver contexts can retain failed iterates; rebuild with the same options.
        model = icepack.models.IceStream(friction=schoofFriction, viscosity=regViscosity)
        ctx.forward_solver = icepack.solvers.FlowSolver(model, dirichlet_ids=[1],
            diagnostic_solver_parameters={"max_iterations": 150, "tolerance": 1e-6})

    def advance(time, dt):
        try:
            substep = round((time - int(timestep) * macro_dt) / dt)
            forcing_step = int(timestep) + substep * (dt / macro_dt)
            result[0] = forecast_member(state, forcing_step, topology=topology,
                icesee_kwargs=icesee_kwargs, _local_dt=dt)
            state.update_from_forecast(result[0])
        except ConvergenceError as exc:
            raise RecoverableForecastError(str(exc)) from exc

    def validate():
        if not np.isfinite(state.registry.pack_owned()).all():
            raise RecoverableForecastError("non-finite physical forecast")

    def report(event):
        if topology.spatial_comm.rank == 0:
            print("[ICESEE] forecast_recovery " + json.dumps(dict(event,
                timestep=int(timestep), ensemble_slot=topology.ensemble_slot)), flush=True)

    advance_with_recovery(start_time=int(timestep) * macro_dt, interval=macro_dt,
        policy=policy, capture=lambda: state.registry.pack_owned().copy(),
        restore=restore, advance=advance, validate=validate,
        comm=topology.spatial_comm, report=report)
    return result[0]


def reactivate_member(
    member_id: int,
    packed_state: np.ndarray,
    *,
    topology: Any,
    icesee_kwargs: Mapping[str, Any],
) -> IcepackNativeState:
    """Rebuild a live native member from a previously packed owned-state
    array -- Option B's bounded-memory round scheduling
    (``StreamingNativeDistributedMemberPool``,
    src/parallelization/distributed_streaming_runtime.py) calls this
    instead of ``initialize_member`` for every activation after the first.

    Unlike ``initialize_member``, this performs NO random perturbation and
    NO Icepack forward solve -- it only reconstructs empty Firedrake
    Functions on the SAME shared, cached ``_SharedPigContext`` (mesh/Q/V are
    never rebuilt; see ``_shared_context``'s per-ensemble-group cache) and
    writes ``packed_state`` into them via the existing, already-tested
    ``DistributedFieldRegistry.unpack_owned``. This must exactly reproduce
    whatever state ``pack_owned()`` last extracted for this member --
    verified by ``test_icepack_native_member_streaming.py``'s
    deactivate/reactivate round-trip tests.
    """

    state = allocate_member(member_id, topology=topology, icesee_kwargs=icesee_kwargs)
    state.registry.unpack_owned(
        np.asarray(packed_state, dtype=np.float64), synchronize=False
    )
    return state


def allocate_member(
    member_id: int,
    *,
    topology: Any,
    icesee_kwargs: Mapping[str, Any],
) -> IcepackNativeState:
    """Empty native storage for one member on the shared, cached context.

    No perturbation and no Icepack solve: a restart allocates every member
    this way and then overwrites all owned values from the checkpoint, so
    resuming never regenerates the initial ensemble.
    """

    ctx = _shared_context(topology, icesee_kwargs)
    return IcepackNativeState(
        thickness=Function(ctx.Q),
        velocity=Function(ctx.V),
        surface=Function(ctx.Q),
        basal_melt=Function(ctx.Q) if _basal_melt_in_state(icesee_kwargs) else None,
        layout_id="icepack-idealized-pig-v1",
    )


def layout_fingerprint(*, topology: Any, icesee_kwargs: Mapping[str, Any]) -> str:
    """Identity of this rank's owned mesh-node numbering.

    A digest of the owned nodal coordinates in owned-DOF order (rounded to
    1 mm), so a restart refuses checkpoint rows written under a different
    mesh numbering instead of silently scattering them onto wrong nodes.
    """

    ctx = _shared_context(topology, icesee_kwargs)
    coords = np.ascontiguousarray(np.round(ctx.coords, 3), dtype=np.float64)
    digest = hashlib.sha256(coords.tobytes())
    digest.update(str(coords.shape).encode())
    return "coords-sha256:" + digest.hexdigest()[:32]


class _IdealizedPigAdapter(IcepackNativeAdapter):
    """Keep application diagnostics outside the model-agnostic runtime."""

    def apply_native_process_noise(self, member, owned, increment, timestep, *, topology, icesee_kwargs):
        if icesee_kwargs.get("basal_melt_noise_basis", "field") == "depth_amplitude":
            state = self._state(member)
            ctx = _shared_context(topology, icesee_kwargs)
            icesee_kwargs = bind_melt_reference(icesee_kwargs, ctx)
            floating = _recompute_floating(ctx, state.surface)
            basis = depth_law_melt(1., floating, ctx.Q, state.surface,
                                   state.thickness, int(timestep)+1, icesee_kwargs)
            increment = np.asarray(increment).copy()
            for name, rows, _ in member.fields.layout.iter_owned_global_intervals():
                if name == "basal_melt_field":
                    increment[rows] *= basis.dat.data_ro
        member.unpack_owned(owned + increment)

    def finalize_native_analysis(self, member, forecast_owned, timestep, *, topology, icesee_kwargs):
        options = dict(icesee_kwargs, analysis_member_id=member.member_id)
        return super().finalize_native_analysis(member, forecast_owned, timestep,
            topology=topology, icesee_kwargs=options)

    def finalize_forecast_noise(self, member, timestep, *, topology, icesee_kwargs):
        from ._melt_forcing import analysis_state_policy
        if analysis_state_policy(icesee_kwargs) != "diagnostic_balance":
            return
        state = self._state(member)
        names = list(icesee_kwargs["vec_inputs"])[:int(icesee_kwargs["num_state_vars"])]
        selected = icesee_kwargs.get("process_noise_vars", names)
        sigmas = icesee_kwargs["sig_Q"]
        geometry_changed = any(name in selected and name in {"h", "s", "u", "v"}
                               and i < len(sigmas) and float(sigmas[i]) != 0
                               for i, name in enumerate(names))
        if geometry_changed:
            finalize_analysis_state(state, None, timestep,
                topology=topology, icesee_kwargs=dict(icesee_kwargs, analysis_member_id=member.member_id), force=True)
        else:
            # Melt-only noise does not change the diagnostic velocity problem.
            ctx = _shared_context(topology, icesee_kwargs)
            icesee_kwargs = bind_melt_reference(icesee_kwargs, ctx)
            floating = _recompute_floating(ctx, state.surface)
            melt = state.basal_melt
            mode = melt_forecast_mode(icesee_kwargs)
            if mode == "estimated_depth_law":
                melt, _ = project_depth_law_melt(melt, floating, ctx.Q, state.surface,
                    state.thickness, int(timestep) + 1, icesee_kwargs)
            elif mode == "estimated_field":
                melt = estimated_melt_forecast(melt, floating, ctx.Q, icesee_kwargs)
            state.update_from_forecast(IcepackForecastFields(thickness=state.thickness,
                velocity=state.velocity, surface=state.surface, basal_melt=melt))
        state.registry.synchronize_ghosts()

    def initialize_native_member(self, member_id, *, topology, icesee_kwargs):
        member = super().initialize_native_member(
            member_id, topology=topology, icesee_kwargs=icesee_kwargs
        )
        from ._mode3_profiles import write_member_profile
        ctx = _shared_context(topology, icesee_kwargs)
        write_member_profile(
            member_id, ctx.h0, ctx.s0, ctx.bed, step=0,
            comm=topology.spatial_comm, icesee_kwargs=icesee_kwargs,
        )
        return member

    def forecast_native_member(self, member, timestep, *, topology, icesee_kwargs):
        super().forecast_native_member(
            member, timestep, topology=topology, icesee_kwargs=icesee_kwargs
        )
        from ._mode3_profiles import write_member_profile
        ctx = _shared_context(topology, icesee_kwargs)
        state = self._state(member)
        write_member_profile(
            member.member_id, state.thickness, state.surface, ctx.bed,
            step=int(timestep) + 1, comm=topology.spatial_comm,
            icesee_kwargs=icesee_kwargs,
        )


def finalize_analysis_state(state, forecast_owned, timestep, *, topology, icesee_kwargs, force=False):
    """Restore admissible geometry and diagnostic velocity after an analysis."""
    from ._melt_forcing import analysis_state_policy
    from ICESEE.src.utils.timestep_contract import observation_cycle_steps
    if analysis_state_policy(icesee_kwargs) != "diagnostic_balance":
        return
    indices = np.asarray(icesee_kwargs.get("obs_index", []), dtype=np.int64)
    count = int(icesee_kwargs.get("number_obs_instants", indices.size))
    if not force and int(timestep) not in observation_cycle_steps(indices[:count]):
        return
    ctx = _shared_context(topology, icesee_kwargs)
    icesee_kwargs = bind_melt_reference(icesee_kwargs, ctx)
    raw = state.thickness.dat.data_ro
    floor = float(icesee_kwargs["hThresh"])
    local = (bool(np.isfinite(raw).all()), float(raw.min()) if raw.size else np.inf,
             float(np.maximum(floor - raw, 0).max()) if raw.size else 0.)
    comm = topology.spatial_comm if topology is not None else None
    metrics = comm.allgather(local) if comm is not None else [local]
    if not all(item[0] for item in metrics):
        raise FloatingPointError("non-finite analyzed thickness")
    correction = max(item[2] for item in metrics)
    limit = icesee_kwargs.get("pig_max_thickness_projection")
    if limit is not None and correction > float(limit):
        raise RuntimeError(f"thickness projection {correction:g} exceeds pig_max_thickness_projection={limit}")
    h = state.thickness.copy(deepcopy=True)
    h.dat.data[:] = np.maximum(float(icesee_kwargs["hThresh"]), h.dat.data_ro)
    s = icepack.compute_surface(thickness=h, bed=ctx.bed)
    zF = mf.flotationHeight(ctx.bed, ctx.Q)
    floating, grounded = mf.flotationMask(s, zF, ctx.Q)
    scale = mf.reduceNearGLBeta(s, ctx.s0, zF, grounded, ctx.Q, icesee_kwargs["GLThresh"])
    beta = icepack.interpolate(ctx.beta0 * scale, ctx.Q)
    u = ctx.forward_solver.diagnostic_solve(
        velocity=state.velocity.copy(deepcopy=True), thickness=h, surface=s,
        beta=beta, fluidity=ctx.A0, uThresh=icesee_kwargs["uThresh"],
        floating=floating, grounded=grounded)
    melt = state.basal_melt
    if melt_forecast_mode(icesee_kwargs) == "estimated_depth_law":
        melt, _ = project_depth_law_melt(melt, floating, ctx.Q, s, h, int(timestep) + 1, icesee_kwargs)
    elif melt_forecast_mode(icesee_kwargs) == "estimated_field":
        melt = estimated_melt_forecast(melt, floating, ctx.Q, icesee_kwargs)
    if icesee_kwargs.get("pig_report_analysis", False):
        import json
        from firedrake import assemble, dx
        volume_added = float(assemble((h - state.thickness) * dx))
        if comm is None or comm.rank == 0:
            print("[ICESEE] physical_analysis " + json.dumps(dict(
                member=icesee_kwargs.get("analysis_member_id"), timestep=int(timestep),
                boundary="noise" if force else "analysis", min_thickness_before=min(item[1] for item in metrics),
                max_thickness_projection=correction, projection_volume_added=volume_added)), flush=True)
    state.update_from_forecast(IcepackForecastFields(
        thickness=h, velocity=u, surface=s, basal_melt=melt))


IDEALIZED_PIG_NATIVE_ADAPTER = _IdealizedPigAdapter(
    initialize_member=initialize_member,
    forecast_member=forecast_member,
    finalize_analysis=finalize_analysis_state,
    reactivate_member=reactivate_member,
    allocate_member=allocate_member,
    layout_fingerprint=layout_fingerprint,
)
