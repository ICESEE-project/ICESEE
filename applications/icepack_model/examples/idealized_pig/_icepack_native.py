# ==============================================================================
# @des: Execution-mode-3 native adapter wiring for the idealized PIG example.
#
#       Bridges the existing whole-array idealized_pig physics in
#       ``_icepack_model.py``/``_icepack_enkf.py`` (mesh setup, the
#       prognostic/diagnostic Icepack step, and the basal-melt forcing
#       schedule) to the mode-3 distributed runtime without ever
#       materializing a full-mesh numpy vector for one ensemble member.
#       Every persistent field stays a Firedrake ``Function`` whose owned
#       degrees of freedom are read/written directly through
#       ``IcepackNativeState``/``DistributedFieldRegistry`` -- see
#       applications/icepack_model/icepack_utils/_distributed_fields.py.
#
#       Mesh, function spaces, the flow solver, and the forcing/boundary
#       fields (bed, SMB, fluidity, friction, inflow thickness) are
#       identical for every ensemble member of one run -- only the
#       perturbed thickness/velocity/surface/basal-melt fields differ.
#       They are therefore built once per process (a rank belongs to
#       exactly one ensemble slot for its whole lifetime under the mode-3
#       process grid) and shared by every member this rank forecasts,
#       instead of being duplicated per member.  The mesh is partitioned
#       across ``topology.spatial_comm`` (only the ranks inside this
#       member's ensemble group), not the world communicator, so each
#       ensemble group's per-rank memory is bounded by its own share of one
#       member's mesh -- never by the full ensemble or the full world size.
# @date: 2026-08-24
# ==============================================================================

from __future__ import annotations

from dataclasses import dataclass
import hashlib
from typing import Any, Mapping

import numpy as np
from firedrake import Function, SpatialCoordinate

import icepack

from ICESEE.applications.icepack_model.examples.idealized_pig._icepack_model import (
    EXPERIMENT_WRONG,
    BasalMeltRate,
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

    mesh_kwargs = dict(icesee_kwargs)
    mesh_kwargs["comm"] = topology.spatial_comm
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
    # (cached alongside the rest of the shared context) and reused by every
    # member's coordinate-keyed initial perturbation (see initialize_member),
    # so a genuine per-DOF physical identity is available without
    # regenerating or broadcasting coordinates per member.
    x_expr, y_expr = SpatialCoordinate(mesh)
    coords = np.stack(
        (
            np.asarray(Function(Q).interpolate(x_expr).dat.data_ro),
            np.asarray(Function(Q).interpolate(y_expr).dat.data_ro),
        ),
        axis=1,
    )

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


def _shared_context(
    topology: Any, icesee_kwargs: Mapping[str, Any]
) -> _SharedPigContext:
    key = id(topology)
    context = _SHARED_CONTEXT_CACHE.get(key)
    if context is None:
        context = _build_shared_context(topology, icesee_kwargs)
        _SHARED_CONTEXT_CACHE[key] = context
    return context


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
    """Same depth-dependent control/warm forcing schedule used by
    ``_icepack_model.run_model``/``_icepack_enkf.generate_true_state``
    (1935-2017 PIG basal-melt scenario timeline).

    ``experiment`` defaults to ``EXPERIMENT_WRONG`` because every native
    member this adapter forecasts is an ensemble/forecast trajectory
    (mirroring ``_icepack_model.run_model``'s own default under the
    2026-09-28 reconciliation) -- the true trajectory is generated once,
    execution-mode-independently, by ``_icepack_enkf.generate_true_state``,
    never by this Mode-3 adapter."""

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
    """Whether this run's declared state composition includes
    ``basal_melt_field`` -- the generic, already-config-driven predicate
    modes 0-2 use (``config/_utility_imports.py``'s ``vec_inputs``), not
    the ``joint_estimation`` bookkeeping flag. ``vec_inputs`` reflects the
    packed vector's actual composition regardless of whether the current
    convention bookkeeps basal_melt_field as a state sub-block or a
    jointly-estimated parameter sub-block (see docs/execution-mode-3-
    design.md's 2026-09-28 reconciliation note)."""

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

    return {
        "Q": ctx.Q,
        "s0": ctx.s0,
        "beta0": ctx.beta0,
        "A0": ctx.A0,
        "bed": ctx.bed,
        "uThresh": icesee_kwargs["uThresh"],
        "GLThresh": icesee_kwargs["GLThresh"],
        "hThresh": icesee_kwargs["hThresh"],
        "water_to_ice": icesee_kwargs.get("water_to_ice", 1.0),
    }


# CORRECTED (2026-09-28, second pass): `joint_estimation` means "are we
# performing joint state-parameter estimation" -- it is NOT a proxy for
# "does basal_melt_field exist" (that question is answered by
# `_basal_melt_in_state`/`vec_inputs`, below). This nudge is genuinely
# joint-estimation-algorithm-specific -- it initializes the ensemble's
# WRONG prior belief about a jointly-estimated parameter, mirroring
# ``_icepack_enkf.initialize_ensemble`` in modes 0-2 -- so its
# `joint_estimation` gate is semantically correct and is kept. Under the
# current Idealized PIG config (joint_estimation=0, basal_melt_field is an
# ordinary state variable, true-vs-wrong divergence comes from
# BasalMeltRate's EXPERIMENT_TRUE/EXPERIMENT_WRONG selector instead), this
# gate is simply False and the nudge does not fire -- so it never operates
# simultaneously with the new experiment mechanism. If a future config
# re-enables joint_estimation for some other jointly-estimated parameter,
# this nudge correctly reactivates as part of that algorithm.
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
    """Perturb thickness and take the initial Icepack step for one member.

    Uses ``_member_initialization_context`` for the same member-keyed seed
    derivation (``initialization_seed``) modes 1/2 use, but draws the
    perturbation itself via ``coordinate_keyed_white_noise`` over this
    rank's OWNED physical coordinates only (``ctx.coords``, O(Nx/P_model)
    memory -- never the global mesh size). This makes the realization a
    deterministic function of (base_seed, member_id, physical coordinate),
    so the SAME mesh node gets the SAME N(0,1) draw regardless of which
    rank owns it or how many spatial ranks the mesh is partitioned across
    -- unlike a plain ``np.random.normal(0, 1, local_size)`` draw, whose
    value at a given array position is only meaningful for one fixed
    partition (confirmed: under real spatial decomposition, that draw
    assigned uncorrelated noise to different physical nodes depending on
    rank count, corrupting scientific R=1-vs-R>1 equivalence -- see
    docs/execution-mode-3-design.md's fix history / this session's Gate-3
    audit). NOT bit-identical to modes 1/2's own (still array-position-
    keyed, and itself not yet fixed for spatial decomposition) draw for
    the same member_id -- that old cross-mode parity was only ever
    incidentally true at R=1 for both, and was never a decomposition-
    invariant guarantee. The scientific formula (h0 + 15*N(0,1), one
    independent standard-normal draw per physical mesh node) is
    unchanged; only which physical node receives which realization is now
    well-defined.
    """

    ctx = _shared_context(topology, icesee_kwargs)

    with _member_initialization_context(icesee_kwargs, member_id) as _local_kwargs:
        noise = coordinate_keyed_white_noise(ctx.coords, _local_kwargs["seed"])
        h_perturbed = ctx.h0.dat.data_ro + 15 * noise

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
    if icesee_kwargs.get("joint_estimation"):
        basal_melt_field = _nudge_basal_melt(basal_melt_field, icesee_kwargs, ctx.Q)

    run_kwargs = _icepack_run_kwargs(ctx, icesee_kwargs)
    h, u, s, _floating, _grounded = Icepack(
        ctx.forward_solver,
        h_p,
        ctx.u0,
        ctx.smb,
        basal_melt_field,
        ctx.bed,
        float(icesee_kwargs["dt"]),
        ctx.h0,
        run_kwargs,
    )

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
) -> IcepackForecastFields:
    """Advance one member's persistent Icepack state by one native step.

    Mirrors ``_icepack_model.run_model``: the inflow boundary thickness is
    always the shared reference ``ctx.h0``, never the member's evolving
    thickness. Unlike an earlier version of this adapter, the basal-melt
    schedule's floating mask is now recomputed fresh from the member's
    CURRENT surface each step (mirroring ``_icepack_model.Icepack``'s own
    internal recomputation), not from a frozen initial ``ctx.floating0``
    -- floating/grounded are pure diagnostics of (surface, bed), never
    persistent DA state, so nothing needs to be packed, stored, or
    checkpointed to keep this correct (see docs/execution-mode-3-
    design.md's 2026-09-28 reconciliation note).
    """

    ctx = _shared_context(topology, icesee_kwargs)
    step = int(timestep)

    floating = _recompute_floating(ctx, state.surface)
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
        float(icesee_kwargs["dt"]),
        ctx.h0,
        run_kwargs,
    )

    return IcepackForecastFields(
        thickness=h,
        velocity=u,
        surface=s,
        basal_melt=basal_melt_field if state.basal_melt is not None else None,
    )


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


IDEALIZED_PIG_NATIVE_ADAPTER = IcepackNativeAdapter(
    initialize_member=initialize_member,
    forecast_member=forecast_member,
    reactivate_member=reactivate_member,
    allocate_member=allocate_member,
    layout_fingerprint=layout_fingerprint,
)
