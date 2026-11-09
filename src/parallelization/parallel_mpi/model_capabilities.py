# ==============================================================================
# @des: Per-application opt-in capability metadata for Stage 4's hierarchical
#       resource planner (resource_plan.py).
#
#       resource_plan.plan_resources() is deliberately model-agnostic: given
#       (world_size, nens, ranks_per_model) it computes a valid topology for
#       ANY of those inputs, including ranks_per_model > 1, without knowing
#       or caring which application asked for it. Whether ranks_per_model
#       > 1 is actually *meaningful* for a given application (does its
#       forecast step know what to do with more than one rank per member?)
#       is a separate question this module answers, following the same
#       registry pattern already established for Mode 3
#       (src/parallelization/distributed_mode3_registry.py) -- a sibling
#       concern, not the same one: that registry gates whether
#       execution_mode 3 exists at all for a model; this one gates whether
#       a *requested* ranks_per_model > 1 is safe under modes 0-2's
#       existing hierarchical topology. Do not conflate the two, and do
#       not force Mode 3 through this registry or this one through Mode
#       3's.
#
#       Every current application is registered here with conservative,
#       verified-honest capabilities (none has been tested with
#       ranks_per_model > 1 in modes 0-2 yet -- see each registration's
#       notes) rather than left unregistered, so a request for
#       ranks_per_model > 1 fails with a clear, named error instead of
#       either a generic KeyError or -- worse -- silently running with
#       multiple ranks a model's forecast step does not actually use.
# ==============================================================================
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class ModelCapabilities:
    """One application's declared hierarchical-topology capabilities.

    ``supports_multi_rank_per_model``: whether this application's forecast
    step is verified to do something meaningful with more than one rank
    cooperating on a single ensemble member (e.g. a spatially partitioned
    mesh). ``False`` does not mean "can never support it" -- it means
    "not yet verified," and ``ranks_per_model > 1`` is refused rather than
    silently accepted for such a model.

    ``requires_one_rank_per_member``: stronger than merely not supporting
    multi-rank -- this application's integration additionally requires
    ``ranks_per_model == 1`` even when ``ranks_per_model`` is left at its
    default (e.g. ISSM's persistent per-rank ``MatlabServer`` identity,
    ``ens_id = rank``, that its own scheduling depends on). Currently
    informational; ``plan_resources``'s own legacy-default policy already
    produces ``ranks_per_model == 1`` for every registered model's actual
    shipped configuration (all run with ``Nens >= world_size`` today), so
    no enforcement branch has been exercised by a real conflict yet.

    ``supports_distributed_state``: whether this application can expose
    genuinely distributed (non-replicated) state across the ranks in one
    model group -- orthogonal to ``StateOwnership`` (see module docstring
    in ``src/utils/state_ownership.py``): this field describes what the
    *model* can do; ``StateOwnership`` describes what a *given run*
    actually does with one rank's data. Do not infer one from the other.

    ``default_ranks_per_model``: an application-preferred default, if any.
    ``None`` means "no preference beyond the planner's own legacy
    default" (every current application).

    ``supports_round_based_reuse``: whether this application's per-member
    state can be safely reused/reinitialized across scheduling rounds
    (mode 2's round loop already assumes this for every current
    application: ``ens_id = color + round_id * subcomm_size_min``).

    ``requires_node_local_model_groups``: whether every rank of one model
    group must be on the same compute node. Set by an application whose
    multi-rank model setup performs collective parallel file I/O that has
    been observed to deadlock when a model group's communicator spans
    nodes. A launch that places one group across nodes is then refused at
    startup (see ``validate_node_local_model_groups``) instead of hanging.
    """

    model_name: str
    supports_multi_rank_per_model: bool = False
    requires_one_rank_per_member: bool = False
    supports_distributed_state: bool = False
    default_ranks_per_model: Optional[int] = None
    supports_round_based_reuse: bool = True
    requires_node_local_model_groups: bool = False
    notes: str = ""


_REGISTRY: dict[str, ModelCapabilities] = {}


def register_model_capabilities(
    model_name: str,
    *,
    supports_multi_rank_per_model: bool = False,
    requires_one_rank_per_member: bool = False,
    supports_distributed_state: bool = False,
    default_ranks_per_model: Optional[int] = None,
    supports_round_based_reuse: bool = True,
    requires_node_local_model_groups: bool = False,
    notes: str = "",
) -> None:
    """Register ``model_name``'s hierarchical-topology capabilities.

    Re-registering the same ``model_name`` overwrites the previous entry
    (mirrors ``distributed_mode3_registry.register_execution_mode_3``'s
    same choice, for the same reason: safe re-import under pytest).
    """
    if not model_name:
        raise ValueError("model_name must be a non-empty string")
    _REGISTRY[model_name] = ModelCapabilities(
        model_name=model_name,
        supports_multi_rank_per_model=bool(supports_multi_rank_per_model),
        requires_one_rank_per_member=bool(requires_one_rank_per_member),
        supports_distributed_state=bool(supports_distributed_state),
        default_ranks_per_model=default_ranks_per_model,
        supports_round_based_reuse=bool(supports_round_based_reuse),
        requires_node_local_model_groups=bool(requires_node_local_model_groups),
        notes=notes,
    )


# Conservative default for any unregistered model: no multi-rank support.
# Never silently treat "unknown" as "supported."
_UNREGISTERED_DEFAULT = ModelCapabilities(model_name="<unregistered>")


def get_model_capabilities(model_name: Optional[str]) -> ModelCapabilities:
    """Return ``model_name``'s registered capabilities, or the
    conservative unregistered default (no multi-rank support) if it has
    not registered any."""
    if not model_name:
        return _UNREGISTERED_DEFAULT
    return _REGISTRY.get(model_name, _UNREGISTERED_DEFAULT)


def validate_ranks_per_model_request(model_name: Optional[str], requested_ranks_per_model) -> None:
    """Fail early and clearly if ``requested_ranks_per_model`` asks for
    more than one rank per model group on a model that has not verified
    it can use them.

    Deliberately narrow: only rejects an *explicit* request for more than
    one rank (an integer > 1, or the literal string ``"auto"``, which the
    caller is responsible for resolving and re-checking if it lands above
    1 for a world_size > nens topology). ``None`` -- the default every
    currently shipped configuration uses -- is never rejected here; it is
    resolved by ``plan_resources``'s own legacy policy, which already
    yields ``ranks_per_model == 1`` for every registered model's actual
    configuration (all run with ``Nens >= world_size``). This keeps
    today's baseline behavior completely unaffected while still refusing
    a genuinely new multi-rank request for a model that has not earned
    it, per Stage 4's "fail early... do not silently pretend multiple
    ranks are useful" requirement.
    """
    if requested_ranks_per_model is None:
        return
    if isinstance(requested_ranks_per_model, str):
        if requested_ranks_per_model.strip().lower() != "auto":
            return  # let plan_resources raise its own clear error for garbage input
        # "auto" is checked by the caller after resolution (it may resolve
        # to 1 for a world_size <= nens topology, which needs no capability
        # at all); this function only judges genuinely explicit requests.
        return
    try:
        requested = int(requested_ranks_per_model)
    except (TypeError, ValueError):
        return  # let plan_resources raise its own clear error for garbage input
    if requested <= 1:
        return

    capabilities = get_model_capabilities(model_name)
    if not capabilities.supports_multi_rank_per_model:
        raise ValueError(
            f"execution_mode 2 was asked for ranks_per_model={requested} for "
            f"model {model_name!r}, but that model has not registered "
            "supports_multi_rank_per_model=True in "
            "src/parallelization/parallel_mpi/model_capabilities.py. "
            "Its forecast step is not verified to use more than one rank "
            "per ensemble member. Use ranks_per_model=1 (or omit it), or "
            "verify and register multi-rank support for this model first."
        )


def node_spanning_model_groups(hosts, group_of_rank):
    """Return ``{group: sorted hosts}`` for every model group whose ranks
    sit on more than one host.

    ``hosts[r]`` is world rank ``r``'s host name and ``group_of_rank(r)``
    its model-group id (``None`` for a rank in no group, e.g. a spare).
    Pure function: every rank evaluates it on the same allgathered input,
    so every rank reaches the same decision.
    """
    group_hosts: dict = {}
    for rank, host in enumerate(hosts):
        group = group_of_rank(rank)
        if group is None:
            continue
        group_hosts.setdefault(group, set()).add(host)
    return {
        group: sorted(names)
        for group, names in sorted(group_hosts.items())
        if len(names) > 1
    }


def validate_node_local_model_groups(model_name, world, group_of_rank, ranks_per_model):
    """Refuse a launch whose model groups span compute nodes when
    ``model_name`` requires node-local groups; otherwise warn on rank 0.

    Collective over ``world`` (one ``allgather`` of host names); must be
    called by every rank of ``world``. Raises the same error on every rank,
    so the job exits cleanly instead of deadlocking later inside a model
    group's collective I/O. Topology is never changed here.
    """
    import socket

    if int(ranks_per_model or 1) <= 1:
        return {}
    hosts = world.allgather(socket.gethostname())
    spanning = node_spanning_model_groups(hosts, group_of_rank)
    if not spanning:
        return spanning
    shown = ", ".join(
        f"group {group} on {'+'.join(names)}" for group, names in list(spanning.items())[:4]
    )
    if len(spanning) > 4:
        shown += f", ... ({len(spanning)} groups in total)"
    message = (
        f"{len(spanning)} model group(s) of {ranks_per_model} ranks span more than "
        f"one compute node ({shown}). Launch so that each node holds a whole "
        f"number of model groups: e.g. srun -N <nodes> --ntasks-per-node=<a "
        f"multiple of {ranks_per_model}> with -n equal to nodes x tasks-per-node."
    )
    if get_model_capabilities(model_name).requires_node_local_model_groups:
        raise RuntimeError(
            f"[ICESEE] model {model_name!r} requires node-local model groups "
            "(its collective model-setup I/O can deadlock across nodes): " + message
        )
    if world.Get_rank() == 0:
        print(f"[ICESEE] WARNING: {message}", flush=True)
    return spanning


def validate_stochastic_method_for_decomposition(
    model_name, resolved_ranks_per_model, icesee_kwargs
):
    """Reject the scientifically ambiguous combination of a genuinely
    distributed model group (``ranks_per_model > 1``) with stochastic
    random fields active under ``random_field_method == "fft"``.

    This is deliberately a *separate* concern from
    ``validate_ranks_per_model_request`` above: that function gates
    whether a model's *forecast step* is verified safe to run distributed
    at all; this one gates whether a *stochastic-field method* is safe to
    combine with a distributed model group once model-level support is
    otherwise fine. FFT (src/run_model_da/_error_generation.py's
    ``sample_periodic_exp_cov``) defines its covariance over flat state-
    vector array-index distance, not physical coordinates -- an
    intentional, documented design (see that function's own docstring:
    "periodic ring distance"), not a bug, and this pass does not change
    it. But Firedrake's own DOF ordering is decomposition-dependent, so
    FFT's realization is not guaranteed independent of
    ``ranks_per_model`` for a genuinely distributed model -- silently
    allowing it there would let a user believe their stochastic ensemble
    is physically reproducible across decompositions when it is not.
    ``random_field_method == "graph"`` (also in ``_error_generation.py``)
    is coordinate-based and, as of Stage 4C's graph fixes
    (coordinate-keyed white noise + deterministic k-NN tie-break), is
    decomposition-invariant to machine precision -- the safe choice here.

    Only fires when stochastic noise is actually active (any nonzero
    ``sig_Q`` entry): a deterministic, no-random-field distributed run
    must not be rejected merely because FFT is the configured default
    (Stage 4C explicitly requires this). An escape hatch,
    ``icesee_kwargs["acknowledge_index_space_stochastic_multirank"]``
    (default ``False``), lets a caller explicitly accept index-space
    semantics for a specific run rather than being silently switched to
    graph -- this function never changes the requested method itself.
    """
    try:
        rpm = int(resolved_ranks_per_model)
    except (TypeError, ValueError):
        return
    if rpm <= 1:
        return

    # The FFT-vs-graph index-space/physical-coordinate ambiguity this
    # guards against only exists for a model whose state is genuinely
    # spatially distributed across a model group's ranks (a real mesh
    # partition whose DOF order is decomposition-dependent, e.g. a
    # Firedrake Function -- see get_mesh_coordinates/prepare_random_field_
    # coordinates). A model registered with supports_distributed_state=
    # False (Lorenz96's replicated small ODE state, ISSM's whole-state
    # MATLAB transfer, flowline's non-distributed state) has no such
    # ambiguity even if ranks_per_model > 1 for it (e.g. Lorenz96's
    # oversubscribed P>Nens resource-planning cases): there is no
    # per-rank spatial DOF slice for FFT's index-space covariance to
    # disagree with. Reject only for a genuinely spatially distributed
    # model -- currently just icepack.
    if not get_model_capabilities(model_name).supports_distributed_state:
        return

    sig_q = icesee_kwargs.get("sig_Q", None)
    if sig_q is None:
        sig_q = []
    try:
        stochastic_active = any(float(s) != 0.0 for s in sig_q)
    except TypeError:
        # Not iterable (e.g. a bare scalar) -- fall back to a direct
        # truthiness/nonzero check rather than `sig_q or []`, which raises
        # ValueError for a numpy array with more than one element.
        try:
            stochastic_active = float(sig_q) != 0.0
        except (TypeError, ValueError):
            stochastic_active = bool(sig_q)
    if not stochastic_active:
        return

    method = str(icesee_kwargs.get("random_field_method", "fft")).strip().lower()
    if method != "fft":
        return
    if bool(icesee_kwargs.get("acknowledge_index_space_stochastic_multirank", False)):
        return

    raise ValueError(
        f"ranks_per_model={rpm} for model {model_name!r} with stochastic "
        "random fields active (a nonzero sig_Q entry) and "
        "random_field_method='fft' is a scientifically ambiguous "
        "combination: FFT defines its covariance over flat state-vector "
        "array-index distance, not physical coordinates, so its "
        "realization is not guaranteed independent of Firedrake's own "
        "decomposition-dependent DOF ordering under ranks_per_model > 1 "
        "(see this model's own notes in model_capabilities.py for the "
        "full investigation). Set random_field_method='graph' for a "
        "decomposition-invariant physical stochastic field, or pass "
        "acknowledge_index_space_stochastic_multirank=True to explicitly "
        "accept index-space (topology-dependent) stochastic semantics "
        "for this run."
    )


# Register capabilities for each adapter; examples must validate their own layouts.
register_model_capabilities(
    "lorenz",
    supports_multi_rank_per_model=False,
    supports_distributed_state=False,
    notes=(
        """Replicated ODE state has no spatial MPI decomposition."""
    ),
)
register_model_capabilities(
    "issm",
    supports_multi_rank_per_model=False,
    requires_one_rank_per_member=True,
    supports_distributed_state=False,
    notes=(
        """ISMIP_Choi uses one persistent MatlabServer per active ICESEE rank.
ISSM model_nprocs controls nested MATLAB/PETSc processes, independently of
ICESEE ranks_per_model. issm_server_count bounds active servers; generic
rounds schedule additional members. Budget all nested processes explicitly.
The whole-member Python/MATLAB interface does not expose distributed owned
state. ISMIP_Choi propagates MATLAB command failures to the MPI driver;
basal_friction_variation failure handling requires separate validation."""
    ),
)
register_model_capabilities(
    "icepack",
    supports_multi_rank_per_model=True,
    supports_distributed_state=True,
    # Firedrake's collective CheckpointFile open (parallel HDF5 -> Open
    # MPI's OMPIO, initialization file on NFS) can deadlock inside
    # MPI_File_open when the model group's communicator spans compute
    # nodes: the shared-file-pointer component is selected inconsistently
    # across the nodes, so one rank waits in a barrier inside
    # mca_sharedfp_lockedfile_file_open while its partner has moved on to
    # an allreduce in mca_common_ompio_set_view (PACE, Open MPI 5.0.9,
    # 2026-10-03, gdb on a 38-rank Nens=19 model_nprocs=2 run; intermittent:
    # isolated probes of the same open often succeed). The whole run then
    # waits forever at the next world collective. Node-local groups never
    # hung across all runs.
    requires_node_local_model_groups=True,
    notes=(
        """Firedrake fields are partitioned across each model group's spatial ranks.
Global packed sizes and owned-row mappings must exclude halo duplicates.
Physical-coordinate nudging and coordinate-keyed graph noise avoid using
partition-local array order as a physical location. Graph neighbor ties use
canonical coordinate keys.
FFT covariance is defined in array-index space and can change with mesh
numbering. Multi-rank stochastic runs therefore require graph noise or explicit
acknowledgment of FFT index-space semantics; deterministic runs are exempt.
Nonlinear solver paths and MPI reduction order can still affect trajectories.
Validate each example and topology rather than assuming bitwise equivalence."""
    ),
)
register_model_capabilities(
    "flowline",
    supports_multi_rank_per_model=False,
    supports_distributed_state=False,
    notes=(
        """The in-process scipy/JAX solve uses whole-member state, not spatial MPI decomposition."""
    ),
)
