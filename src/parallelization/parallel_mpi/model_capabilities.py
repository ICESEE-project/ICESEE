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
    """

    model_name: str
    supports_multi_rank_per_model: bool = False
    requires_one_rank_per_member: bool = False
    supports_distributed_state: bool = False
    default_ranks_per_model: Optional[int] = None
    supports_round_based_reuse: bool = True
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


# --- Conservative registrations for every application ICESEE currently
# ships. None has been run/verified with ranks_per_model > 1 in modes 0-2
# in this codebase's history, so every one is registered supports_multi_
# rank_per_model=False -- an honest "not yet," not a permanent "never."
register_model_capabilities(
    "lorenz",
    supports_multi_rank_per_model=False,
    supports_distributed_state=False,
    notes=(
        "3-variable ODE state has no useful spatial decomposition at all "
        "(see the mode-3 native adapter's identical reasoning, "
        "distributed_native_adapter.py) -- there is nothing for a second "
        "rank per member to own."
    ),
)
register_model_capabilities(
    "issm",
    supports_multi_rank_per_model=False,
    requires_one_rank_per_member=True,
    supports_distributed_state=False,
    notes=(
        "Python/MATLAB boundary is one persistent MatlabServer per active "
        "ICESEE rank (ISMIP_Choi); model_nprocs is ISSM's own internal "
        "MATLAB/PETSc parallelism, entirely separate from and orthogonal "
        "to ICESEE-level ranks_per_model (see "
        "_issm_native.py/mode3_runner.py's identical finding for mode 3's "
        "whole-member adapter). Future distributed-state support would need "
        "a genuinely new Python<->MATLAB boundary, not just a topology "
        "change -- confirmed unavailable in mode 3's own investigation "
        "(ISSM's internal PETSc solve gathers to one complete-member HDF5 "
        "write inside ISSM's own compiled core, invisible to and "
        "unmodifiable from ICESEE); supports_distributed_state stays False.\n"
        "Stage 4D.2 (mode 2, ISMIP_Choi): bounded persistent-server "
        "round-robin scheduling for Nens > (active servers), reusing the "
        "SAME generic mode-2 rounds/color/subcomm_size_min machinery "
        "Icepack's Stage 4C rounds regression already exercises (no "
        "separate scheduler written) -- new "
        "resource_plan.plan_resources(..., max_model_groups=...) parameter, "
        "wired to ICESEE's own ``issm_server_count`` config key. Confirmed "
        "with real MATLAB/ISSM runs (ISMIP_Choi, ~444x112/3347-vertex "
        "mesh): Nens=1 baseline, Nens=2/Nserver=1 (one persistent server "
        "correctly processes two logical members sequentially, no state "
        "leakage -- member states verified distinct and finite), Nens=4/"
        "Nserver=2 (2 rounds, all 4 members distinct/finite/correctly "
        "assigned), all clean (zero orphan MATLAB/mpiexec/issm.exe "
        "processes on exit). Nserver=1-vs-Nserver=2 scientific equivalence "
        "for identical Nens: MACHINE PRECISION (~1e-13 to 1e-16 relative, "
        "every state variable, both members) -- cleaner than Icepack's "
        "analogous Stage 4C result precisely because ISSM's state is "
        "always whole-member (never spatially decomposed across ICESEE "
        "ranks), so there is no decomposition-dependent nonlinear-solver-"
        "path noise for server count to introduce. requires_one_rank_per_"
        "member stays True in this registry despite this: it was never "
        "actually enforced anywhere in the codebase (confirmed by grep -- "
        "MatlabServer._check_conditions only ever rejected Nens < size, "
        "never blocked Nens > size), so it was already effectively dead/"
        "inert, and flipping capability-registry flags was out of Stage "
        "4D.2's explicit scope; a future pass should either start "
        "enforcing it accurately or retire it.\n"
        "Nested-MPI resource safety (Stage 4D.1's own predicted risk, "
        "confirmed with real evidence in Stage 4D.2): each active server's "
        "``mpiexec -np model_nprocs`` nested ISSM solve is launched "
        "independently, uncoordinated with every other active server's "
        "nested job on the same host, AND with icesee_da_full_parallel.py's "
        "own pre-existing adaptive model_nprocs scaling (observed inflating "
        "a configured model_nprocs=2 to 6-9 on an otherwise-idle host in "
        "real runs) -- report_nested_mpi_resource_usage (mat2py_utils.py, "
        "called from run_da_issm.py right after icesee_mpi_init) prints a "
        "non-fatal diagnostic (issm_server_count x configured model_nprocs "
        "vs os.cpu_count()) rather than blocking; real 18-nested-rank "
        "oversubscription on one machine completed correctly in testing, "
        "just with more contention, so a hard block was not justified.\n"
        "A real, reproducible failure-propagation bug was found and fixed "
        "this stage (not merely round-robin-specific, but directly "
        "relevant to it): generate_true_state/generate_nurged_state/"
        "initialize_ensemble/run_model_inverse/initialize_model's own "
        "MATLAB-command-send check in ISMIP_Choi's _issm_enkf.py/"
        "_issm_model.py caught a failed MATLAB command, killed the server, "
        "then silently returned None or fell through to read a stale/"
        "missing output file instead of propagating the failure -- "
        "reproduced directly (a real ISSM transient-solve error left a "
        "subsequent generate_nurged_state call hanging forever sending "
        "commands to the now-dead server). All five now re-raise, letting "
        "the shared mode-2 driver's own top-level handler reach "
        "comm_world.Abort(1) instead; verified both the pre-fix hang and "
        "the post-fix clean abort with real runs. run_model's own MAIN "
        "per-timestep forecast path already re-raised correctly before "
        "this stage (not touched). basal_friction_variation's "
        "_issm_enkf.py/_issm_model.py has the IDENTICAL bug pattern in six "
        "places, including its own main run_model forecast path (unlike "
        "ISMIP_Choi's, which was already safe there) -- confirmed by "
        "inspection, NOT fixed this stage (untested this session, out of "
        "scope); flag before relying on that application's failure "
        "handling."
    ),
)
register_model_capabilities(
    "icepack",
    supports_multi_rank_per_model=True,
    supports_distributed_state=True,
    notes=(
        "Three Stage 4C continuations established: Firedrake's mesh "
        "partition is genuinely distributed (not replicated) across a "
        "model group's ranks (supports_distributed_state=True, "
        "accurate). Distributed forecast, distributed EnKF analysis "
        "(uneven partitions, multiple analysis cycles), real spare "
        "ranks, and real ensemble rounds all work -- verified with real "
        "Icepack runs (P=2/4/6, Nens=1/2/4, ranks_per_model=2), fixing "
        "several real, confirmed, non-Icepack-specific bugs in the "
        "shared modes-0-2 machinery along the way (per-rank-local vs "
        "global state size passed to EnKF_fully_parallel_IO; "
        "icesee_get_index's local-nd-as-inter-variable-stride bug under "
        "uneven partitions; forecast write-back only persisting "
        "sub_rank==0's slice; an inference-plugin coordinate provider "
        "invoked as a single-rank Firedrake collective; two "
        "ranks_per_model==1-branch spare-rank COMM_NULL crashes; a "
        "distributed-ensemble-init branch that silently never "
        "initialized round-1+ members). Two genuinely Icepack-specific "
        "spare-rank gaps (run_da_icepack.py building a mesh "
        "unconditionally; prepare_random_field_coordinates reached by "
        "spares under random_field_method=graph) were also fixed.\n"
        "Stochastic decomposition invariance -- the one remaining gate "
        "-- was investigated and partially closed. FFT (the default "
        "random_field_method) defines its covariance over flat "
        "array-index distance ('periodic ring distance', per "
        "sample_periodic_exp_cov's own docstring) -- an intentional, "
        "documented, index-space model, not a bug; it is NOT physically "
        "decomposition-invariant by design and was deliberately left "
        "unchanged. The 'graph' method's adjacency/smoothing already "
        "used real physical coordinates, but its initial white-noise "
        "draw (`np.random.randn(n)`) was array-index-ordered -- a real "
        "implementation gap relative to graph's own coordinate-based "
        "intent. Fixed with a narrowly scoped, authorized change: "
        "coordinate_keyed_white_noise "
        "(src/run_model_da/_error_generation.py) keys the initial draw "
        "to (already member/variable/init-vs-process-noise-namespaced "
        "seed, quantized physical coordinate) instead of array "
        "position -- graph construction, adjacency, smoothing, "
        "correlation-length interpretation, sig_Q semantics, and FFT "
        "are all untouched. Verified exactly permutation-invariant in "
        "isolation (machine precision), and verified against real "
        "Icepack R=1-vs-R=2 coordinate-mapped noise dumps: correlation "
        "improved from ~0.0-0.4 (pre-fix, genuinely different physical "
        "fields) to ~0.99998-1.0 (post-fix). A small residual "
        "(~0.1-1.5% relative, ~1e-1 absolute) remains, root-caused to "
        "k-NN tie-breaking in the graph's adjacency construction for "
        "this test mesh's perfectly regular grid (exactly-equidistant "
        "neighbor candidates resolved differently depending on the "
        "order cKDTree receives coordinates in) -- confirmed directly "
        "(1 of 81 nodes' neighbor set differed after coordinate "
        "mapping) and out of this fix's authorized scope (k-NN "
        "definition was explicitly not to be changed); expected to be "
        "smaller or absent on a less perfectly symmetric production "
        "mesh. This is orders of magnitude smaller than the pre-fix "
        "discrepancy and does not itself block declaring the noise "
        "generator decomposition-invariant.\n"
        "However, a full real-Icepack R=1-vs-R=2 DA-trajectory "
        "comparison (coordinate-mapped, graph method, same seed) still "
        "shows large (100-500% relative) discrepancies for h/u/v -- "
        "traced to a SEPARATE, Icepack-application-specific, "
        "out-of-scope bug: _icepack_enkf.py's initialize_ensemble "
        "constructs the deterministic h_nurge_ic/u_nurge_ic initial-"
        "condition bump via LOCAL ARRAY-INDEX slicing "
        "(h0.dat.data_ro[:h_indx], where hdim = h0.dat.data_ro.size is "
        "this rank's own local DOF count) rather than physical-"
        "coordinate selection -- decomposition-dependent for the same "
        "reason the graph noise was, but in deterministic model-adapter "
        "code, not the shared stochastic-field module this pass was "
        "authorized to touch. (smb, which the nurge does not touch, "
        "matched to ~0.6-0.7% relative -- consistent with the graph fix "
        "alone working correctly; h/u/v, which the nurge does touch, did "
        "not.) supports_multi_rank_per_model stays False: full noisy-DA "
        "R=1-vs-R=2 scientific equivalence is blocked by this newly "
        "identified, separate, Icepack-specific nurge-construction bug, "
        "not by the stochastic-field infrastructure, which this pass "
        "brought to a decomposition-invariant (graph) or honestly-"
        "documented-as-not (FFT) state.\n"
        "A fourth continuation closed the graph method's one remaining "
        "gap: its k-NN adjacency construction broke exact-distance ties "
        "non-deterministically (cKDTree's own tie order depends on "
        "input array order) -- reproduced exactly (1 of 81 real Icepack "
        "mesh nodes affected) and fixed narrowly in _build_knn_adjacency "
        "(src/run_model_da/_error_generation.py): query a small margin "
        "of extra candidates beyond k, then select the final k by "
        "(distance, canonical quantized-coordinate key) instead of "
        "whatever order cKDTree returned them in. k, the distance "
        "metric, weighting, smoothing, correlation-length interpretation, "
        "and normalization are all unchanged. Verified: zero neighbor-set "
        "mismatches (was 1/81) on the real Icepack mesh, and the R=1-vs-"
        "R=2 graph noise comparison (both init and process noise, all "
        "four variables) improved from ~0.9999985-1.0 correlation / "
        "~0.1-1.5% residual to **exact machine precision** (correlation "
        "1.0 to 1e-16, max-abs-diff ~1e-14-1e-16). The graph stochastic-"
        "field infrastructure is therefore now fully decomposition-"
        "invariant, not merely 'almost' -- confirmed, not assumed.\n"
        "This isolates the full-DA discrepancy (found in the third "
        "continuation) entirely to the nurge-construction bug, with "
        "clean corroborating evidence: re-running the same real R=1-vs-"
        "R=2 full-DA comparison with the tie-break fix in place left "
        "smb (never touched by the nurge) at the same small ~0.6% "
        "residual as before (PETSc/solver-accumulation noise, unrelated "
        "to graph), while h/u/v (touched by the nurge) were completely "
        "unchanged at their prior ~300-500% relative discrepancy -- "
        "proving the graph fix is unrelated to and does not mask the "
        "nurge issue. That issue's exact recovered semantics: "
        "_icepack_enkf.py's initialize_ensemble unpacks 'x' and 'Lx' "
        "(physical x-coordinate and domain length) but never uses "
        "either in its actual bump construction, which instead slices "
        "h0.dat.data_ro[:h_indx] -- this rank's own local-array-order "
        "DOFs. Empirically, this does NOT correspond to a clean "
        "physical rule even at R=1 (the array-slice's x-range was "
        "0-3125 out of a full 0-5000 domain, only 59% overlapping with "
        "the actual lowest-x nodes) -- so 'preserve exact legacy R=1 "
        "numerics' and 'implement the evidently-intended physical "
        "threshold' are two different, mutually exclusive outcomes, not "
        "resolved from the code/config alone (nurged_entries_percentage "
        "'percentage of nurged entries in the state vector' can be read "
        "either way). Per this task's own instruction to stop rather "
        "than guess at genuinely ambiguous scientific intent, this was "
        "investigated and reported, not silently fixed -- "
        "supports_multi_rank_per_model stayed False pending that "
        "decision (see the accompanying report's alternatives A/B/C).\n"
        "A fifth continuation resolved that decision: Option A (physical "
        "x-coordinate-based nudging), explicitly authorized over "
        "B/C. _physical_nudge_expr (in both the actually-executed "
        "synthetic_ice_stream/_icepack_enkf.py and the shared fallback "
        "icepack_utils/_icepack_enkf.py -- call_model() prefers the "
        "former but falls back to the latter, so both needed the fix) "
        "replaces the array-index slice with a UFL conditional over "
        "physical x: xi = x/Lx, taper = -amplitude*(1 - xi/threshold) "
        "for xi <= threshold else 0, applied per-rank from each rank's "
        "own local mesh coordinates (no gather needed -- interpolate() "
        "is already a purely local operation per Firedrake DOF). "
        "nurged_entries_percentage keeps its config key name but now "
        "means 'fraction of physical x-domain length', not 'fraction of "
        "flattened state-vector index range' -- an intentional, "
        "authorized correction, not a legacy-preserving change (the old "
        "R=1 behavior was itself accidental/decomposition-dependent, "
        "confirmed above). Only h/u/v are touched (v via the same "
        "u-taper, since u0.dat.data_ro[:,1] is the v-component); smb is "
        "untouched, verified bit-for-bit identical with the nudge "
        "toggled on/off (test_icepack_physical_nudge.py::"
        "test_smb_unaffected_by_nudge_amplitude).\n"
        "Full validation (test_icepack_physical_nudge.py, real Icepack/"
        "Firedrake, nx=4,ny=4): (1) taper formula itself -- exact linear "
        "ramp, correctly masked, monotonic, collapses to zero for "
        "amplitude=0 or threshold<=0. (2) Deterministic (sig_Q=0) R=1-vs-"
        "R=2: machine precision (<1e-8 relative, per-variable, as a set "
        "-- R=1 and R=2 assign the same physical node to different flat "
        "array positions, a benign Firedrake-partitioner DOF-ordering "
        "artifact confirmed separately and unrelated to this fix; valid "
        "comparison is per-variable-block-as-a-set, not raw flat index) "
        "through the FULL forecast+analysis cycle, not just t=0. (3) "
        "Noisy (graph method, default nonzero sig_Q) initial ensemble "
        "mean R=1-vs-R2: also machine precision (<1e-8 relative) -- the "
        "coordinate-keyed noise and the physical nudge compose cleanly. "
        "(4) Full DA cycle (2 analysis updates): relative residual grows "
        "from machine precision at t=0 to ~1e-5-1e-3 (h/u) and up to "
        "~2-3% (v, smb -- small-magnitude fields where PETSc/Firedrake's "
        "nonlinear diagnostic-solve convergence path differs slightly "
        "under different decompositions) by the final step -- expected "
        "solver-path noise, NOT a decomposition-invariance defect: "
        "confirmed by the deterministic-only case (2) staying at machine "
        "precision through the same number of analysis cycles, isolating "
        "the residual's source to the nonlinear solve's interaction with "
        "the (now decomposition-invariant) stochastic perturbation, not "
        "to the nudge or the DA algorithm itself. (5) Real spare-rank "
        "(P=6,Nens=2,ranks_per_model=2) and rounds (P=4,Nens=4,"
        "ranks_per_model=2) regressions both pass with graph method, the "
        "physical nudge, and inference_plugin_enabled=true (the "
        "production default, not disabled) -- see "
        "test_icepack_multirank_analysis.py.\n"
        "supports_multi_rank_per_model is now True: the nurge-"
        "construction bug that was the sole blocker is fixed and "
        "verified (deterministic + noisy-init at machine precision, full "
        "DA cycle within expected solver-noise tolerance, spare ranks "
        "and rounds both real-Icepack-verified with the production "
        "inference plugin active). The one remaining caveat is the "
        "pre-existing, orthogonal FFT-vs-graph gap documented above: FFT "
        "is not decomposition-invariant by design (index-space "
        "covariance), so validate_stochastic_method_for_decomposition "
        "(this module) rejects ranks_per_model>1 combined with an active "
        "stochastic field (nonzero sig_Q) and random_field_method='fft' "
        "at ICESEE-mpi-init time, with a message recommending "
        "random_field_method='graph' (the decomposition-invariant "
        "method validated above) or, for a user who explicitly wants "
        "FFT's index-space semantics under multi-rank anyway, "
        "acknowledge_index_space_stochastic_multirank=True. This guard "
        "does NOT fire for a deterministic (sig_Q all-zero) run even "
        "though fft is params.yaml's configured default method -- "
        "verified directly (test_model_capabilities.py) -- so it never "
        "blocks a legitimate deterministic multi-rank run. It also only "
        "fires for a model with supports_distributed_state=True (a real "
        "spatial DOF partition, currently only icepack): a first version "
        "of this guard fired for every model with ranks_per_model>1, "
        "which broke Lorenz96's already-passing real-MPI P>Nens "
        "oversubscribed-resource-planning regression suite (Lorenz's "
        "small replicated ODE state has no spatial DOF ordering for "
        "FFT's index-space covariance to disagree with) -- caught by "
        "running the full regression suite after this change and fixed "
        "by gating on supports_distributed_state."
    ),
)
register_model_capabilities(
    "flowline",
    supports_multi_rank_per_model=False,
    supports_distributed_state=False,
    notes=(
        "State is a small (2*NX+1) in-process scipy/JAX solve with no "
        "spatial MPI decomposition of any kind, same category as Lorenz96."
    ),
)
