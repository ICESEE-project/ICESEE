# 0005 — State ownership: `replicated` vs `distributed`

## Status

Accepted. Implemented in `src/utils/state_ownership.py`
(`StateOwnership`, `resolve_state_ownership`, `combine_member_state`,
`ensure_state_array`), consumed by
`src/parallelization/_mpi_generate_true_wrong_state.py`,
`src/parallelization/_mpi_ensemble_intialization.py`, and
`src/utils/tools.py` (`icesee_get_index`). Validated end-to-end (real
`mpirun`, real driver, not only synthetic tests) for the `replicated` case
via Lorenz-96 across `P ∈ {1,2,3,4,8,10,16}`, `Nens=4`
(`src/tests/test_lorenz96_mode2_p_nens_matrix.py`). The `distributed` case
is validated statically only (Icepack's `nd` is confirmed a genuine
per-rank Firedrake/PETSc dof count by direct code inspection) — no real
multi-rank Icepack run has been executed. This ADR should not be read as
claiming that validation; keep the language precise until it is done.

## Context

A model communicator (the ranks cooperating on one ensemble member) can
relate to that member's state in two fundamentally different ways, and
ICESEE's shared code was, until this ADR, silently assuming only one of
them:

**`replicated`**: every rank in the model communicator computes/holds an
identical, complete copy of the state. Lorenz-96 is the concrete example —
`nd` (its reported state size) is a fixed configuration constant
(`num_state_vars`), the same value on every rank regardless of how many
ranks share that member, because the model does no spatial decomposition
of its own across those ranks.

**`distributed`**: each rank owns a disjoint, contiguous partition of the
state. Icepack is the concrete example — `nd` is `h0.dat.data.size`, a
genuine Firedrake-owned (non-ghost) local dof count that differs by rank
whenever the member's mesh is partitioned across more than one rank.

ICESEE's pre-existing shared code (`global_shape = sum(dim_list)`,
`subcomm.gather(...) + vstack/hstack`, and `icesee_get_index`'s
rank-indexed, cumulative-offset `dim_list` handling) was written only for
the `distributed` case. It happens to be correct for Icepack, by
coincidence of what Icepack's adapter puts in `nd`, not because the code
understood the distinction. For Lorenz-96 at `ranks_per_model > 1`, this
silently produced a `global_shape` inflated by the model-communicator size
(summing an already-full, identical value from every rank) and, separately,
concatenated identical per-rank copies instead of recognizing them as
redundant replicas — confirmed via a live reproduction (`ValueError: could
not broadcast input array from shape (3,) into shape (6,)`) once a
multi-rank-per-member configuration was actually exercised (it had not
been, for any application, before this investigation).

## Decision

Introduce an explicit `StateOwnership` descriptor
(`distribution: "replicated" | "distributed"`, `global_size`, `local_size`,
`local_offset`) and two operations keyed on it:

- `resolve_state_ownership(icesee_kwargs, subcomm, local_size=None)` —
  determines global/local size and offset from a model communicator's
  reported per-rank size(s).
- `combine_member_state(ownership, subcomm, data, root=0)` — assembles one
  member's global state from per-rank contributions: a no-op (no
  communication) for `replicated`, a gather+concatenate for `distributed`.

The **default is `distributed`**, not `replicated`. This is a deliberate
compatibility choice, not an arbitrary one: every model adapter shipped
before this ADR already relies on the pre-existing sum/gather arithmetic
being correct, and it *is* correct for every adapter whose `nd` is already
a genuine local partition (Icepack; ISSM trivially, since it never has
more than one rank per member — see ADR 0004/Stage 3 ISSM findings). A
model adapter opts into `replicated` explicitly, in its own configuration
(`state_distribution: replicated` in `applications/lorenz_model/examples/
lorenz96/params.yaml`), rather than the shared default changing under
adapters that were already correct. Lorenz-96 is the one adapter that
needed this: at `ranks_per_model > 1` it was never correct before (it
deadlocked/crashed), so opting in changes nothing that previously worked.

A related, generically-scoped normalization, `ensure_state_array`
(`np.atleast_1d`), addresses a separate but adjacent issue: a model
adapter can return a bare NumPy scalar for a single-element variable block
(`hdim == 1`) rather than a one-element array — the same logical quantity,
different NumPy representation. This is not an ownership/distribution
concept (a `distributed` partition is already guaranteed array-shaped by
construction); it is applied once, uniformly, wherever per-variable
adapter output is consumed, rather than by checking the calling model's
name.

## Consequences

- `icesee_get_index` (`src/utils/tools.py`) must treat `replicated` state
  identically across every rank (index as if `comm.Get_rank() == 0`
  unconditionally) rather than partition-index it — implemented alongside
  this ADR, since the two had to be consistent for either to be correct.
- Any new model adapter must declare `state_distribution` explicitly if it
  is `replicated`; silence means `distributed`, and an adapter whose `nd`
  is not a genuine local partition will reproduce the original Lorenz-96
  failure mode if run with `ranks_per_model > 1` without declaring it.
- This ADR does not implement `ranks_per_model` or block-contiguous rank
  placement (ADR 0003, still Proposed) — it only makes the state-size
  arithmetic those features will depend on correct for both ownership
  kinds. ADR 0003 should build on this, not duplicate it.
- The ISSM shared-infrastructure leak identified during this work
  (`EnKF_parallel_io.py`'s `_requires_shared_finalization` hardcoding
  `model_name == "issm"`) is unrelated to state ownership specifically and
  is tracked separately (ADR 0002 / model-adapter boundary), not resolved
  here.
