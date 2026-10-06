# 0004 — Distributed, rank-local state ownership (no global materialization)

## Status

Proposed for modes 0-2's forecast/model-adapter boundary. Already designed
in detail for a fully spatially distributed execution mode in
`docs/execution-mode-3-design.md` ("Status and purpose": mode 3 is a
proposed architecture, not yet accepted by the runtime configuration).
This ADR is the short decision record; `docs/execution-mode-3-design.md`
remains the authoritative detailed design for the mode-3 case.

## Context

Mode 2's *analysis* step already avoids materializing a global `Nx x Ne`
matrix (ADR 0001). The *forecast*/model layer has not been held to the same
standard: `idealized_pig`'s model code still moves whole per-member state
vectors, and no current execution mode partitions a single member's spatial
state across multiple ranks in a way ICESEE's shared code understands
directly (today that's entirely Firedrake's/ISSM's own internal domain
decomposition, invisible to ICESEE).

## Decision (proposed)

Extend the "local ownership, no global materialization" principle already
proven in ADR 0001 to the model-adapter boundary (ADR 0002): the adapter
exposes only a rank's locally owned state slab plus enough layout metadata
(global offsets/ownership ranges) for ICESEE's analysis code to address it,
never a global numpy array. Mode 3's two-dimensional (ensemble x spatial)
topology in `docs/execution-mode-3-design.md` is the long-term target for
models whose per-member state itself must be spatially distributed; this
ADR governs the same principle applied incrementally to modes 0-2 first.

## Consequences

- No global state array should be introduced anywhere in new shared code
  as a "temporary convenience" — if a global view is genuinely needed
  (e.g. for a diagnostic), it must be an explicit, clearly bounded, opt-in
  operation, not the default data path.
- This directly supersedes the (already-dead) `O(P*Nx*Ne)` gather/broadcast
  pattern in `src/parallelization/_parallel_i_o.py`
  (`gather_and_broadcast_data_default_run`) as a design direction — that
  code should be removed rather than reused as a starting point.
- Applies equally to modes 0/1 once/if their analysis path is migrated off
  the legacy `Nx x Nx`-scale `EnKF.py` — see the roadmap discussion in the
  audit for sequencing.
