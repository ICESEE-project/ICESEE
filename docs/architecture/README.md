# Architecture Decision Records

This directory records *why* major ICESEE architectural decisions were made,
not just what the code currently does (that belongs in `docs/*.md` and the
Wiki). Each ADR is numbered, short, and has one of these statuses:

- **Proposed** — a documented direction, not yet fully implemented or agreed.
- **Accepted** — implemented and in active use; changing it should start a
  new ADR that supersedes this one, not a silent rewrite.
- **Superseded by NNNN** — replaced; kept for history.

## Index

| ADR | Title | Status |
|---|---|---|
| [0001](0001-memory-bounded-ensemble-analysis.md) | Memory-bounded ensemble analysis (mode 2) | Accepted |
| [0002](0002-model-adapter-boundary.md) | Model-adapter boundary between orchestration and model internals | Proposed |
| [0003](0003-hierarchical-mpi-communicator-topology.md) | Hierarchical ensemble/model MPI communicator topology | Proposed |
| [0004](0004-distributed-state-ownership.md) | Distributed, rank-local state ownership (no global materialization) | Proposed |
| [0005](0005-state-ownership-replicated-vs-distributed.md) | State ownership: `replicated` vs `distributed` | Accepted |

## Writing a new ADR

Copy the shape of an existing one: context (what problem/evidence motivated
this), decision, consequences (including what it rules out), and status.
Reference the audit or benchmark evidence that motivated the decision where
one exists, rather than asserting performance claims without it.
