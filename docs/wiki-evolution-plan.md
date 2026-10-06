# Wiki evolution plan

The [ICESEE Wiki](https://github.com/ICESEE-project/ICESEE/wiki) is the
detailed user/developer/HPC knowledge base. README stays the concise entry
point, `docs/*.md` and `docs/architecture/` hold engineering-contract and
decision-record detail, and the Wiki holds user/developer/HPC how-to
content — see this repository's documentation-architecture discussion for
the full split. This plan maps current Wiki pages to
where the scalability architecture work belongs, without rewriting existing
pages wholesale. It is a plan, not an executed edit — the Wiki is a
separate git repository (`ICESEE.wiki.git`) and no push access/authorization
was exercised to produce this plan.

## Current Wiki structure (as of this review)

| Page | Covers today |
|---|---|
| Home | Overview, quick start, links |
| 1. Installation | Install instructions |
| 2. Usage | Running applications, EnKF variants |
| 3. Guide to Integrating Models | Hook contract, `icesee_kwargs`, directory layout — **this is the current model-adapter documentation** |
| 4. Build ICESEE as a package | PyPI packaging |
| 5. Development Notes | Project structure, namespace packages |
| 6. Common Issues and Solutions | Troubleshooting |
| 7 / 7.1. ISSM MATLAB Installation (macOS) | Platform-specific ISSM setup |
| 8. Short Tutorials | Worked examples |

## Proposed additions (new pages, not rewrites of the above)

1. **"9. Execution Modes"** — accurately document modes 0/1/2 as they
   actually behave (not aspirationally): serial, partial-parallel, and
   full-parallel, including the honest current caveat that mode 2's
   forecast-step concurrency depends on the `Ne`/`P` ratio (see the
   scalability audit, §E) until ADR 0003's `ranks_per_model` work lands.
   Update this page's mode-2 description once that lands, rather than
   writing the target state as if already true.
2. **"10. MPI Architecture"** — `COMM_WORLD` → ensemble communicator →
   model communicator, `P ~= Ne * ranks_per_model`, explicitly documenting
   the non-divisible case, cross-referencing ADR 0003 and
   `docs/architecture/`. Include the "4 ensembles + 4 MPI ranks means ~1
   model rank per ensemble member, i.e. no intra-model MPI decomposition"
   example verbatim — this is exactly the ambiguity the audit could not
   resolve from the repository alone, and is worth making unmissable.
3. **"11. Memory Model"** — local vs. ensemble vs. global state, what's
   static/immutable vs. per-timestep, chunked HDF5 access; largely draws
   from the already-written `docs/large-scale-execution.md` and ADR 0001,
   repackaged for the Wiki's user-facing audience rather than duplicated
   from scratch.
4. **"12. Large Multidimensional Data & Very Large Domains"** — explicit
   statement that the ~37 GB `idealized_pig` input dataset is a current
   test case, not an architectural ceiling; what changes (and what doesn't)
   as `Nx`/`Ny`/dataset size grow toward continental-scale domains.
5. **"13. Parallel I/O"** — file lifetime (the audit's batch-window
   findings), chunking, why files should not be reopened every timestep;
   cross-reference `docs/execution-mode-2-development.md`.
6. **"14. HPC Launch Patterns"** — worked `Nens`/`-n`/`ranks_per_model`/
   `ntasks-per-node` examples once ADR 0003 lands; until then, an explicit
   warning about the current strided rank-assignment behavior so users
   aren't surprised by node-locality effects.
7. **"15. Profiling and Performance"** — how to enable the instrumentation
   in `docs/performance-instrumentation.md` and read its output.
8. **"16. Adding Models"** — do not replace Wiki page 3; add a
   cross-reference to `docs/model-adapter.md`'s conceptual target section
   once/if the `local_state`/`set_local_state` capabilities are actually
   implemented, so page 3 stays the source of truth for the concrete API
   and doesn't drift into aspirational content.

## Sequencing

Pages 1 (Execution Modes) and 2 (MPI Architecture) should be written
**after**, not before, ADR 0003's rank-assignment fix lands and is
validated — writing them now would either document the current strided/
ambiguous behavior as if intentional, or describe the target as if already
true. Pages 3-5 above (Memory Model, Large Data, Parallel I/O) can be
written sooner since they describe the already-implemented and validated
mode-2 analysis design (ADR 0001).
