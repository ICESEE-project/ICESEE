# Testing ICESEE

## Why this document exists

A real crashed mode-2 MPI run was observed to exit with status 0 and no
visible indication of failure short of reading console output closely: the
exception was caught, printed, and then silently swallowed with no
re-raise and no `MPI_Abort`. A test harness (or CI step) checking only the
process return code would have reported that run as a pass. This document
exists so that mistake is not repeated in either the runtime error-handling
convention or in how tests are written against it.

**A test is successful only if the program actually reached normal
completion.** For any MPI-launched test, that means checking, at minimum:

1. process return code (`0`, and specifically checked, not assumed),
2. absence of error/traceback markers in stdout+stderr, even when (1)
   passes — a crash can still (in principle, or under a launcher/environment
   quirk) fail to produce a nonzero code, so this is a genuinely independent
   check, not a restatement of (1),
3. the expected completion marker/result (e.g. a specific output dataset,
   an expected number of processed timesteps),
4. expected output dimensions/shape,
5. numerical sanity where applicable (finite values, at minimum).

A bounded `timeout` on every subprocess-launched MPI test is not optional:
without one, a hang looks identical to "still running" forever, which is
worse than an ordinary failure because it blocks the whole test run
instead of reporting anything.

## ICESEE's MPI failure convention

`icesee_da_full_parallel.py`'s top-level exception handler is the
reference implementation: on an unrecoverable error, print full
diagnostics (rank number, exception, traceback) from every rank that hit
it, attempt only rank-local, non-collective cleanup (e.g. a checkpoint
save from rank 0), then call `comm_world.Abort(1)`. It deliberately does
**not** attempt a `Barrier()`, a collective file close, or any other
collective operation, because an exception on one rank does not guarantee
any other rank is still reachable in a matched state — see ADR discussion
in `CHANGELOG.md` under the Stage-3 stabilization entry, and
`src/tests/test_mpi_failure_handling.py` for the regression test proving
this pattern prevents a healthy rank from hanging forever at a barrier a
failed rank will never reach.

Any new collective-heavy MPI code path in this codebase should follow the
same shape: matched collectives only where every intended participant is
guaranteed reachable; on failure, diagnose loudly and abort, never try to
gracefully coordinate a shutdown across a communicator that may already be
in an inconsistent state.

## Test tiers

### Level 1 — Unit tests

Pure-Python logic with no MPI launch required (may still import `mpi4py`
and use `MPI.COMM_SELF`/small in-process fake communicators to exercise
multi-rank *arithmetic* without real parallelism). Fast, run on every
change.

Example: `src/tests/test_state_ownership.py` — `StateOwnership`
resolution arithmetic, `ensure_state_array` scalar/array normalization,
using a lightweight fake communicator to simulate several ranks' worth of
`allgather`/`gather` results without launching real MPI processes.

### Level 2 — Synthetic MPI tests

Real multi-rank jobs launched via `mpirun`/`mpiexec` as a subprocess, but
against a minimal purpose-built worker script, not a full model/driver.
Used to isolate one MPI-architecture property (a collective's
safety/matching, a communicator topology, deliberate failure handling)
from everything else that a real model run would also exercise.

Example: `src/tests/test_mpi_failure_handling.py` — proves one rank
raising an exception terminates the whole job (no hang) with a nonzero
result and a visible traceback, using a trivial two-line worker, plus a
negative control (a healthy job must still exit 0) so the test suite is
demonstrated to distinguish PASS from FAIL, not just report failure
unconditionally.

**Launcher selection for any Level 2/3+ test**: use
`src/tests/_mpi_launcher.py`'s `find_compatible_mpi_launcher()`, never a
bare `shutil.which("mpirun")`. A sandbox/HPC environment can have more
than one MPI installation on `PATH` (e.g. one bundled inside a
PETSc/Firedrake build, another from Homebrew or the system package
manager); `shutil.which` takes whichever is first with no compatibility
guarantee. A mismatch between the launcher and the Python's `mpi4py` build
does not raise — it was observed to silently run the job as N disconnected
single-rank processes instead of one real N-rank job (each identifying as
rank 0, warned only via a "suspicious MPI execution environment" message
on stderr). That would make a test pass for the wrong reason: no real
communicator, so no real collective to hang on or fail on.
`find_compatible_mpi_launcher()` verifies this functionally (launches a
real 2-rank rank-identification probe) rather than trusting a version
string, because a bundled build can genuinely report the same vendor name
as the correct one while still being ABI-incompatible.

### Level 3 — Lightweight real-model integration

The actual ICESEE driver, against a small/fast model (Lorenz-96), across
representative MPI topologies. This is where model-adapter and driver code
gets exercised end-to-end, cheaply enough to run routinely.

Example: `src/tests/test_lorenz96_mode2_p_nens_matrix.py` — the real
`run_da_lorenz96.py` entry point under mode 2, across
`P ∈ {1,2,3,4,8,10,16}` at `Nens=4` (including every `P > Nens`,
multi-rank-per-member case), each validated against all five checks in
"Why this document exists" above, plus the `replicated`-ownership
invariant (`global_size == local_size` regardless of model-communicator
size) read directly from the driver's own `[ICESEE][ownership]`
diagnostic line.

### Level 4 — Large-model integration

Reduced (not full-scale) Icepack/ISSM cases through the real driver. Not
yet built out for the `Nens < P` / multi-rank-per-member path — Level 3
established the pattern; Level 4 needs it repeated against a genuinely
`distributed`-ownership adapter (Icepack), which Level 3 cannot cover
(Lorenz-96 is `replicated`). Do not claim Icepack multi-rank correctness
from Level 1-3 results — architecture/static compatibility is established
(see ADR 0005), real MPI integration validation is not.

### Level 5 — Performance/scaling

Strong, weak, ensemble, I/O, and memory scaling studies. Explicitly out of
scope until the architecture is stabilized (Levels 1-4 genuinely green) —
see the audit and Stage 2/3 reports. **Do not use Level 5 to discover
Level 1/2 correctness bugs**: a scaling run that hangs or silently
truncates is not scaling data, it is an undiagnosed Level 1-3 bug wearing
a scaling experiment's clothes. Every bug found and fixed during Stage 3
(state-size arithmetic, two distinct collective mismatches, a shared-file
write race, a missing context wrapper, a wrong-key/wrong-shape noise-field
call) would have been indistinguishable from "this configuration doesn't
scale" if first encountered inside a scaling run instead of a bounded,
assertion-rich Level 3 test.

## Running the suite

```
PYTHONPATH=<repo parent>:<repo root> pytest src/tests -q
```

CI (`.github/workflows/ci.yml`) currently runs a curated subset directly by
filename rather than the whole directory; keep that list in sync when
adding a new fast, CI-appropriate test file.
