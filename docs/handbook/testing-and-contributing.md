# Testing and contributing safely

ICESEE tests must demonstrate the intended execution, not merely a successful launcher exit. A former error handler could print a crash and return zero; a mismatched MPI launcher could start multiple independent rank-zero processes. The current regression strategy therefore combines exit checks, error-marker checks, completed outputs, exact shapes, finite values, ownership checks and member-wise comparisons.

## Prepare an isolated test environment

Use the approved scientific environment or the CI dependencies, not an arbitrary system Python. The repository directory is itself the `ICESEE` package, and code uses both `ICESEE.src...` and `src...` imports. In a checkout directory named `ICESEE`, put both its parent and itself on `PYTHONPATH`:

```bash
export ICESEE_ROOT="$(pwd -P)"
export PYTHONPATH="$(dirname "$ICESEE_ROOT"):$ICESEE_ROOT${PYTHONPATH:+:$PYTHONPATH}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python -c 'import ICESEE; print(ICESEE.__file__)'
```

For a worktree with another basename, make a temporary parent containing a symlink named `ICESEE` and add that parent instead; do not rename the maintained working tree just to satisfy imports. Every test launching an application must give it a unique temporary `--data_path`. `src/tests/conftest.py` detects deleted/recreated application output directories and fails the test session if a test damages them.

## Useful focused checks

From the repository root, the following modules exercise central logic without requiring an Icepack/MATLAB production dataset:

```bash
python -m pytest -q \
  src/tests/test_state_ownership.py \
  src/tests/test_timestep_contract.py \
  src/tests/test_random_streams.py \
  src/tests/test_observation_error_terms.py \
  src/tests/test_analysis_stabilization.py \
  src/tests/test_distributed_analysis.py \
  src/tests/test_distributed_checkpoint.py \
  src/tests/test_distributed_fields.py \
  src/tests/test_native_process_noise.py
```

These tests may still import mpi4py or serial scientific libraries. “No MPI launch” does not mean “no dependencies.” Report skips separately. This focused list is a starting point, not the entire CI suite or a substitute for target-model verification.

Use the checked launcher helper before real MPI tests:

```bash
python - <<'PY'
from ICESEE.src.tests._mpi_launcher import find_compatible_mpi_launcher
launcher = find_compatible_mpi_launcher()
assert launcher is not None, 'No compatible MPI launcher'
print(launcher)
PY
```

Then run an appropriate bounded subset in the same approved environment:

```bash
python -m pytest -q -rs \
  src/tests/test_resource_plan_mpi_topology.py \
  src/tests/test_mpi_failure_handling.py \
  src/tests/test_lorenz96_mode2_p_nens_matrix.py \
  src/tests/test_enkf_parallel_io_topology.py \
  src/tests/test_enkf_readonly_analysis.py
```

Mode-2 topology tests need parallel HDF5/h5py. The per-launch timeouts inside tests are part of the contract; a scheduler allocation or outer process timeout should bound the entire test step too. A skipped real-MPI test means its behavior was not exercised. CI's MPI job explicitly rejects skipped MPI cases.

The full authoritative selection is [`.github/workflows/ci.yml`](../../.github/workflows/ci.yml): fast deterministic modules run separately from real-MPI modules and Firedrake/ISSM/data-dependent validation. Keep new tests in the appropriate tier rather than copying a stale count into the documentation. The workflow has pull-request and manual triggers; do not assume a direct push to every experimental branch automatically runs all jobs.

## What each validation layer establishes

| Layer | Evidence | Does not establish |
|---|---|---|
| Pure/helper tests | Indexing, shapes, stream identities, product algebra, policy validation | Correct MPI participation or model physics |
| Synthetic MPI tests | Actual communicators, collectives, failure aborts, ownership and shard contracts | Valid glacier inputs or scientific priors |
| Lightweight real-model integration | End-to-end drivers with complete outputs over a topology matrix | Distributed mesh correctness just because the model is replicated |
| Reduced target model | Actual Icepack/ISSM bridge, mesh and geometry behavior, model output/restart | Large-workload resource sufficiency or general scaling |
| Representative production case | That specific workload and environment completed and passed its protocol | Universal production readiness or scientific validity of unrelated settings |
| Scientific evaluation | Agreed physical/error diagnostics, independent observations or synthetic truth, sensitivity and uncertainty | General numerical equivalence for changed code/decompositions |

## Modes 1 and 2 parity

The paired harness creates fresh runs from the same full YAML configuration, forces float64/full retention, fingerprints source/config/output, verifies the launcher, and compares members without loading an entire history into memory. From the repository root:

```bash
python scripts/benchmarks/run_execution_mode_parity.py \
  --config docs/handbook/examples/lorenz-smoke.yaml \
  --application applications/lorenz_model/examples/lorenz96/run_da_lorenz96.py \
  --output /absolute/new/parity-case \
  --ranks 2 --nens 8 --mpiexec /absolute/compatible/mpirun
```

The harness reads the supplied YAML directly to make its mode files; use a fully resolved YAML rather than a relative `extends` chain that would move when copied. Its default pointwise gate is `atol=1e-10`, `rtol=1e-8`; float32 or changed scientific settings need a separate justified protocol. `--bitwise` is a fixed-platform diagnostic, not a portable promise across reduction orders.

The smaller CI plotting harness has a distinct gate: pointwise `1e-7` absolute/relative tolerance through the assimilation window, then per-variable full-history norm bounds. It also checks the truth, observations and physical clock. Describe the gate actually used; do not report one harness's tolerances as another's. Compensating member errors can preserve the ensemble mean, so mean-only agreement is inadequate.

For a target-model change, locate the first divergent stage: generated prior, model-only forecast, process-noise increment, active observation ordering/factor, transform, post-analysis control, geometry/inversion, or storage. `execution_parity_trace` is opt-in diagnostic output for this purpose and can consume extra storage. Preserve the first failing case before applying another change.

## Restart verification

Compare uninterrupted and interrupted/resumed cases with the same source/environment/configuration. Restore a checkpoint before an analysis event, continue across it, and verify final members, observation count, time grid and physical diagnostics. Check that the first resumed cycle is exactly the recorded `next_step`; no analysis/noise event should be duplicated. Corrupt or omit a disposable checkpoint shard in a test copy and ensure validation rejects it or selects a prior valid state. Never damage the production tree to test recovery.

Native PIG restart policy, generic Mode-2 rolling restart, and the default shared native policy are distinct. A unit test of checkpoint serialization is not a full adapter resume test. Layout/policy changes and cross-installation comparisons need their own recorded tolerance and protocol.

## Contribution review

Make one coherent change on a branch based on an exact baseline. State the triggering bug or desired behavior, the scientific/operational consequence, and what remains outside the tested scope. Include the smallest reproducing case, selected tests, meaningful comparison output and any updated documentation. Preserve existing user data and unrelated edits. Do not reformat large scientific modules while repairing one collective or row mapping.

If a change affects seeds, observation-error factors, priors, bounds, timing, state ordering or solver behavior, identify it explicitly as a numerical/scientific change. If it claims only storage/scheduling equivalence, provide same-input member-wise parity and, where relevant, restart evidence. A faster runtime with changed initialization is not a controlled speedup.

Keep architecture decisions current: mark a design accepted only for the scope actually implemented, and preserve unresolved limits. Refresh the generated key/API indexes with `python scripts/docs/build_reference.py`, review links, and update affected wiki guidance. The documentation should let a student locate the changed implementation and understand why its validation is sufficient for the intended use.
