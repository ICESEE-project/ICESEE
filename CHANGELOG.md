# Changelog

All notable changes to ICESEE are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/);
versioning follows the version currently declared in `pyproject.toml`
(pre-1.0, so minor versions may include breaking changes).

This file starts tracking changes going forward; it is not a reconstruction
of prior git history. See `git log` for the full history before this file
existed.

## [Unreleased]

### Changed
- `idealized_pig` (icepack): restored the gate on the per-timestep flowline
  diagnostic (`applications/icepack_model/examples/idealized_pig/_icepack_model.py`)
  so it only runs at the configured save steps instead of every timestep for
  every ensemble member. No change to written output.
- Execution mode 2 (`src/parallelization/EnKF_parallel_io.py`,
  `src/parallelization/_mpi_ensemble_intialization.py`): fixed a
  collective-mismatch deadlock risk in the file-backed ensemble I/O. The
  HDF5 batch window is now only ever opened/extended by a call every rank
  in the shared communicator reaches (`read_forecast`, and two explicit
  priming calls added at ensemble-init time); the subcommunicator-root-only
  write path (`write_forecast`) now asserts the window is already open
  instead of possibly attempting the collective itself. Also raised the
  minimum batch window to 2 timesteps so the paired read/write for a
  timestep always shares one window instead of forcing a rebuild on most
  steps. No change to the ensemble-analysis mathematics.

- Added `src/utils/state_ownership.py`: a lightweight `StateOwnership`
  descriptor (`replicated` vs `distributed`) and `resolve_state_ownership`/
  `combine_member_state` helpers, used wherever a model communicator's
  per-rank state sizes are combined into a member's global state. Default
  is `distributed` (preserves the exact arithmetic every existing model
  adapter already relied on); Lorenz-96 opts into `state_distribution:
  replicated` in its own `params.yaml`, since every rank in a Lorenz-96
  model communicator computes/holds the identical full state rather than a
  partition of it.
- `src/parallelization/_mpi_generate_true_wrong_state.py`,
  `src/parallelization/_mpi_ensemble_intialization.py`,
  `src/utils/tools.py` (`icesee_get_index`): fixed the state-size/indexing
  bug this replicated-vs-distributed distinction exists to prevent.
  `global_shape = sum(dim_list)` and `subcomm.gather(...)` +
  `vstack`/`hstack` previously assumed every rank's reported state size was
  a disjoint partition to concatenate; for a replicated model (every rank
  reports the same full size) this double/triple-counted the state and
  duplicated data instead of assembling a partition. `icesee_get_index`
  separately assumed `dim_list` was always a per-rank partition array
  indexable by `comm.Get_rank()`, which is also wrong for replicated state.
  Deleted the dead, byte-for-byte-unused duplicate function
  `generate_true_wrong_state_full_parallel` (never called anywhere) while
  fixing its live twin, to remove the ambiguity of two copies of the same
  bug.
- `src/parallelization/_mpi_generate_true_wrong_state.py`: fixed three
  further pre-existing bugs found while making the `Nens < P` (multi-rank-
  per-member) path reach the state-size fix above end-to-end for the first
  time: (1) the true/nurged-state HDF5 file was opened with
  `driver='mpio'` scoped to a member's subcommunicator but only entered by
  that subcommunicator's root rank -- a collective-mismatch hang identical
  in kind to the mode-2 I/O fix above, fixed by opening with `COMM_SELF`
  since the write is genuinely single-rank; (2) every ensemble member
  independently generated **and wrote** the true/nurged trajectory to the
  same shared file, racing every other member to create the same HDF5
  dataset name (`ValueError: Unable to synchronously create dataset (name
  already exists)`) -- the true/nurged trajectory is a single reference
  the whole ensemble is compared against (matching how the `Nens >= P`
  branch already treats it), not a per-member quantity, so only member 0
  now generates and writes it; (3) a vestigial `comm_world.bcast` of a
  shape that was never read by anything downstream, which became a
  genuine rank-subset collective mismatch once (2) restricted the
  surrounding block to member 0 -- removed rather than rescoped.
- `src/parallelization/_mpi_ensemble_intialization.py`: the `Nens < P`
  ensemble-initialization branch called `model_module.initialize_ensemble`
  directly, without the `_member_initialization_context` wrapper its
  sibling `Nens >= P` branch already uses to provide a lightweight
  `statevec_ens` shape proxy; any adapter reading
  `icesee_kwargs["statevec_ens"]` (e.g. Lorenz-96) raised a bare
  `KeyError`. Fixed to use the same wrapper.

## Stage 3 stabilization (continued)

The prior entry's "known issues" are resolved below; the `Nens < P` path
now completes genuinely end-to-end (see the new Level-3 test) rather than
hanging, false-passing, or crashing partway through.

### Changed
- `src/run_model_da/icesee_da_full_parallel.py`
  (`icesee_model_data_assimilation_full_parallel`): fixed the top-level
  exception handler, which was the root cause of a crashed run being able
  to report success. It performed an unconditional `comm_world.Barrier()`
  plus further collectives (`enkf_parallel_io.close()` on a collectively-
  opened HDF5 handle, a `finalize_stack()` recovery-view build) that
  assumed every rank reached the handler in a matched state -- not
  guaranteed, since an exception is often raised by a rank subset -- and
  then, after printing a message, fell through and returned normally with
  no re-raise and no abort, so a fully crashed job (every rank raising)
  could still exit 0. Fixed to the standard MPI failure pattern: print
  full diagnostics from every rank that hit the handler (never suppress
  the original exception), attempt only a rank-local (rank 0, no
  collective) best-effort checkpoint save, then `comm_world.Abort(1)`
  immediately. The automatic collective "recovery view" build on error was
  removed rather than made conditional -- it cannot be done safely without
  assuming the rest of the communicator is reachable; a user can still
  build one explicitly (`finalize_stack`) from a stopped run's partial
  output.
- `src/parallelization/_mpi_ensemble_intialization.py` (`Nens < P` branch):
  fixed three further bugs blocking genuine completion, found via the new
  Level-3 test (`src/tests/test_lorenz96_mode2_p_nens_matrix.py`):
  1. `hdim` was computed as `initial_data[key].shape[0] //
     num_state_vars` -- dividing one variable's own already-per-variable
     block length by `num_state_vars` a second time (for Lorenz-96: a
     length-1 block // 3 = 0 via integer floor division). `hdim` does not
     depend on the variable and must come from the member's total state
     size (`icesee_kwargs["global_shape"] // num_state_vars`), matching
     every other correct use of this quantity in this codebase; hoisted
     out of the per-variable loop accordingly.
  2. `generate_enkf_field` reads the per-variable block length from
     `icesee_kwargs["noise_dim"]`, not `"hdim"` (confirmed against the
     correct usage in `generate_initial_member_increment`, same file) --
     the wrong key name meant the function always saw `noise_dim=None`
     and raised `ValueError: generate_enkf_field requires 'hdim'.`
  3. `ii_sig` (the variable index) was passed as `None` instead of this
     variable's actual index in the loop; with `ii_sig=None`,
     `generate_enkf_field` generates every variable's noise in one call
     (`(hdim * num_vars,)`) rather than just this one variable's `(hdim,)`
     slice, raising `ValueError: non-broadcastable output operand with
     shape (1,) doesn't match the broadcast shape (3,)`. Fixed to pass
     each variable's index from `enumerate(key_list)`.
- `src/utils/state_ownership.py`: added `ensure_state_array` (`np.atleast_1d`),
  a generic (not model-name-specific) normalization for the case a model
  adapter returns a bare NumPy scalar for a single-element variable block
  (`hdim == 1`) rather than a one-element array -- the same logical
  quantity, different NumPy representation. Applied in
  `_mpi_ensemble_intialization.py`'s `Nens < P` branch, which previously
  raised `IndexError: tuple index out of range` calling `.shape[0]` on
  such a scalar (Lorenz-96's `initialize_ensemble` returns `u0b[0]`, a
  `numpy.float64`).
- `src/parallelization/_mpi_generate_true_wrong_state.py`: added an
  `[ICESEE][ownership] distribution=... local_size=... global_size=...
  model_comm_size=...` diagnostic line so the replicated/distributed
  arithmetic is directly observable in a run's own output, not only
  inferable from code.

### Added
- `src/tests/test_state_ownership.py` (Level 1): unit tests for
  `StateOwnership` resolution, `combine_member_state`, and
  `ensure_state_array`, including the replicated invariant
  (`global_size == local_size` regardless of model-communicator size) and
  the distributed invariant (`global_size == sum(local_size_r)` for
  disjoint partitions), using a lightweight in-process fake communicator
  (no real MPI launch needed for this tier).
- `src/tests/_mpi_launcher.py`: shared `find_compatible_mpi_launcher()`
  helper for any test that launches real MPI jobs via subprocess. Not a
  version-string check -- a functional one (launches a real 2-rank
  rank-identification probe) -- because a PETSc/Firedrake-bundled `mpirun`
  can report the same vendor name as a separately-built one on `PATH`
  while being ABI-incompatible with this Python's `mpi4py`; confirmed live
  while building this helper (`shutil.which("mpirun")` picked the
  incompatible one, silently running a "2-rank" job as 2 disconnected
  singletons instead of raising).
- `src/tests/test_mpi_failure_handling.py` (Level 2) +
  `src/tests/parallel_mpi/_intentional_failure_worker.py` /
  `_healthy_worker.py`: proves one rank raising an exception terminates
  the whole job without hanging, with a nonzero result and a visible
  traceback (the regression guard for the exception-handling fix above),
  plus a negative control (a healthy job must still exit 0) so the test
  demonstrably distinguishes PASS from FAIL.
- `src/tests/test_lorenz96_mode2_p_nens_matrix.py` (Level 3): the real
  `run_da_lorenz96.py` mode-2 driver across
  `P ∈ {1,2,3,4,8,10,16}, Nens=4`, each validated against exit code,
  absence of error/traceback markers, the replicated-ownership invariant
  read from the new diagnostic line, exact expected output shape
  (`(3, 1001)` -- a full 1000-timestep run, not a silent truncation), and
  finite values. All 7 cases pass genuinely as of this stage.
- `docs/testing.md`: the test-tier hierarchy (Level 1-5) and the PASS/FAIL
  determination rules above, formalized.
- `docs/architecture/0005-state-ownership-replicated-vs-distributed.md`:
  Accepted (implemented) ADR for the ownership distinction.
- `paper/notes/replicated-vs-distributed-state.md`: architecture-only
  design note (no performance claims) on the same distinction, for the
  manuscript's model-independent-architecture section.
- `.github/workflows/ci.yml`: added the three new fast test files to the
  CI pytest step.

### Known issues (not yet fixed)
- CI's existing `scripts/ci/run_lorenz96_ci.py` selects its MPI launcher
  via a bare `shutil.which("mpirun") or shutil.which("mpiexec")`, the same
  pattern confirmed unreliable while building `src/tests/_mpi_launcher.py`
  (see above). Not fixed in this stage -- it did not block this stage's
  objective (the new tests use the compatibility-checked helper instead) --
  but CI environments with more than one MPI installation on `PATH` should
  be audited before trusting that script's mode-1/mode-2 CI runs.
- `EnKF_parallel_io.py`'s `_requires_shared_finalization` still hardcodes
  `model_name == "issm"` in shared runtime code rather than an adapter
  capability (identified in the prior Stage 3 investigation, not this
  stabilization stage). Still not fixed; still recorded.
- `ranks_per_model` / block-contiguous communicator topology (ADR 0003) is
  still not implemented. The `Nens < P` path being genuinely functional
  now (this stage) removes the blocker that made that work untestable;
  Stage 3D-3H is the recommended next stage.
