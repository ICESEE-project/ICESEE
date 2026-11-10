# Idealized PIG runtime validation — 10 November 2026

## Current status

The initial-state and observation/process-noise clock defects are corrected.
Modes 1–3 completed the five-year, 32-member PIG simulations with aligned
2/3/4-year analyses. Output verification and native restart checks passed;
post-run verifier allocation timeouts are recorded separately below.
See [the corrected-clock validation](README.md#corrected-clock-validation)
for the complete measured results and limits.

The measurements below predate the corrected state-column clock. They establish
the stated runtime/output properties for their recorded revisions, not timing
correctness of those older trajectories.

## Earlier runtime validations

Native Mode 3 completed the 32-member, five-year PACE validation in job
15125921.2: 100 forecast steps, three analyses, exit zero in 34 minutes 24
seconds. Final checkpoint fields and all 33 truth/member profile files were
checked. The batch launcher is `run_idealized_1.sbatch`; see
[MODE3_PACE.md](MODE3_PACE.md) for fetching, submitting and reading results.

Modes 1 and 2 subsequently completed the same 32-member, five-year case in
job 15142488 after the history-layout and finalization repairs. Mode 1 exited
zero in 42:01 and Mode 2 in 1:43:05. Both final fields/means and all 66 profile
files were checked. Per-variable final member relative L2 differences were at
most 1.54e-12; no new acceptance tolerance was inferred from this measurement.
Mode 2 still spends substantial time writing forecasts, so completion does not
establish efficient scaling. The earlier failed outputs remain preserved.

All three modes passed earlier concurrent one-year, two-member diagnostics.
Native Mode 3 also resumed from checkpoint 14, crossed analysis event 15 and
reproduced uninterrupted final fields. Native PIG execution is selected with
`--execution_mode=3`; no separate opt-in is required and it never falls back to
Mode 2.
Initial generic perturbations in Modes 1/2 and unperturbed native initialization
differ, so PIG member-wise equivalence across all three modes is not claimed.
The production check is not a controlled strong/weak scaling result.

## Direct Mode 3 selection validation

Job 15369255 completed the 24-member/eight-rank, five-year PIG
configuration on one PACE node using `--execution_mode=3` alone. Slurm reported
exit zero in 47:09 (MPI step 46:48). The final checkpoint records completed step
99, next step 100, model time 5.0 and three completed analyses. All eight shards
cover every owned field of all 24 members without gaps or duplicate ownership.
All final values are finite: thickness 30..4181.62 m, velocity components within
approximately -3867..974 in model units, surface 3.13..2632.13 m and basal melt
0..20. All 25 profile files have 415 points and four written columns at years
0, 2, 3 and 4. Solver settings and scientific configuration were not changed
while removing the additional selection gate. This validates completion and
inspected output integrity for this PACE case, not a controlled scaling result
or independent scientific accuracy.

## Cross-installation comparison

Job 15371824 repeated the 24-member/eight-rank case in a second PACE installation.
It completed all 100 steps and three analyses with exit zero in 46:14 (MPI45:54).
All eight final shards, five fields for every member, and all 25 profiles passed
output checks. Runtime/configuration content, input dataset values, generated
observations and error arrays matched the reference; truth profiles were exact.

The runs are not bitwise equivalent. Maximum final thickness difference was
0.211592 m (RMSE0.0151454 m, relative L2 8.57429e-6); maximum member profile
difference was0.0135394 m. Larger differences appeared after assimilation.
Python/Open MPI versions differed; the cause was not isolated. These measurements
do not define a scientific acceptance tolerance.

## Runtime repairs

- The Mode-2 initialization reuse decision is made once on world rank zero
  and broadcast, preventing ranks from entering different collectives while
  step-zero files are being created.
- Analysis reads close writable HDF5 windows and reopen immutable forecasts
  read-only. Each rank writes a private serial analysis shard; bounded assembly
  retains the existing finalization and atomic publication contract. File locking
  is retained.
- Mode 1 no longer constructs an unused state-by-state process covariance.
  The global vector dimension is still normalized after ownership resolution.
- Icepack progress messages report the actual MPI communicator size. External
  MATLAB/ISSM process resource estimates are separate and unchanged.
- Single-rank Firedrake contexts use its managed `COMM_SELF`; distributed model
  groups keep their own communicators. This avoids PyOP2 teardown callbacks
  attached to expired single-rank split communicators.
- Serial observation IDs are mapped by physical coordinates to distributed
  model rows before routing observations or stochastic increments.
- Native forecasts apply the existing configured process-noise generator,
  member IDs, seeds and schedule. No initial thickness perturbation was added.
  Deterministic AR history is replayed on restart. Policy identity distinguishes
  these checkpoints from older runs without native process noise.
- Native persistent and streamed analysis paths accept the same bounded row
  controls as Modes 1/2, including the configured inflation.
- Deterministic graph transition matrices are cached with a 64-MiB per-process
  limit. Keys include coordinate contents, neighbor count and tolerance.
  Random fields, seeds and normalization are unchanged.

The native noise bridge uses a model-group leader's canonical stochastic
workspace and scatters increments to owned rows. Physical member state is
never gathered. This bridge has an O(global state size) workspace with a
256-MiB allocation guard; it is not fully distributed random-field generation.
Its temporary AR files live in node-local scratch.

Compilation caches use node-local scratch, isolated by job and Slurm step.
XDG/SQLite state is rank-private. Persisting compilation caches on a shared
project filesystem is not required for these validations.

## Scientific settings and solver behavior

Scientific `features_debug` changes are retained. The runtime repair does not
raise the PETSc iteration cap or change physical parameters, uncertainty
amplitudes, observation-error method or configured inflation values.
The legacy Icepack options `max_iterations` and `tolerance` are not translated
to different effective PETSc settings.

The round-off-spread guard remains active. A future configuration with no
defined ensemble spread is refused collectively rather than producing a
large transform from numerical noise.

## Evidence

PACE job 15114944 completed all three concurrent short PIG cases with exit zero.
Every truth/member profile file was present and finite.

PACE job 15115483 completed the concurrent one-year cases and native restart
check with exit zero. Modes 1/2 retained finite `270615 x 2 x 21` histories,
with maximum absolute member-wise difference `4.05e-11`. Native fields were
finite at model time 1.0 with two completed analyses. Restarting an isolated
copy from checkpoint 14 gave maximum final difference `4.09e-12` and relative
L2 difference `1.64e-16` from the uninterrupted run. The original output was
preserved.

The current compact initialization input was 29.2 MB. Historical statements
about reading a 37-GB spin-up file do not describe that input. Timings include
setup and use different allocations; they are not controlled speedup estimates.

## Profile files

Under the selected `--data_path` (normally `_modelrun_datasets`), the files are:

- `hs_profiles_true.h5`
- `hs_ensemble_profiles0.h5`, `hs_ensemble_profiles1.h5`, and subsequent member IDs

These are HDF5 files, not folders or plain text. The diagnostic tests produced
415 sampling points and four profile columns. `written_steps` and
`profile_times`, where present, distinguish saved times from unwritten columns.
Names without an extension may belong to older runs.

## Reproducing the checks

Set `ICESEE_SOURCE_ROOT` to the reviewed scientific checkout and
`ICESEE_SPACK_ROOT` to the installed dependency environment. The scripts refuse
existing output targets and record source status.

- `pace/test_modes_1_3.sbatch`: short concurrent diagnostics.
- `pace/test_modes_1_3_year.sbatch`: one-year cases plus isolated mid-run restart.
- `run_idealized_1.sbatch`: Mode-3 launcher submitted from the example directory;
  defaults to 24 members, eight ranks on one node, 48 GB and two hours. It reuses
  the retained dependency installation and saved external YAML when present.
- `pace/validate_production_modes.sbatch`: concurrent 32-member, five-year
  confirmation; job 15125921 produced the mixed terminal results above.

Verify the intended remote with `git remote -v`; inherited worktrees can have
`origin` pointing to the parent environment repository. Scientific source is
published to ICESEE-project/ICESEE. Do not reset or replace another active
checkout to migrate this candidate.

## Outstanding acceptance checks

Measure remaining Mode-2 forecast I/O, validate additional spatial topologies,
compare native scientific configurations with
matched initial ensembles, and perform controlled repeated strong/weak scaling.
The validated native case used eight model groups, four members per group and
two spatial ranks per group. Final thickness was 30..4181.29 m; all five fields
were finite. Finiteness alone does not establish scientific accuracy.
The synthetic example's shared whole-member native binding supports one
spatial rank per model and rejects resume until it has a validated restart
binding. That binding is not evidence of general spatial decomposition for
every model.

## Full ensemble history for plotting

The current PIG launcher retains every native cycle and finalizes a standalone
`icesee_ensemble_data.h5` after successful completion. Export reads bounded
chunks from closed shards, aligns state rows to the serial truth/observation
numbering and publishes the aggregate atomically. Checkpoints remain for restart.
For five years with `dt=0.05`, checkpoint 99 is only the final snapshot; the
aggregate has 101 time columns including the initial ensemble. Existing rolling
archives cannot supply pruned intermediate snapshots. This storage change does
not introduce a different filtering method or solver/physical configuration.
See the README for the exact schema, retention requirements and validation scope.


## Exporting a failed run's available history

The common batch launcher attempts serial finalization after a failed MPI step,
using only consecutive committed checkpoints. Its `exit-status` remains the
simulation's original failure code; `finalizer-exit-status` and `finalizer.log`
report recovery separately. For a manually launched or forcibly interrupted run,
stop all its processes and run the standalone exporter with `--allow-partial`.
The README gives the exact command and the partial-file schema.

A partial file has `export_complete=True`, `complete=False`, explicit
requested/completed counts and a `t` axis ending at the last committed time.
Initial-only recovery is supported. Missing middle snapshots, missing shards,
invalid numbering and nonfinite values are errors. A complete published file
is protected from partial replacement; `--output` can select another target.
The notebook and generic output reader can use the partial file without the
end-of-run summary. Hard node/allocation termination may prevent automatic
recovery, so manual finalization remains necessary in that case.

The historical full-history run exposed an extra initial forecast and a
one-cycle observation offset. Its initial thickness matched truth at `dt=.05`
instead of zero, and observations at 2/3/4 years updated states at
2.05/3.05/4.05. Those saved values and timestamps remain unchanged.

The corrected clock copies initial physical fields without advancing, analyzes
observation column `j` after cycle `j-1`, and applies configured process noise
before that analysis. Restart metadata and native noise identities distinguish
this contract from older output. Choose a fresh data path when upgrading;
initial-time observations require an unimplemented pre-forecast analysis stage
and are rejected. See the README's time-grid section for details.

PACE job **15432397** completed the one-year, two-member tests in Modes 1–3,
a four-rank noise ownership/replay check, and a native restart. All steps exited
zero. Analyses occurred at 0.5, 0.75 and 1.0 years, including the final state.
Restart from checkpoint 13 through 19 crossed two analysis events; its final
relative L2 difference from the uninterrupted run was 6.57e-13 (maximum absolute
difference 2.99e-9), within the declared `rtol=1e-9`, `atol=1e-7` checks.
This replay used copied completed output, not an interrupted partial-profile test.

The three five-year, 32-member simulations in job **15433074** also exited zero:

| Execution mode | MPI ranks / ranks per model | Simulation step | Elapsed |
| --- | --- | --- | --- |
| 1 | 8 / 1 | 15433074.2 | 42:15 |
| 2 | 8 / 1 | 15433074.0 | 1:38:31 |
| 3 | 16 / 2 | 15433074.1 | 37:11 |

These concurrent runs shared a two-node allocation; these times are operational
measurements, not a controlled scaling comparison. Each output has ensemble
shape `(270615, 32, 101)` at `dt=0.05`, three analyses at years 2/3/4
(columns 40/60/80), and 33 finite truth/member profile files of shape `(415, 4)`
with written times 0/2/3/4. Native output retains the initial condition and all
100 committed forecast steps; checkpoint 99 records next step 100 and time 5.

The original allocation **timed out during post-run verification**, after all
three simulations had finished successfully and the Mode 1–2 output scans had
passed. Verification-only job **15439331** also timed out in repeated virtual
HDF5 reads. These allocation failures are retained in the record. A subsequent
read-only row-slab reader reproduced every short-case check and completed the
remaining full native field/mean/profile/manifest checks. Together those checks
took 67.37 seconds wall time; this is reader validation, not model runtime.
No simulation was rerun to obtain that result.

All final Mode 1–2 members were compared in bounded slabs: the largest per-field
relative L2 difference was 9.12e-13. Native initial thickness and surface exactly
matched truth at zero; initial velocity differences were at most 9.10e-13.
The native five-year thickness range was 29.94–4180.72 m, with final values
30.00–4180.68 m; velocity components stayed within approximately -7013 to 1550
in model units. Assimilation produced transient basal-melt values -0.3664 to
20.0441, returning to 0–20 at the final state; no unreviewed clipping was added.
These checks establish timing, output integrity and the stated restart behavior,
not scientific accuracy. Generic and native initial perturbation paths differ,
so the three modes are not claimed to have identical member trajectories.

To reproduce the prepared cases from this example directory (requires a PACE
allocation and existing lab inputs):

```bash
export ICESEE_SOURCE_ROOT="$(git rev-parse --show-toplevel)"
export ICESEE_SPACK_ROOT="$HOME/r-arobel3-0/ICESEE-Spack"
sbatch pace/test_modes_1_3_year.sbatch
# After checking the short case:
sbatch pace/validate_production_modes.sbatch
```

Both scripts run `tools/verify_timing_outputs.py` on their closed outputs. This
checker targets these fixed observation/profile schedules and geometry; it is
not a general validator for arbitrary PIG configurations.


Job 15384953 (`da6f9dde`, 24 members/eight ranks) completed the full five-year
history path in 1:29:27 with every step exiting zero. The aggregate contains
270615 state rows, 24 members and 101 time columns. Independent reconstruction
compared every one of its 655,970,760 values and every mean exactly with the
checkpoint shards; all 25 profiles passed shape, written-time and finite-value
checks. Export took 2485.64 seconds in that baseline. The revised contiguous
reader/recovery implementation is verified separately on the closed archive;
a new PACE performance result is not implied. See the README for the measured
scope of the historical storage test and the separate corrected-clock validation.
