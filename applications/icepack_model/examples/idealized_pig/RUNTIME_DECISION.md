# Idealized PIG runtime validation — 9 November 2026

## Current status

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
