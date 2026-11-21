# Bringing an application into ICESEE

An application defines a physical problem and binds it to ICESEE's estimation and execution machinery. It includes the input preparation, model initialization, state layout, observation interpretation, forecast step, constraints and output/restart policy. We bring an application on board by making these relationships explicit before selecting a larger process layout. A model-family registry entry and a runnable scientific application are different pieces of that work.

This guide follows runtime source `3748e416d3a8120cc860119bd7224c6f9ad41ff7`. The existing [adapter guide](model-integration.md) gives the shared contracts; [Icepack](coupling-icepack.md) and [ISSM](coupling-issm.md) show how different backends satisfy them. The [templates](templates/README.md) are deliberately incomplete, import-safe scaffolds. Their physics raises `NotImplementedError`; they are not installed model backends or a claim of a working assimilation experiment.

## Choose the kind of integration

| Intended change | Main work | What it does not establish |
| --- | --- | --- |
| A new scenario for an existing example | Reviewed configuration/inputs and a fresh run identity | New forcing or observation coverage still needs scientific validation; a shape-compatible file is insufficient. |
| A new example in an existing model family | An example package, entrypoint, initialization and example hooks | The family-level fallback may select another implementation unless discovery is checked. |
| A new backend/model family | A new family package, `MODEL_CONFIG` entry and the selected runners' hook contracts | A registry status of `supported` is a dispatch decision, not completed testing. |
| A whole-member native application | Persistent member/context adapter, one-rank spatial policy and explicit analysis policy | Native naming does not remove the complete-member memory requirement or add restart support. |
| A spatial native application | Owned field layout, communicator-aware physics, observation routing, finalization and a validated native runner | A lower-level reference adapter alone does not provide a production lifecycle. |

A first integration usually starts with a model-only forecast and a tiny serial assimilation. This isolates units, state order and observation meaning before MPI ownership complicates diagnosis. Mode 1 and Mode 2 can then test the same realized prior and observations through different orchestration/storage paths. Mode 3 is an additional capability, not a mandatory renaming of the generic array adapter.

## Define the scientific and array contracts

The application specification records the prognostic variables, diagnostic variables and parameters, each variable's units and global/local dimensions, and how these enter `vec_inputs`. It identifies which quantities the observations measure, their uncertainty model, spatial support and requested/mapped times. We also specify the model timestep, initial ensemble construction, process noise, constraints, inversion schedule and acceptance metrics.

The ordered packed state must mean the same thing in initialization, forecast results, observation indexing, analysis controls and output. A useful initial fixture gives each variable and member a distinct known pattern. An array that has the expected length but exchanges two blocks is scientifically wrong. Some shared row controls require equal block sizes, while other paths support `var_nd` or native segmented layouts. The existence of one flexible loader option does not make every downstream operation compatible with unequal blocks.

`state_distribution: replicated` means the model-communicator ranks hold the same complete state. `distributed` means their owned pieces are disjoint. Ghost copies are excluded from the latter. This distinction controls dimension assembly; it does not select a physical model or a filter. Stable member IDs, rather than rank IDs, identify initialization, stochastic streams, filenames and restored state.

## Resolve the example and hook module

[`SupportedModels`](../../applications/supported_models.py) accepts the family keys `icepack`, `issm`, `lorenz` and `flowline` in this snapshot. It lowercases the requested `model_name`, requires that family to be present with status `supported`, then uses the **current directory basename** to try:

```text
<model_name>_model.examples.<directory_basename>._<model_name>_enkf
```

For example, from `applications/icepack_model/examples/my_case`, `model_name: icepack` selects `icepack_model.examples.my_case._icepack_enkf` when that module exists. The fallback module for Icepack is `icepack_model.icepack_utils._icepack_enkf`. Directory names need valid Python identifiers. The `applications/` directory is placed on the import path by the loader. Example packages should include their package files and dependencies explicitly rather than depend on a coincidental shell path.

`example_name` in configuration is useful application metadata; it is not the expression used by this legacy discovery function. A new backend `coastal` follows the same naming convention with `coastal_model/examples/my_case/_coastal_enkf.py`, alongside a family fallback and its explicit `MODEL_CONFIG` entry. A caller can supply a custom `model_config` to `SupportedModels`, but that does not automatically change every runner's normal construction of the registry. Production integration needs a coherent registry and entrypoint.

The entrypoint therefore establishes its example directory before configuration/dispatch and records the selected module's actual `__file__`. Output-path conventions also belong to the backend: current ISSM code builds some Python–MATLAB paths relative to its application workspace and uses shared file IPC/`Models` directories. A unique output path alone does not isolate concurrent ISSM runs; its coupling guide specifies the additional workspace boundary. A missing example module can lead to fallback. An exception inside an existing imported example module is reported as an actual failure; hiding it by choosing another backend would make diagnosis unreliable. Native classes and layout objects should use consistent `ICESEE...` package imports to avoid duplicate Python class identities.

## Implement the generic boundary against active callers

The following signatures are present in the PIG example and its active generic forecast path. They are a concrete starting point rather than a universal formal class interface.

| Signature | Expected relationship to the runner |
| --- | --- |
| `initialize_model(**icesee_kwargs)` | Entrypoint/model setup. PIG returns a tuple of model fields/spaces/solver; synthetic Icepack returns an updated context dictionary. The entrypoint adapts the return deliberately. |
| `initialize_ensemble(ens, **icesee_kwargs)` | Returns this member's initialized blocks in the declared layout. The context can include a lightweight `statevec_ens` shape proxy; it is not necessarily a materialized ensemble. |
| `forecast_step_single(ensemble=None, **icesee_kwargs)` | Receives one packed member through the active `advance_ensemble_member` path and returns updated variable blocks. `k`/`ens_id` identify cycle/member. Exactly one forecast interval belongs here. |
| `generate_true_state(**icesee_kwargs)` | Writes initial column zero and subsequent reference states using the application's clock. A supplied file-backed target can be written in place and return `None`. |
| `generate_nurged_state(**icesee_kwargs)` | Optional unassimilated/wrong reference, with the same clock and target discipline. Preserve the repository's spelling in hook names. |

Optional hooks have caller-specific conventions. For example, the serial runner calls `post_analysis_update(ensemble_vec, **icesee_kwargs)` and expects the updated context and ensemble; its covariance override is called with `Cov_Obs_fun(sig_obs=..., nd=..., icesee_kwargs=...)`. Generic inversion paths call `inverse_step_single(ensemble=member, **member_kwargs)`. Optional `Obs_fun` and `JObs_fun` are passed into the filter where supported. A new implementation follows those actual consumers, rather than assuming a same-named hook is used identically in every mode. PIG does not define every optional hook in this table.

The [legacy hook scaffold](templates/legacy-hooks.py.in) preserves the core signatures and describes the required return/storage behavior without inventing physics. The [driver scaffold](templates/application-driver.py.in) keeps configuration loading inside `main`; importing it to inspect its interface does not load the cleanup-capable configuration module. The [YAML scaffold](templates/params.yaml.in) illustrates existing configuration sections and state dimensions. Its backend, physical meaning and numerical values remain unimplemented and unvalidated.

## Introduce parallel execution in stages

Mode 0 establishes the basic hook/clock contract with one process. Mode 1 distributes ensemble work while retaining in-memory orchestration. Mode 2 uses file-backed state and bounded analysis; the application must honor file-backed reference targets and coherent HDF5 writer/collective ownership. Neither mode automatically makes a full-member model callback spatially distributed. Multi-rank model solves require compatible communicator use and explicit ownership metadata in the selected backend.

Every solver collective in a model group needs all its participants in the same order. A world-root-only truth writer must not open a file collectively on an unreachable multi-rank subgroup. A spare rank has no model mesh and cannot construct one on `COMM_NULL`. Tests of more ranks than members, fewer ranks than members and non-divisible scheduling expose assumptions that a balanced two-member example can hide. The [MPI standard](https://www.mpi-forum.org/docs/mpi-5.0/mpi50-report.pdf) and [HDF5 collective requirements](https://support.hdfgroup.org/documentation/hdf5/latest/collective_calls.html) define communication/I/O semantics; the source identifies which communicator each call actually uses.

For a whole-member native route, `register_native_application(model_name, *, application, adapter, setup, analysis_cycle=None, supported_features=())` binds an adapter and setup to the shared lifecycle. The generic policy requires `model_nprocs=1` and refuses requested resume without a validated binding. A declared feature and an analysis policy must agree: unsupported localization/inference/inversion flags are rejected by the shared route. Synthetic Icepack is an example of an explicit file-backed native analysis policy, not a hidden dispatch back to Mode 2.

For a spatial native route, the [native adapter scaffold](templates/native-adapter.py.in) shows the required method signatures. `initialize_native_member` constructs persistent native fields on `topology.spatial_comm`. `forecast_native_member` advances those same fields. `observe_native_member` evaluates owned observation rows in the declared canonical order. `finalize_native_analysis` implements the chosen post-analysis policy and synchronization. Additional inversion, allocation/reactivation and checkpoint-restoration capabilities are only exposed after implementation and testing. A physics projection callback is optional in a concrete wrapper; its absence must be described, as in current native PIG.

`DistributedFieldRegistry` packs/unpacks owned field pieces, while `NativeDistributedMember` carries `member_id`, fields and opaque model context. Observation IDs require a mapping when reference and distributed meshes have different numbering. A checkpoint needs layout/policy/progress identity and a restoration route for both packed fields and necessary auxiliary model state. A complete-field NumPy adapter cannot acquire a no-global-member property solely by wrapping it in native class names.

Mode-3 dispatch looks up registration using `model_name` and the current directory basename. The explicit binding is `register_execution_mode_3(model_name, runner, *, notes='', application=None)`. It can infer an application from the runner module name, but an explicit `application='my_case'` avoids ambiguity. Entrypoints normally import their `mode3_runner`; dispatch also attempts an example-specific bootstrap import. An unregistered case fails with `NotImplementedError`. Registration alone is not proof that a runner supplies all scientific policies or restart semantics.

## Verify before expanding the workload

The acceptance sequence is model-only forecast, distinguishable-member initialization, known-state observation mapping, tiny analysis comparison, actual MPI execution, checkpoint/recovery, and a representative physical experiment. At the numerical comparison stage we use identical realized priors, observation/error factors, seeds and time grids; means alone cannot detect compensating member errors. Recovery crosses an analysis event and checks observation counts as well as final fields. An actual interrupted-model test is distinguished from replaying completed output or a synthetic checkpoint failure.

Existing infrastructure gates provide a reusable baseline:

```bash
cd "$ICESEE_SOURCE_ROOT"
python -m pytest -q src/tests/test_state_ownership.py src/tests/test_timestep_contract.py src/tests/test_distributed_adapter.py src/tests/test_distributed_fields.py
python -m pytest -q src/tests/test_distributed_native_runtime.py src/tests/test_distributed_checkpoint.py src/tests/test_mode3_checkpoint_recovery.py
python -m pytest -q src/tests/test_resource_plan_mpi_topology.py src/tests/test_mode3_checkpoint_recovery_mpi.py src/tests/test_distributed_history_export_mpi.py
```

These tests do not contain the new application's physics. The new case adds its own fixtures and backend-dependent checks. MPI tests use a functionally compatible launcher; backend or launcher skips are recorded as missing coverage. The template syntax can be checked without importing any model/configuration code:

```bash
python - <<'PY'
import ast
from pathlib import Path
import yaml
for path in Path('docs/handbook/templates').glob('*.py.in'):
    ast.parse(path.read_text(), filename=str(path))
yaml.safe_load(Path('docs/handbook/templates/params.yaml.in').read_text())
print('Scaffold syntax checked; no model was run')
PY
```

The final handoff contains the source/configuration/input manifest, exact entrypoint and launch topology, dependency versions, state/observation/unit diagram, backend lifecycle, validation results and a capability matrix. The [capability template](templates/capabilities.md.in) distinguishes implemented, tested and unvalidated combinations. Frequent mistakes are treating `example_name` as the legacy import selector, silently using the family fallback, counting ghosts, assuming all native routes have restart, changing initialization during a parity comparison, and importing application configuration merely to inspect defaults beside valuable output. A new application is ready for a larger run when these contracts and its scientific acceptance criteria are reviewable, rather than when its registry entry first imports.
