# ICESEE production and development handbook

This handbook explains how to operate, investigate, and extend ICESEE while the principal developer is away. Start with a small reproducible run, retain the exact source and environment used for every experiment, and treat a completed calculation and a scientifically acceptable calculation as separate decisions. Students should be able to identify the code responsible for a result, reproduce a failure on a small case, and recover a stopped run without reconstructing the developer's working environment.

The implementation baseline is the published `features_debug` commit `3748e416d3a8120cc860119bd7224c6f9ad41ff7`, inspected on 10 November 2026. A branch name can move. Record the full commit and any local diff. The repository package metadata declares version `0.1.9`; the earlier paper and archive use `v1.0.0`. These identifiers have not been reconciled into a new release. Use the commit as the unambiguous identifier until maintainers resolve release metadata. The handbook is documentation of an inspected implementation, not certification of every model/configuration combination.

## Reading routes

| Reader and immediate task | Read in this order |
|---|---|
| New student running a first case | [Student tutorial](student-tutorial.md), [configuration](configuration.md), [operations](operations.md) |
| Student submitting an ice-sheet calculation | [Deployment](deployment.md), the example README, [operations](operations.md), [numerics](numerics.md) |
| Developer diagnosing a hang or unexpected state | [Architecture](architecture.md), [testing](testing-and-contributing.md), [module and API index](api-reference.md) |
| Developer integrating another model | [Model integration](model-integration.md), [state ownership](architecture.md#state-ownership), [configuration](configuration.md) |
| Researcher interpreting the updated developments | [Numerics](numerics.md), [architecture](architecture.md), the technical-report manuscript and its evidence registry |
| Local team deciding whether a case can proceed | [Handover status and acceptance checklist](production-readiness.md), [operations](operations.md), the selected application's evidence |

The existing [PIG README](../../applications/icepack_model/examples/idealized_pig/README.md) describes that application's input data, submit scripts, profiles, and inspected runtime results. Read it together with this handbook. A generic parameter default does not establish that an application implements the requested feature.

## Concepts to learn first

An **ensemble member** is a scientific realization identified by a stable integer, not an MPI rank. A **model communicator** contains ranks collaborating on one forecast solve. An **ensemble communicator** connects groups that hold corresponding pieces of different members. **Owned** field entries belong to this rank; **ghost** entries mirror neighbors and must not become extra analysis rows. A **state column** is the initial state or one completed model interval, while a **forecast cycle** is the zero-based interval advancing between two state columns.

`execution_mode` selects a runner. It does not choose the observation-error model, localization, physical parameters, or an appropriate MPI installation. Modes 1 and 2 have a same-input parity contract; Mode 3 is registered per application and has several distinct capability boundaries. Do not infer mathematical equivalence merely from selecting the same filter name.

## Production handover contract

Every production case should have an immutable run directory containing its source identity, resolved configuration, environment, input checksums, launch command, scheduler and MPI logs, output verification, and scientific acceptance decision. A dedicated working directory may be reused only for a deliberately supported restart. Fresh imports of the configuration loader can clean that directory; see [operations](operations.md#protecting-run-directories).

Choose local people to fill the operator, scientific reviewer, and software maintainer roles before sustained production. This handbook defines their decisions without assigning people who have not accepted responsibility. Two students should be able to reproduce the smoke case and restore the same representative checkpoint independently. Record that exercise with the run evidence. Operational independence comes from practicing the procedure as well as reading it.

## Where to find the implementation

The [architecture map](architecture.md#code-map) identifies the active drivers, analysis, storage, and adapter modules. The generated [API index](api-reference.md) covers Python modules outside the test suite and lists source symbols rather than presenting an artificial universal adapter. The [configuration key index](configuration-key-index.md) records actual loader fallback expressions with source locations. These generated files can be refreshed with `python scripts/docs/build_reference.py` from the repository root; review the diff when the code changes.

This handbook intentionally preserves historical names in commands. For example, `lorenz96` is the directory name of the lightweight example, but its current three-variable equations are those of the Lorenz-63 system. Likewise, `nurged` is an existing filename/hook spelling for the unassimilated comparison trajectory. Preserve these API spellings when running code and use accurate scientific terminology when writing about them.
