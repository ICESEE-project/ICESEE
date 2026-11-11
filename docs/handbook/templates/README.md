# Application scaffolds

These files support the [onboarding guide](../application-onboarding.md). They
document a starting interface, rather than supplying an executable model.
The `.py.in` bodies raise `NotImplementedError` where physics is required,
and none registers a model or imports application configuration at module
scope. The YAML family is intentionally absent from the current registry.

| File | Intended use |
| --- | --- |
| [legacy-hooks.py.in](legacy-hooks.py.in) | Core generic signatures based on the PIG hooks, with explicit clock and file-backed return contracts. |
| [application-driver.py.in](application-driver.py.in) | Entrypoint structure with configuration deferred to `main` and an unimplemented setup boundary. |
| [native-adapter.py.in](native-adapter.py.in) | Persistent native callback signatures; optional capabilities remain absent. |
| [params.yaml.in](params.yaml.in) | Existing YAML section/key structure with illustrative dimensions and unresolved physical meaning. |
| [capabilities.md.in](capabilities.md.in) | Evidence record separating implementation from tested behavior for every advertised mode. |

A new example first chooses its real family/directory names and records its
scientific/state/observation specification. The developer copies only the
needed scaffolds into that example, implements the actual model setup and
callbacks, and checks how the selected runners consume their returns. The
generic hooks and a native adapter are distinct boundaries; copying both does
not connect them automatically. Native registration follows an implemented,
tested lifecycle, rather than serving as a substitute for it.

`example_name` is metadata. Legacy discovery uses the current directory's
basename and the family key to construct the hook module path. The new
directory must therefore be a valid Python identifier and the entrypoint must
establish the intended directory before dispatch. Optional hooks are added
only when their active call sites and capabilities are understood.

The numerical YAML values describe an interface example and are not a
scientific recommendation. In particular, timestep conventions differ between
existing Icepack examples, and the new driver must establish one `dt`, `nt`
and time grid. Fresh configuration loading can remove an existing output
directory. Importing the scaffolds for inspection does not load that module;
executing a completed driver requires a new reviewed run path or an explicitly
supported resume route.

An AST/YAML syntax check establishes only that the scaffolds parse. Model-only,
observation, numerical-parity, MPI, recovery and physical-validation gates in
the onboarding guide are still required. No test result for these incomplete
templates is a new coupled-model result.
