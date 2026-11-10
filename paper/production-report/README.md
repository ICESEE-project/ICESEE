# ICESEE technical report and preprint source

This report documents distributed execution, recoverable histories, state-column timing and production operations at ICESEE development commit `3748e416d3a8120cc860119bd7224c6f9ad41ff7`. It is a new technical-development narrative with a student operations component, distinct from the original model-description paper. It does not assign a new software release or claim a new filter merely because computation is streamed.

The [rendered technical report](ICESEE-technical-report.pdf) is provided alongside its editable `main.tex` source. It includes its bibliography and vector diagrams, with no downloaded figure or publisher-class dependency. Build with a current TeX Live installation:

```bash
./build.sh
```

The result is `build/main.pdf`. Run twice to resolve the table of contents, citations and cross-references. The handover's rendered PDF was checked visually; building successfully alone does not establish layout quality.

Read the [student/developer handbook](../../docs/handbook/README.md) for detailed configuration, API, deployment and operating procedures. `evidence-summary.json` and `measured-results.csv` retain the public numerical summary; the handover evidence bundle contains exact source records, hashes, scoped claims, logs and an audit script. Raw PACE histories and input data are not included in this source directory.

## Scope of evidence

The five-year PIG simulation steps completed in all three modes, but the enclosing simulation and verification-only allocations timed out during output scans. Remaining closed-output checks passed later with a bounded slab reader. Unequal process counts and concurrent execution mean the reported elapsed times are operational measurements, not a scaling curve. Different native/generic initial perturbation policies prevent a three-mode scientific-equivalence claim.

The new documentation checks exercised the exact source with 184 focused tests, 12 real local MPI tests and the existing modes-0/1/2 three-variable Lorenz integration. The path is named `lorenz96`, while its actual equations are Lorenz-63. New container builds, full Icepack production runs and PACE allocations were not performed for this handover.

## Before public submission

1. Compare the novelty statement to the final revised version of the earlier paper, including its August 2026 methodology additions.
2. Reconcile package version, citation/archive metadata and license metadata with the intended release; the report currently identifies an exact commit rather than inventing a new version.
3. Finalize actual authors, affiliations, contributions, funding and acknowledgements.
4. Archive the exact code, environment definitions, configuration, inputs or access restrictions, output shards, validators and figure-generation code under a persistent identifier. An older ICESEE DOI does not archive these developments.
5. Add controlled scaling and scientific evaluation if making those claims. The report states concrete designs and acceptance limits; unperformed experiments remain future work.
6. Recheck venue requirements, references, equations, final figures and the code/data availability statement immediately before submission.

No arXiv or journal submission is made by this documentation update.
