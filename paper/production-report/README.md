# ICESEE technical report and preprint source

This report explains the four parallelization modes, ensemble formulation, model coupling, related literature, recoverable histories and observation timing at ICESEE development commit `3748e416d3a8120cc860119bd7224c6f9ad41ff7`. It is a new technical-development narrative with a student operations component, distinct from the original model-description paper. It does not assign a new software release or claim a new filter merely because computation is streamed.

The [rendered technical report](ICESEE-technical-report.pdf) is provided alongside its editable `main.tex` source. Its bibliography is embedded in the LaTeX, with verified BibTeX metadata also provided in `references.bib`. It includes vector diagrams, with no downloaded figure or publisher-class dependency. Build with a current TeX Live installation:

```bash
./build.sh
```

The result is `build/main.pdf`. The script performs three passes to resolve the contents, author-year citations and cross-references. The delivered PDF is visually reviewed after each substantive revision. The scientific narrative follows the author's first-person plural style, with primary-source citations at the points where the underlying ideas are introduced.

The [parallelization chapter](../../docs/handbook/parallelization-modes.md) traces each mode through initialization, forecasts, observation and state-row work, collectives, storage and restart. The [annotated literature](../../docs/handbook/literature.md) provides a reading route through the methods and software context. [The endpoint experiment](experiments/README.md) reproduces the small two-rank Modes 1–3 checks. The [sample-run route](../../docs/handbook/sample-runs.md) includes the new coarse Icepack case and compact ISSM exercises.

Detailed [inversion and physical inference](../../docs/handbook/inversion-and-inference.md), [current developments](../../docs/handbook/current-developments.md), [Icepack coupling](../../docs/handbook/coupling-icepack.md), [ISSM coupling](../../docs/handbook/coupling-issm.md), [new-application onboarding](../../docs/handbook/application-onboarding.md) and [adapter templates](../../docs/handbook/templates/README.md) accompany the report. They distinguish source-inspected paths from the numerical cases actually exercised.

Read the [student/developer handbook](../../docs/handbook/README.md) for detailed configuration, API, deployment and operating procedures. `evidence-summary.json` and `measured-results.csv` retain the public numerical summary; the handover evidence bundle contains exact source records, hashes, scoped claims, logs and an audit script. Raw PACE histories and input data are not included in this source directory.

## Scope of evidence

The five-year PIG simulation steps completed in all three modes, but the enclosing simulation and verification-only allocations timed out during output scans. Remaining closed-output checks passed later with a bounded slab reader. Unequal process counts and concurrent execution mean the reported elapsed times are operational measurements, not a scaling curve. Different native/generic initial perturbation policies prevent a three-mode scientific-equivalence claim.

The new documentation checks exercised the exact source with 184 focused tests, 12 real local MPI tests and the existing modes-0/1/2 three-variable Lorenz integration. The path is named `lorenz96`, while its actual equations are Lorenz-63. The new coarse synthetic Icepack case was executed in Modes 0–3 (one rank for Mode 0, two for Modes 1–3); all completed ten forecasts and two analyses with finite `(140,4,11)` histories. Modes 1/2 differ by at most `2.01e-10`. Native output uses state snapshots and sharded history rather than the generic history export. New container builds, full Icepack production runs and PACE allocations were not performed for this handover.

The [ISSM basal-friction record](issm-basal-results.json), [measurement table](issm-basal-measurements.csv) and [verification guide](issm-basal-verification.md) cover actual fixed, joint-EnKF and hybrid treatments in Mode 1, each with one Python rank and one internal solver process. All produced finite `(20082,4,11)` histories, consistent means and aligned initial/terminal reference states. Only hybrid invoked the inverse: four member calls at year two, with exact coefficient handoff. All four optimizers stopped on `dxmin` during line search; completion is not evidence of gradient-tolerance convergence, objective improvement or parameter skill. Choi-specific inverse/native capabilities remain source-inspected.

## Before public submission

1. Compare the novelty statement to the final revised version of the earlier paper, including its August 2026 methodology additions.
2. Reconcile package version, citation/archive metadata and license metadata with the intended release; the report currently identifies an exact commit rather than inventing a new version.
3. Finalize actual authors, affiliations, contributions, funding and acknowledgements.
4. Archive the exact code, environment definitions, configuration, inputs or access restrictions, output shards, validators and figure-generation code under a persistent identifier. An older ICESEE DOI does not archive these developments.
5. Add controlled scaling and scientific evaluation if making those claims. The report states concrete designs and acceptance limits; unperformed experiments remain future work.
6. Recheck venue requirements, references, equations, final figures and the code/data availability statement immediately before submission.

No arXiv or journal submission is made by this documentation update.
