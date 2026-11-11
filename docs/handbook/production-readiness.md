# Production handover status and acceptance checklist

This checklist records what was inspected or exercised during the 10 November 2026 documentation handover and what still requires a local decision. It authorizes no new scientific claims and does not certify every ICESEE application as production ready. Use the pinned source, [operations procedure](operations.md) and application-specific validation record together. A case is accepted by its recorded configuration and evidence, not by a branch name or a completed checklist for another case.

## Inspected and exercised

| Item | Status at this handover | Practical meaning |
|---|---|---|
| Exact source baseline | **Inspected:** `features_debug` at `3748e416d3a8120cc860119bd7224c6f9ad41ff7` | `main` differs; installers and environment activation can import another checkout. Check the actual import path and SHA in every batch environment. |
| Student smoke procedure | **Exercised by the documentation workflow:** serial 8-member, 50-cycle case; four analyses at 0.2, 0.4, 0.6, 0.8; complete `(3,8,51)` member history | This establishes that the supplied small tutorial is executable and its checked outputs are consistent. It does not establish that students have independently practiced it or that distributed glacier physics is valid. |
| Local focused and MPI regression checks | **Exercised in the handover environment:** 184 focused tests and 12 real-MPI tests passed, with no skips in those selections | These are selected contract/topology/recovery gates. The exact selections and environment are in the handover evidence package; they are not all possible model/backend tests. |
| Existing Lorenz integration | **Exercised:** Modes 0/1/2, 30 members, 1,000 cycles with the repository comparator | A lightweight full-driver check. Its historical path says Lorenz96 but implements Lorenz-63 equations. |
| Terminal-observation Lorenz fixture | **Executed locally:** two-rank Modes 1/2/3, eight members, 50 forecasts and five analyses including column 50 | Finite `(3,8,51)` histories and consistent means; Modes 1/2 complete-member comparison passes `rtol=atol=1e-7`. Mode 3 has completion/integrity evidence; see [sample runs](sample-runs.md). |
| Coarse Icepack teaching fixture | **Executed locally:** Mode 0 on one Python rank and Modes 1–3 on two ranks; four members, ten forecasts and two analyses | Finite `(140,4,11)` packed state series and positive thickness; Modes 1/2 maximum absolute difference `2.007709554163739e-10`. Native snapshots/shards have a distinct output contract. This is neither multi-rank spatial validation nor a skill/scaling result; see [sample runs](sample-runs.md). |
| Compact basal-friction ISSM treatments | **Executed:** fixed, joint-EnKF and hybrid in Mode 1, one Python rank and one internal solver process | Finite `(20082,4,11)` histories, consistent means and two aligned analyses in every treatment. Hybrid made four inverse calls with exact coefficient handoff, all stopping on line-search `dxmin`; optimizer convergence/skill is not established. Choi remains a separate source-inspected case. See [sample runs](sample-runs.md). |
| Recorded five-year PIG calculation | **Retained evidence:** 32 members, 100 cycles, three analysis times; all three simulation steps exited zero | Parent/verification allocations timed out during later reads; remaining native output checks subsequently passed a bounded standalone scan. Record these statuses separately. |
| Recorded PIG restart | **Retained evidence:** short checkpoint replay crossed analysis events and agreed within the declared comparison tolerance | It does not establish all interruption, mesh repartitioning, cross-environment or other-adapter resume combinations. |
| Full/partial history procedures | **Documented and contract-tested in selected paths** | Preserve backing shards/manifests and use closed files. Full history cannot be recovered after rolling retention has removed it. |
| Updated Spack/container deployment instructions | **Repository inspection and procedure documentation** | No new full Spack build, container build or host/container MPI qualification was performed in this handover. |

The manuscript and evidence registry provide exact run identifiers, tolerances, selected logs and hashes. Do not promote these statements into a new strong/weak-scaling result or three-mode member-wise scientific equivalence: the production timing record uses unequal resources, concurrent execution and different native initialization.

## Required before a student team depends on a new production case

| Gate | Current disposition | Go/no-go decision |
|---|---|---|
| Local primary/backup operator and scientific/software reviewer roles | **Open:** procedures define roles; no people have been assigned by this document | Record accepted local responsibility and a reachable escalation route. Do not depend on the absent developer as the only operator. |
| Independent student practice | **Open:** automated smoke/restart evidence does not demonstrate human rehearsal | Two operators independently reproduce the small case, locate its evidence and inspect/restore the target-model checkpoint using approved resources. |
| Exact target inputs and access | **Case-dependent:** PIG paths/inputs are documented; a new portable full input archive was not deposited here | Confirm every rank can read the intended mesh, spin-up, observations and forcing, with checksums and backup access. Stop a case with missing/ambiguous inputs. |
| Installed source and numerical stack | **Case-dependent:** source/Spack/container pins are distinct | Functional MPI probe, correct import SHA, correct h5py capability and reduced target-model acceptance must pass in the actual batch environment. |
| Physical acceptance criteria | **Open for each new scientific configuration** | Declare units, admissible field ranges, geometry constraints, innovations/spread and parameter/holdout diagnostics before the expensive run. A finite file is insufficient. |
| Native melt extrema in retained PIG evidence | **Requires scientific interpretation:** transient analyzed values approximately -0.3664 to 20.0441; final values 0 to 20 in configured model units | Do not infer a universally valid melt bound or silently clip to obtain acceptance. The reviewer must decide whether transient values are physically permissible for the specific experiment and document any changed policy. |
| Storage and restart retention | **Case-dependent** | Include raw checkpoints, output copies/VDS sources, means, temporary export and quota/inodes. Preserve restart history and verify the moved archive before deleting originals. |
| Source/environment change during production series | **Not automatically accepted** | Keep the known working baseline and rerun selected parity/restart/target-model gates before promotion. |

These open gates are work for the local operating team and scientific reviewer; they are not requests to interrupt the absent developer for every routine job. Unchanged accepted cases can follow the documented procedures. A case that fails an integrity gate should remain unaccepted while unaffected approved work continues.

## Required before public preprint or journal submission

| Item | Current disposition |
|---|---|
| Release identity and license metadata | **Unresolved:** package metadata `0.1.9`/MIT conflicts with citation metadata `1.0.0` and BSD 2-Clause license text. Preserve source commit identity until maintainers reconcile the records. |
| Persistent archive for these developments | **Not newly deposited:** the earlier paper's archive is not the archive of this later snapshot. Deposit exact new source/configurations, input data or access arrangements, results and figure-generation code. |
| Author participation, order and affiliations | **To be finalized from actual contributions** before submission. The handover draft is not an authorship agreement. |
| New scaling or scientific skill claims | **Not established by handover checks:** require controlled experimental designs and predeclared metrics described in the manuscript. |
| Venue and submission | **Draft stage:** no arXiv/GMD/JAMES submission was made by this documentation handover. |

Keep a dated copy of the completed case-specific acceptance record. During the six months, maintain a run ledger and periodic backup/restart practice as described in operations; this checklist itself does not schedule or monitor future work.
