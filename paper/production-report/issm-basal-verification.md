# Compact ISSM backend verification

We completed the fixed-friction, joint-EnKF and hybrid teaching profiles with four members, ten forecast cycles and two analyses, using ISSM 4.24 and MATLAB R2025b. The tested runtime source is `02ac5496d3d6cacf0ded737cd821e3e6f43c9d42`; exact source, configuration, input and output hashes are retained in [the result manifest](issm-basal-results.json). [The measurements CSV](issm-basal-measurements.csv) identifies the artifact behind each value.

| Profile | Date | Elapsed seconds | Member inverse calls |
|---|---|---:|---:|
| Fixed friction | 10 November 2026 | 118.751 | 0 |
| Joint EnKF | 11 November 2026 | 109.871 | 0 |
| Hybrid | 10 November 2026 | 153.071 | 4 |

Each case used Mode 1 with one Python rank. Every outgoing ISSM command used one internal process under the explicit `ICESEE_PERFORMANCE_TEST: true` resource-retention policy. That flag invokes an existing resource branch; it is not a benchmark claim. Historical uncapped attempts used different internal process counts and remain clearly labeled in the manifest.

Every accepted member history has shape `(20082, 4, 11)` and its mean has shape `(20082, 11)`. The histories and means are finite, every stored column contains data, and stored means match recomputed means exactly. Truth and background have nonzero initial and terminal columns matching MATLAB exports `*_0` and `*_10` exactly. Analysis records align years 1 and 2 with state/observation columns 5 and 10. The EnKF case additionally passed a recursive scan of 31 HDF5 files and 233 datasets. Missing-observation NaNs in `hu_obs` and transposed `R` follow matching masks; active error values are finite and nonnegative.

The hybrid case called all four member inverses at year two. Each returned 3347 finite coefficient values, changed its packed coefficient array, and matched the stored final coefficient exactly. All four optimizers stopped on `dxmin` during line search with 29 printed objective/gradient rows. We retain the rounded printed objective/gradient values without calling the last trial the accepted final objective. These observations do not establish gradient-tolerance convergence, objective improvement, calibrated uncertainty or physical skill. BAMG nodal correspondence and the documented inverse geometry/velocity handoffs remain separate scientific validation questions.

The existing one-year prior spin-up was retained. State-column export and delayed-inverse activation were corrected without tuning forecast physics or the inverse objective. These cases establish the stated local execution, storage and handoff contracts; they do not provide multi-rank ICESEE scaling or coupled-restart evidence. Elapsed times include startup and completion checks and are not a controlled performance comparison.
