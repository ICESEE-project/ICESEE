# Reading the numerical methods in context

The references below provide a route from ensemble estimation to the ice-sheet applications that motivate ICESEE. Each answers a different question. A method paper explains a statistical construction; a model paper explains the forecast model and its assumptions. The behavior of a particular ICESEE run additionally depends on its source revision, selected analysis path, configuration, inputs and numerical environment.

The [numerical chapter](numerics.md) gives the implemented factor, controls and time convention. The readings provide context for those choices without replacing the source-level definition. In particular, an ensemble-space transform is not necessarily an ETKF or LETKF, and a compact covariance calculation does not resolve finite-ensemble sampling limitations.

## Estimation and ensemble analysis

| Reference | Question it helps answer | Connection to ICESEE |
| --- | --- | --- |
| [Carrassi et al. (2018), Data assimilation in the geosciences](https://doi.org/10.1002/wcc.535), *WIREs Climate Change*, 9(5), e535. | How do models, observations and uncertainty enter an evolving state estimate? | Establishes the estimation vocabulary used to distinguish initial, model and observation uncertainty. It is a useful starting point before following configuration keys. |
| [Evensen (2003), The Ensemble Kalman Filter](https://doi.org/10.1007/s10236-003-0036-9), *Ocean Dynamics*, 53, 343–367. | How can ensemble statistics support a practical analysis without evolving a full forecast covariance? | Provides numerical and implementation context. ICESEE's exact factor construction, truncation and streamed products still need the definitions in the numerical chapter. Disk-resident members and ensemble formulations have precedents. |
| [Burgers et al. (1998), Analysis Scheme in the Ensemble Kalman Filter](https://doi.org/10.1175/1520-0493%281998%29126%3C1719%3AASITEK%3E2.0.CO%3B2), *Monthly Weather Review*, 126(6), 1719–1724. | What role do observation perturbations play in a stochastic ensemble analysis? | Helps distinguish explicit observation-error factors from ICESEE's forecast-dependent legacy factor. The presence of a configured observation standard deviation alone does not make those policies equivalent. |
| [Tippett et al. (2003), Ensemble Square Root Filters](https://doi.org/10.1175/1520-0493%282003%29131%3C1485%3AESRF%3E2.0.CO%3B2), *Monthly Weather Review*, 131(7), 1485–1490. | Why can several deterministic transforms represent the same analysis covariance? | Gives context for the separate deterministic filter choices. Distributing the existing stochastic calculation changes its execution, rather than selecting one of these algorithms. |
| [Sakov and Oke (2008), A deterministic formulation of the ensemble Kalman filter](https://doi.org/10.1111/j.1600-0870.2007.00299.x), *Tellus A*, 60(2), 361–371. | What approximation underlies the DEnKF anomaly update? | Explains the small-correction, half-gain construction. Its approximation should be described explicitly when comparing it with a square-root or stochastic analysis. |

## Sampling error and local analysis

| Reference | Question it helps answer | Connection to ICESEE |
| --- | --- | --- |
| [Gaspari and Cohn (1999), Construction of correlation functions in two and three dimensions](https://doi.org/10.1002/qj.49712555417), *Quarterly Journal of the Royal Meteorological Society*, 125(554), 723–757. | What mathematical properties underlie commonly used compact correlation tapers? | Provides the correlation-function context for localization. A taper's mathematical form and the scientifically appropriate localization radius are different choices. |
| [Hamill et al. (2001), Distance-Dependent Filtering of Background Error Covariance Estimates in an Ensemble Kalman Filter](https://doi.org/10.1175/1520-0493%282001%29129%3C2776%3ADDFOBE%3E2.0.CO%3B2), *Monthly Weather Review*, 129(11), 2776–2790. | Why do weak, distant sample correlations become problematic in a finite ensemble? | Motivates assessing localization and ensemble size even when analysis avoids dense covariance storage. It does not provide one radius that is appropriate for every ice-sheet observation network. |
| [Hunt et al. (2007), Efficient data assimilation for spatiotemporal chaos](https://doi.org/10.1016/j.physd.2006.11.008), *Physica D*, 230(1–2), 112–126. | What defines a local ensemble transform Kalman filter? | Supplies a specific deterministic local algorithm for comparison. ICESEE's grouped local stochastic transforms retain their own factor and innovation policy; execution Mode 3 is not an LETKF selector. |

## Ice-sheet estimation and forecast models

| Reference | Question it helps answer | Connection to ICESEE |
| --- | --- | --- |
| [Bonan et al. (2014), An ETKF approach for initial state and parameter estimation in ice sheet modelling](https://doi.org/10.5194/npg-21-569-2014), *Nonlinear Processes in Geophysics*, 21(2), 569–582. | How can ensemble methods address uncertain initial conditions and parameters in a flowline model? | Provides a shallow-ice twin-experiment example and a distinct ETKF formulation. Its experiment is prior context, rather than evidence for the current PIG runtime or physical estimates. |
| [Gillet-Chaulet (2020), Assimilation of surface observations in a transient marine ice sheet model using an ensemble Kalman filter](https://doi.org/10.5194/tc-14-811-2020), *The Cryosphere*, 14(3), 811–832. | How can transient surface observations constrain a marine flowline state and basal conditions? | Connects observation type and model evolution to state and parameter estimation. The idealized geometry and estimation assumptions matter when designing a new experiment. |
| [Choi et al. (2025), Estimation of the state and parameters in ice sheet model using an ensemble Kalman filter and Observing System Simulation Experiments](https://doi.org/10.5194/tc-19-5423-2025), *The Cryosphere*, 19, 5423–5444. | How does observing-system design affect estimation in a two-dimensional setting? | Provides context for twin and observing-system experiments, including sensitivity to observation density and prescribed uncertainty. More observations are not a universal guarantee of greater skill. |
| [Shapero et al. (2021), icepack: a new glacier flow modeling package in Python, version 1.0](https://doi.org/10.5194/gmd-14-4593-2021), *Geoscientific Model Development*, 14(7), 4593–4616. | Which glacier-flow models and variational abstractions underpin the Icepack forecast backend? | Separates the model's description from ICESEE's assimilation and runtime layer. The model-paper version does not identify the library versions used by a later production run. |

Together, these readings support a sequence of questions: what uncertainty is represented, how observations update it, how localization changes that calculation, and which physical quantities the chosen model and observation network can constrain. Numerical equivalence between execution modes answers a computational question within that sequence. Physical accuracy, parameter identifiability and uncertainty calibration require their own experiment design and evidence.

## Parallel computation and data layout

| Reference | Question it helps answer | Connection to ICESEE |
| --- | --- | --- |
| [MPI Forum (2025), MPI standard version 5.0](https://www.mpi-forum.org/docs/mpi-5.0/mpi50-report.pdf) | Which processes participate in a communicator and its collective operations? | Clarifies the distinction between model groups and ensemble analysis participants. The standard version is distinct from the installed MPI implementation version. |
| [Amdahl (1967), Validity of the single processor approach](https://doi.org/10.1145/1465482.1465560) | What limits speedup for a fixed workload? | Provides context for root analysis, communication, synchronization and I/O costs that do not decrease with independent forecast work. |
| [Gustafson (1988), Reevaluating Amdahl's law](https://doi.org/10.1145/42411.42415) | How does speedup interpretation change when the workload grows? | Motivates an explicit definition of ensemble, spatial and observation scaling. It does not turn unequal-resource operational runs into a measured weak-scaling study. |
| [h5py documentation, Virtual Datasets](https://docs.h5py.org/en/3.15.0/vds.html) | How can several source datasets appear as one logical array? | Explains why the Mode-2 history view needs its backing files. This documentation version does not identify the h5py version used by the recorded experiments. |

The [parallelization chapter](parallelization-modes.md) connects these concepts to the actual ownership, collectives and storage used by each execution mode.
