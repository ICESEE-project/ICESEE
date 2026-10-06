# Forty-year basal-friction method comparison

This directory defines a matched comparison of:

- **WBF:** a prescribed, deliberately incorrect basal-friction field;
- **EBF:** joint EnKF estimation of basal friction, geometry, and velocity;
- **IBF:** member-wise ISSM friction inversion coupled to EnKF state and bed
  estimation.

The experiment is independent of the principal long-duration IBF experiment.
All three cases use the same 40-member initial ensemble, four-year
pre-assimilation stabilization, `dt = 0.1 yr`, forcing, truth, and observation
realizations. Surface elevation and horizontal velocity are observed annually
from years 2 through 24. Bed elevation is observed at years 2, 6, 10, 14, 18,
and 22, with a terminal survey at year 24. Assimilation then stops, leaving a
16-year free forecast through year 40.

Only the friction treatment changes among the three child profiles.

This separate numerical protocol was required for a fair comparison. In
preliminary matched WBF and EBF runs started directly with the `dt = 0.2 yr`
step used by the principal 100-year IBF experiment, the diagnostic
stress-balance velocity underwent a sharp initialization transient. Although
the EBF analysis changed basal friction, those increments did not reliably
alter the following forecast trajectory. The IBF workflow made this defect
less apparent because its member-wise inversion re-solves the stress balance
after an analysis. A coupling preflight showed that a common four-year
stabilization at `dt = 0.1 yr`, followed by forecasts at the same step, allowed
EBF friction increments to propagate into the dynamic state. The stabilization
and smaller step are therefore applied identically to WBF, EBF, and IBF; they
are numerical controls, not additional information supplied to IBF. We call
the observed behavior an initialization transient or velocity collapse, not
numerical diffusion, because the tests diagnose the symptom and coupling
failure but do not establish a diffusion mechanism.

## Commands

Run these commands from `applications/issm_model/examples/ISMIP_Choi`:

```bash
mpirun -np 8 python run_da_issm.py --Nens=40 --model_nprocs=1 \
  -F rebutal_experiments/method_comparison_40yr/param_wbf.yaml

mpirun -np 8 python run_da_issm.py --Nens=40 --model_nprocs=1 \
  -F rebutal_experiments/method_comparison_40yr/param_ebf.yaml

mpirun -np 8 python run_da_issm.py --Nens=40 --model_nprocs=1 \
  -F rebutal_experiments/method_comparison_40yr/param_ibf.yaml
```

## Suggested manuscript introduction

> To isolate the effect of basal-friction treatment, we performed an
> additional controlled 40-year experiment using identical initial ensembles,
> forcing, observations, and EnKF settings. Preliminary WBF and EBF runs
> initialized directly at the 0.2-year step used in the principal experiment
> exhibited a sharp diagnostic-velocity initialization transient, and EBF
> friction increments did not reliably propagate into the following forecast.
> All three comparison cases therefore use the same four-year stabilization
> and subsequent forecasts at a 0.1-year step; a coupling preflight verified
> propagation of EBF coefficient increments under this common protocol.
> Surface elevation and horizontal
> velocity were assimilated annually from years 2–24. Bed-elevation surveys
> were assimilated at years 2, 6, 10, 14, 18, and 22, with a terminal survey
> at year 24. Assimilation was then discontinued, and all experiments were
> propagated freely through year 40. The three cases differed only in their
> treatment of basal friction: prescribed wrong basal friction (WBF), joint
> EnKF friction estimation (EBF), and friction recovered by member-wise
> inversion within the hybrid EnKF framework (IBF). This targeted experiment
> is used only for the method-comparison figure; the principal IBF results use
> the long-duration configuration described separately.
