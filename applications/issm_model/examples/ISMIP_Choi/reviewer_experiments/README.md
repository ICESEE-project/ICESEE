# Reviewer friction-recovery experiments

All profiles inherit the same ensemble initialization, physical forcing, and
stable EnKF formulation from `../params.yaml`. Surface elevation and horizontal
velocity are observed every 2 years through year 24, and thickness is inferred
rather than observed directly. The fixed-friction and EnKF-only profiles use
one sparse bed snapshot at year 2. The revised hybrid is also the dense-bed
sensitivity: grounded-bed surveys at years 2, 8, 14, 20, and 24 on cross-flow
tracks spaced 10 km apart. Each 50-year run includes a 26-year
observation-free forecast.

1. `wrong_friction_fixed.yaml`: fixed wrong-friction control. The EnKF updates
   the state and bed, while friction has zero process noise and is restored to
   its forecast value after every analysis.
2. `friction_enkf_only.yaml`: augmented-state EnKF recovery of bed and friction,
   with 6 km and 4 km localization radii, respectively.
3. `friction_inversion_hybrid.yaml`: EnKF recovery of the state and bed followed
   by member-wise ISSM velocity inversion for friction. The revised hybrid run
   writes to `_reviewer_friction_inversion_hybrid_v3`. The original sparse-bed
   run and v2 dense-bed run remain available for direct leakage, observation-
   density, and regularization comparisons.
4. `friction_inversion_hybrid_low_prior_tuned.yaml`: tuned low-prior hybrid.
   It preserves the successful 85% thickness and -80 m grounded-bed initial
   conditions, assimilates annual surface/velocity observations through year
   30, and retains sparse bed surveys at years 2, 8, 14, 20, and 24. Friction
   remains fixed during the geometry-only analyses before year 6, after which
   member-wise inversion uses bounds of 1500--4500. The localized bed increment
   cap is reduced from 30 m to 20 m to suppress the observed survey-corridor
   overshoot patches.

Run from the `ISMIP_Choi` directory, replacing the MPI layout as needed:

```bash
# The DA initializer reads ./bed_kriging_results.h5; regenerate it before a
# run whenever the kriging formulation, observation geometry, or Nens changes.
python generate_bed_kringing.py \
  --data-path _reviewer_friction_inversion_hybrid \
  --Ne 40 --stride-km 10 --track-half-width-km 2.5 \
  --background-length-km 40 \
  --output-file bed_kriging_results.h5

mpiexec -n 60 python run_da_issm.py -F reviewer_experiments/wrong_friction_fixed.yaml
mpiexec -n 60 python run_da_issm.py -F reviewer_experiments/friction_enkf_only.yaml
mpiexec -n 60 python run_da_issm.py -F reviewer_experiments/friction_inversion_hybrid.yaml

# Recommended tuned hybrid sensitivity. This creates its own output directory;
# do not resume it from an older profile with a different observation schedule.
mpirun -np 8 python run_da_issm.py --Nens=40 --model_nprocs=1 \
  -F reviewer_experiments/friction_inversion_hybrid_low_prior_tuned.yaml
```

Report RMSE for thickness, velocity, bed, and grounded friction (also grounded
excluding the grounding-line band), plus centerline grounding-line error. Use
the identical masks for all profiles. Compare both the final analysis at year 24
and the free-forecast endpoint at year 50.
