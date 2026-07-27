# Rebuttal WBF / EBF / IBF experiments

`param.yaml` is the single, self-contained source of common settings. It does
not extend another configuration. Each run extends it directly:

- `param_ibf.yaml`: EnKF state/bed update followed by ISSM friction inversion;
- `param_wbf.yaml`: EnKF state/bed update with fixed wrong friction;
- `param_ebf.yaml`: augmented-state EnKF friction recovery without inversion.

All runs span 50 years, assimilate annual surface/velocity observations from
year 2 through year 40, and use sparse bed surveys at years 2, 8, 14, 20, 24,
30, and 36. The wrong initial bed is constructed across the full mesh; bed
observations and EnKF bed increments remain restricted to grounded ice.

Before an expensive ensemble run, generate only the initial states and inspect
the full-domain signed bed difference:

```bash
mpirun -np 1 python run_da_issm.py --Nens=1 --model_nprocs=1 \
  -F rebutal_experiments/param_ic_check.yaml

python rebutal_experiments/check_initial_prior.py \
  _modelrun_datasets_rebuttal_ic_check
```

Then run from the `ISMIP_Choi` directory:

```bash
mpirun -np 8 python run_da_issm.py --Nens=40 --model_nprocs=1 \
  -F rebutal_experiments/param_ibf.yaml

REBUTTAL_MPI_RANKS=8 REBUTTAL_NENS=40 REBUTTAL_MODEL_NPROCS=1 \
  bash rebutal_experiments/run_experiments.sh
```

The corrected profiles write to new `*_fullbed_v2` directories. They do not
resume or overwrite the earlier runs whose floating bed was preserved from
the reference state.
