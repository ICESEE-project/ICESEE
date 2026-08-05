# Running friction experiments

Run each of the three friction-treatment configurations from the same directory:

```bash
mpirun -np 8 python run_da_issm.py \
  --Nens=40 --model_nprocs=1 \
  -F friction_experiments/seaward_gl_prior/param_ibf.yaml
```

```bash
mpirun -np 8 python run_da_issm.py \
  --Nens=40 --model_nprocs=1 \
  -F friction_experiments/seaward_gl_prior/param_wbf.yaml
```

```bash
mpirun -np 8 python run_da_issm.py \
  --Nens=40 --model_nprocs=1 \
  -F friction_experiments/seaward_gl_prior/param_ebf.yaml
```

The first run uses ISSM's inversion-based friction update (IBF), the
second retains a fixed, uncorrected friction field (WBF), and the third
updates friction through EnKF cross-covariances only (EBF).

## Visualize the comparison

```bash
plot_filter_comparison.m
```