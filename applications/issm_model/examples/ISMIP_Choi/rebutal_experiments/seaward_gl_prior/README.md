# Seaward-grounding-line prior experiment

This robustness experiment uses the accepted rebuttal configuration but
starts the model grounding line on the seaward side of the truth. It is the
opposite-sign grounding-line challenge to the landward prior shown in the
existing long simulation.

The displacement is created physically: a compact, smooth thickness
perturbation is added across both sides of the grounding zone. It first
repairs any local retreat caused by the heterogeneous background prior and
then grounds a controlled strip of shelf. The bed, surface, base, thickness,
mask, and hydrostatic floating geometry remain internally consistent. No
plotted contour is translated. Distance is measured along flow from the
locally upstream GL front, so the lateral arms of the U-shaped GL cannot
spread the perturbation across the whole shelf.

First generate only the initial state:

```bash
mpirun -np 8 python run_da_issm.py \
  --Nens=40 --model_nprocs=1 \
  -F rebutal_experiments/seaward_gl_prior/param_ic_check.yaml
```

Then audit its direction and magnitude:

```bash
python rebutal_experiments/check_initial_prior.py \
  _modelrun_datasets_rebuttal_ic_seaward_gl_v1 \
  --expect-seaward-gl
```

Do not launch the full IBF run until the audit confirms:

1. `prior-truth` centerline GL offset is positive but no more than 40 km;
2. floating-to-grounded vertices exceed grounded-to-floating vertices;
3. the geometry fields remain smooth and the bed-error statistics remain
   comparable with the accepted baseline.

After approval, the full hybrid run is:

```bash
mpirun -np 8 python run_da_issm.py \
  --Nens=40 --model_nprocs=1 \
  -F rebutal_experiments/seaward_gl_prior/param_ibf.yaml
```
