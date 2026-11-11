#!/usr/bin/env python3
"""Rebuild the Lorenz teaching trace from its retained numerical CSV."""
from pathlib import Path
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def main():
    here = Path(__file__).resolve().parent
    with (here/'lorenz-traces.csv').open() as stream:
        rows = list(csv.DictReader(stream))
    fig, axes = plt.subplots(3, 1, figsize=(7.1, 6.8), sharex=True,
                             constrained_layout=True)
    for axis, variable in zip(axes, ('x', 'y', 'z')):
        series = [r for r in rows if r['variable'] == variable]
        t = np.array([float(r['time']) for r in series])
        for member in range(8):
            axis.plot(t, [float(r[f'member_{member}']) for r in series],
                      color='#81b5d7', lw=.7, alpha=.6,
                      label='Members' if member == 0 else None)
        axis.plot(t, [float(r['truth']) for r in series], color='#222222',
                  lw=1.6, label='Synthetic truth')
        axis.plot(t, [float(r['background']) for r in series], color='#b36d15',
                  lw=1.2, ls='--', label='Background')
        axis.plot(t, [float(r['mean']) for r in series], color='#00659b',
                  lw=1.6, label='Ensemble mean')
        observed = [r for r in series if r['observation']]
        axis.scatter([float(r['time']) for r in observed],
                     [float(r['observation']) for r in observed],
                     marker='o', s=22, facecolor='white', edgecolor='#a02945',
                     zorder=5, label='Observations')
        axis.set_ylabel(variable)
        axis.grid(alpha=.2)
        axis.spines[['top', 'right']].set_visible(False)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='outside upper center', ncol=3, frameon=False)
    axes[-1].set_xlabel('Model time')
    out = here/'figures'; out.mkdir(exist_ok=True)
    fig.savefig(out/'lorenz-traces.pdf')
    fig.savefig(out/'lorenz-traces.png', dpi=160)
    plt.close(fig)


if __name__ == '__main__':
    main()
