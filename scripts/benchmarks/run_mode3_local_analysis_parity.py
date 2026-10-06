#!/usr/bin/env python3
"""Real-MPI parity gate for mode-3 distributed grouped-local analysis."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from mpi4py import MPI


REPOSITORY = Path(__file__).resolve().parents[2]
if str(REPOSITORY) not in sys.path:
    sys.path.insert(0, str(REPOSITORY))

from src.parallelization.distributed_local_analysis import (  # noqa: E402
    apply_distributed_local_patches,
    iter_distributed_local_patches,
)
from src.utils.localization import compute_X5_from_matrices  # noqa: E402


def _interval(count: int, rank: int, size: int) -> tuple[int, int]:
    quotient, remainder = divmod(int(count), int(size))
    start = rank * quotient + min(rank, remainder)
    return start, start + quotient + (rank < remainder)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--targets", type=int, default=31)
    parser.add_argument("--observations", type=int, default=17)
    parser.add_argument("--nens", type=int, default=6)
    parser.add_argument("--radius", type=float, default=0.13)
    parser.add_argument("--chunk-rows", type=int, default=4)
    parser.add_argument("--atol", type=float, default=1.0e-12)
    parser.add_argument("--rtol", type=float, default=1.0e-12)
    args = parser.parse_args()

    comm = MPI.COMM_WORLD
    rank, size = comm.Get_rank(), comm.Get_size()
    if min(args.targets, args.observations, args.nens, args.chunk_rows) <= 0:
        parser.error("counts and chunk rows must be positive")

    target_coords = np.linspace(0.0, 1.0, args.targets)[:, None]
    obs_coords = np.linspace(0.025, 0.975, args.observations)[:, None]
    observation_ids = 1000 + 3 * np.arange(args.observations, dtype=np.int64)
    rng = np.random.default_rng(20260821)
    yprime = rng.normal(size=(args.observations, args.nens))
    yprime -= yprime.mean(axis=1, keepdims=True)
    eta = rng.normal(scale=0.15, size=yprime.shape)
    eta -= eta.mean(axis=1, keepdims=True)
    innovations = rng.normal(size=yprime.shape)
    forecast = rng.normal(size=(args.targets, args.nens))

    target0, target1 = _interval(args.targets, rank, size)
    obs0, obs1 = _interval(args.observations, rank, size)
    local_forecast = forecast[target0:target1].copy()
    local_analysis = local_forecast.copy()
    for rows, patches in iter_distributed_local_patches(
        target_global_rows=np.arange(target0, target1, dtype=np.int64),
        target_coordinates=target_coords[target0:target1],
        observation_ids=observation_ids[obs0:obs1],
        observation_coordinates=obs_coords[obs0:obs1],
        yprime=yprime[obs0:obs1],
        eta=eta[obs0:obs1],
        innovations=innovations[obs0:obs1],
        radius=args.radius,
        spatial_comm=comm,
        target_chunk_size=args.chunk_rows,
    ):
        apply_distributed_local_patches(
            local_analysis[rows], local_forecast[rows], patches
        )

    gathered = comm.gather((target0, target1, local_analysis), root=0)
    passed = None
    if rank == 0:
        distributed = np.empty_like(forecast)
        for start, stop, values in sorted(gathered):
            distributed[start:stop] = values
        reference = forecast.copy()
        for target_index, coordinate in enumerate(target_coords):
            neighbors = np.flatnonzero(
                np.linalg.norm(obs_coords - coordinate, axis=1) <= args.radius
            )
            if neighbors.size:
                order = neighbors[np.argsort(observation_ids[neighbors])]
                transform = compute_X5_from_matrices(
                    yprime[order], eta[order], innovations[order], args.nens
                )
                reference[target_index] = forecast[target_index] @ transform
        difference = np.abs(distributed - reference)
        maximum_error = float(np.max(difference, initial=0.0))
        passed = bool(np.allclose(
            distributed, reference, rtol=args.rtol, atol=args.atol
        ))
        print("Mode-3 distributed grouped-local stochastic-analysis parity")
        print(f"  spatial ranks:       {size}")
        print(f"  target rows:         {args.targets}")
        print(f"  observation rows:    {args.observations}")
        print(f"  ensemble members:    {args.nens}")
        print(f"  target chunk rows:   {args.chunk_rows}")
        print(f"  max absolute error:  {maximum_error:.12g}")
        print(f"  result:              {'PASS' if passed else 'FAIL'}")
    passed = comm.bcast(passed, root=0)
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
