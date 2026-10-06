#!/usr/bin/env python3
"""2-D MPI parity gate for the mode-3 grouped-local runtime bridge."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from mpi4py import MPI


REPOSITORY = Path(__file__).resolve().parents[2]
if str(REPOSITORY) not in sys.path:
    sys.path.insert(0, str(REPOSITORY))

from src.parallelization.distributed_adapter import (  # noqa: E402
    DistributedAnalysisTargets,
    DistributedObservationCoordinates,
    DistributedStateLayout,
)
from src.parallelization.distributed_local_runtime import (  # noqa: E402
    apply_grouped_local_analysis_local,
)
from src.parallelization.distributed_runtime import (  # noqa: E402
    LocalMemberEnsemble,
    members_for_ensemble_slot,
)
from src.parallelization.distributed_topology import (  # noqa: E402
    create_distributed_topology,
)
from src.utils.localization import compute_X5_from_matrices  # noqa: E402


def _interval(count: int, rank: int, size: int) -> tuple[int, int]:
    quotient, remainder = divmod(int(count), int(size))
    start = rank * quotient + min(rank, remainder)
    return start, start + quotient + (rank < remainder)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ensemble-groups", type=int, default=2)
    parser.add_argument("--spatial-ranks", type=int, default=2)
    parser.add_argument("--nens", type=int, default=4)
    parser.add_argument("--targets", type=int, default=31)
    parser.add_argument("--observations", type=int, default=17)
    parser.add_argument("--radius", type=float, default=0.13)
    parser.add_argument("--chunk-rows", type=int, default=4)
    parser.add_argument("--atol", type=float, default=1.0e-12)
    parser.add_argument("--rtol", type=float, default=1.0e-12)
    args = parser.parse_args()

    topology = create_distributed_topology(
        MPI.COMM_WORLD,
        ensemble_groups=args.ensemble_groups,
        spatial_ranks=args.spatial_ranks,
    )
    if args.nens < args.ensemble_groups:
        parser.error("nens must be at least the number of ensemble groups")

    rng = np.random.default_rng(20260822)
    state = rng.normal(size=(args.targets, args.nens))
    forecast_observations = rng.normal(size=(args.observations, args.nens))
    observed = rng.normal(size=args.observations)
    target_coords = np.linspace(0.0, 1.0, args.targets)[:, None]
    obs_coords = np.linspace(0.025, 0.975, args.observations)[:, None]
    obs_ids = 1000 + 3 * np.arange(args.observations, dtype=np.int64)

    target0, target1 = _interval(
        args.targets, topology.spatial_rank, topology.spatial_ranks
    )
    obs0, obs1 = _interval(
        args.observations, topology.spatial_rank, topology.spatial_ranks
    )
    member_ids = members_for_ensemble_slot(
        args.nens, topology.ensemble_groups, topology.ensemble_slot
    )
    local = LocalMemberEnsemble(
        DistributedStateLayout(args.targets, target0, target1),
        {member: state[target0:target1, member] for member in member_ids},
    )
    target = DistributedAnalysisTargets(
        "state",
        np.arange(target1 - target0, dtype=np.int64),
        np.arange(target0, target1, dtype=np.int64),
        target_coords[target0:target1],
        ("state",),
    )
    metadata = DistributedObservationCoordinates(
        "state", obs_ids[obs0:obs1], obs_coords[obs0:obs1]
    )
    local_observations = {
        member: forecast_observations[obs0:obs1, member]
        for member in member_ids
    }
    updated = apply_grouped_local_analysis_local(
        local,
        target,
        metadata,
        local_observations,
        observed[obs0:obs1],
        radius=args.radius,
        topology=topology,
        icesee_kwargs={"Nens": args.nens},
        error_mode="legacy_prior_anomalies",
        target_chunk_size=args.chunk_rows,
    )

    payload = [
        (member, target0, target1, updated.members[member])
        for member in member_ids
    ]
    gathered = topology.world.gather(payload, root=0)
    passed = None
    if topology.world_rank == 0:
        distributed = np.empty_like(state)
        for rank_payload in gathered:
            for member, start, stop, values in rank_payload:
                distributed[start:stop, member] = values
        reference = state.copy()
        yprime = (
            forecast_observations
            - forecast_observations.mean(axis=1, keepdims=True)
        )
        eta = yprime
        innovations = observed[:, None] - forecast_observations
        for row, coordinate in enumerate(target_coords):
            selected = np.flatnonzero(
                np.linalg.norm(obs_coords - coordinate, axis=1) <= args.radius
            )
            if selected.size:
                order = selected[np.argsort(obs_ids[selected])]
                transform = compute_X5_from_matrices(
                    yprime[order], eta[order], innovations[order], args.nens
                )
                reference[row] = state[row] @ transform
        maximum_error = float(np.max(np.abs(distributed - reference), initial=0.0))
        passed = bool(np.allclose(
            distributed, reference, atol=args.atol, rtol=args.rtol
        ))
        print("Mode-3 2-D MPI grouped-local runtime parity")
        print(
            f"  process grid:       {args.ensemble_groups} x "
            f"{args.spatial_ranks}"
        )
        print(f"  ensemble members:   {args.nens}")
        print(f"  target rows:        {args.targets}")
        print(f"  observation rows:   {args.observations}")
        print(f"  max absolute error: {maximum_error:.12g}")
        print(f"  result:             {'PASS' if passed else 'FAIL'}")
    passed = topology.world.bcast(passed, root=0)
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
