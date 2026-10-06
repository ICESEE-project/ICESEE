#!/usr/bin/env python3
"""Exercise non-selectable mode-3 forecast and analysis paths under real MPI.

This is an early state-only gate, not an ICESEE execution mode.  It compares
distributed Lorenz forecasts and a stochastic-EnKF analysis transform with
their serial whole-member equivalents. Complete members are reconstructed only
by the explicit parity helper and never by the mode-3 runtime itself.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from mpi4py import MPI


REPOSITORY = Path(__file__).resolve().parents[2]
if str(REPOSITORY) not in sys.path:
    sys.path.insert(0, str(REPOSITORY))

from applications.lorenz_model.lorenz_utils.distributed_adapter import (  # noqa: E402
    LorenzDistributedAdapter,
    lorenz_rk4_step,
)
from src.parallelization.distributed_runtime import (  # noqa: E402
    apply_ensemble_transform_local,
    forecast_local_members,
    initialize_local_members,
    reconstruct_global_member_for_parity,
)
from src.parallelization.distributed_analysis import (  # noqa: E402
    StochasticAnalysisProducts,
    contiguous_observation_layout,
    ensemble_transform_from_products,
    iter_local_ensemble_row_blocks,
    validate_observation_partition,
)
from src.parallelization.distributed_topology import (  # noqa: E402
    create_distributed_topology,
)


def _initial_ensemble(number_of_members: int) -> np.ndarray:
    """Return deterministic, distinct initial states for the parity gate."""

    member = np.arange(number_of_members, dtype=float)
    return np.vstack(
        (
            1.0 + 0.10 * member,
            2.0 - 0.05 * member,
            3.0 + 0.02 * member,
        )
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ensemble-groups", type=int, default=2)
    parser.add_argument("--nens", type=int, default=4)
    parser.add_argument("--steps", type=int, default=20)
    parser.add_argument("--analysis-step", type=int, default=6)
    parser.add_argument("--dt", type=float, default=0.01)
    parser.add_argument("--atol", type=float, default=1.0e-13)
    parser.add_argument("--rtol", type=float, default=1.0e-13)
    args = parser.parse_args()

    world = MPI.COMM_WORLD
    if args.ensemble_groups <= 0 or args.nens <= 0 or args.steps < 0:
        parser.error("ensemble groups and members must be positive; steps cannot be negative")
    if world.Get_size() % args.ensemble_groups:
        parser.error("MPI size must be divisible by --ensemble-groups")

    topology = create_distributed_topology(
        world,
        ensemble_groups=args.ensemble_groups,
    )
    initial = _initial_ensemble(args.nens)
    icesee_kwargs = {
        "Nens": args.nens,
        "nd": 3,
        "u0b": initial[:, 0],
        "distributed_initial_ensemble": initial,
        "sigma_96": 10.0,
        "beta_96": 8.0 / 3.0,
        "rho_96": 28.0,
        "dt": args.dt,
    }
    adapter = LorenzDistributedAdapter()
    local = initialize_local_members(adapter, topology, icesee_kwargs)

    serial = initial.copy()
    observation_rows = np.arange(3, dtype=np.int64)
    observation_values = np.array([1.4, 2.2, 3.1], dtype=float)
    obs_layout = contiguous_observation_layout(observation_rows, topology)
    validate_observation_partition(obs_layout, topology.spatial_comm)
    maximum_error = 0.0
    first_failure: tuple[int, int, float] | None = None
    for timestep in range(args.steps + 1):
        local_records = {}
        for member_id, local_state in local.members.items():
            reconstructed = reconstruct_global_member_for_parity(
                local_state,
                local.layout,
                topology.spatial_comm,
            )
            if topology.is_spatial_root:
                local_records[member_id] = reconstructed

        records = world.gather(
            local_records if topology.is_spatial_root else None,
            root=0,
        )
        if world.Get_rank() == 0:
            distributed: dict[int, np.ndarray] = {}
            for record in records:
                if record:
                    distributed.update(record)
            if tuple(sorted(distributed)) != tuple(range(args.nens)):
                raise RuntimeError("distributed schedule did not reconstruct every member")
            matrix = np.column_stack([distributed[i] for i in range(args.nens)])
            difference = np.abs(matrix - serial)
            step_error = float(np.max(difference, initial=0.0))
            maximum_error = max(maximum_error, step_error)
            if first_failure is None and not np.allclose(
                matrix,
                serial,
                rtol=args.rtol,
                atol=args.atol,
            ):
                row, member = np.unravel_index(np.argmax(difference), difference.shape)
                first_failure = (timestep, member, float(difference[row, member]))

        if timestep < args.steps:
            local = forecast_local_members(
                adapter,
                local,
                timestep,
                topology,
                icesee_kwargs,
            )
            if world.Get_rank() == 0:
                serial = np.column_stack(
                    [lorenz_rk4_step(serial[:, i], icesee_kwargs) for i in range(args.nens)]
                )

            if timestep + 1 == args.analysis_step:
                local_observations = {
                    member_id: adapter.observe_local_member(
                        local_state,
                        obs_layout.global_row_ids,
                        layout=local.layout,
                        topology=topology,
                        icesee_kwargs=icesee_kwargs,
                    )
                    for member_id, local_state in local.members.items()
                }
                products = StochasticAnalysisProducts.zeros(args.nens)
                observation_blocks = iter_local_ensemble_row_blocks(
                    local_observations,
                    args.nens,
                    topology.ensemble_comm,
                    row_chunk_size=2,
                )
                local_values = observation_values[
                    obs_layout.owned_start:obs_layout.owned_stop
                ]
                value_offset = 0
                for _, forecast_observations in observation_blocks:
                    width = forecast_observations.shape[0]
                    products.add_chunk(
                        forecast_observations,
                        local_values[value_offset:value_offset + width],
                        error_mode="legacy_prior_anomalies",
                    )
                    value_offset += width
                products = products.reduced(topology.spatial_comm)
                transform = ensemble_transform_from_products(products)
                local = apply_ensemble_transform_local(
                    adapter,
                    local,
                    transform,
                    timestep + 1,
                    topology,
                    icesee_kwargs,
                    row_chunk_size=2,
                )

                if world.Get_rank() == 0:
                    serial_products = StochasticAnalysisProducts.zeros(args.nens)
                    serial_products.add_chunk(
                        serial,
                        observation_values,
                        error_mode="legacy_prior_anomalies",
                    )
                    serial_transform = ensemble_transform_from_products(
                        serial_products
                    )
                    serial = serial @ serial_transform

    passed = None
    if world.Get_rank() == 0:
        passed = first_failure is None
        print("Mode-3 Lorenz forecast and stochastic-analysis parity")
        print(f"  MPI process grid: {topology.ensemble_groups} x {topology.spatial_ranks}")
        print(f"  ensemble members: {args.nens}")
        print(f"  compared states:  {args.steps + 1}")
        print(f"  analysis step:    {args.analysis_step}")
        print(f"  max absolute error: {maximum_error:.12g}")
        if first_failure is not None:
            step, member, error = first_failure
            print(f"  first mismatch: step={step}, member={member}, max_abs={error:.12g}")
        print(f"  result: {'PASS' if passed else 'FAIL'}")
    passed = world.bcast(passed, root=0)
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
