#!/usr/bin/env python3
"""Interrupted/restarted mode-3 trajectory parity across MPI process grids."""

from __future__ import annotations

import argparse
from pathlib import Path
import shutil
import sys
import tempfile

from mpi4py import MPI
import numpy as np


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from src.parallelization.distributed_adapter import contiguous_state_layout
from src.parallelization.distributed_checkpoint import (
    load_distributed_checkpoint,
    save_distributed_checkpoint,
)
from src.parallelization.distributed_runtime import (
    LocalMemberEnsemble,
    apply_ensemble_transform_local,
    members_for_ensemble_slot,
)
from src.parallelization.distributed_topology import create_distributed_topology


class _IdentityFinalizer:
    """Minimal adapter surface required by the transform runtime."""

    @staticmethod
    def finalize_local_analysis(
        local_forecast,
        local_analysis,
        member_id,
        timestep,
        *,
        layout,
        topology,
        icesee_kwargs,
    ):
        del local_forecast, member_id, timestep, layout, topology, icesee_kwargs
        return np.asarray(local_analysis)


def _arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nens", type=int, default=7)
    parser.add_argument("--global-size", type=int, default=103)
    parser.add_argument("--cycles", type=int, default=8)
    parser.add_argument("--checkpoint-cycle", type=int, default=3)
    parser.add_argument("--source-ensemble-groups", type=int, default=2)
    parser.add_argument("--target-ensemble-groups", type=int, default=1)
    parser.add_argument("--row-chunk-size", type=int, default=11)
    parser.add_argument("--atol", type=float, default=2.0e-12)
    parser.add_argument("--rtol", type=float, default=2.0e-12)
    return parser.parse_args()


def _initial(member_id: int, start: int, stop: int) -> np.ndarray:
    rows = np.arange(start, stop, dtype=np.float64)
    return 0.2 * member_id + np.sin(0.017 * rows) + 1.0e-4 * rows**2


def _forecast(values, member_id: int, cycle: int, start: int) -> np.ndarray:
    rows = np.arange(start, start + values.size, dtype=np.float64)
    return (
        (1.0 + 2.0e-4 * (cycle + 1)) * values
        + 5.0e-4 * (member_id + 1)
        + 2.0e-5 * (cycle + 1) * rows
    )


def _transform(number_of_members: int, cycle: int) -> np.ndarray:
    """Stable dense right transform exercising ensemble communication."""

    transform = np.eye(number_of_members, dtype=np.float64)
    strength = 0.015 / (cycle + 2)
    for member in range(number_of_members):
        other = (member + 1) % number_of_members
        transform[member, member] -= strength
        transform[other, member] += strength
    return transform


def _advance_distributed(local, topology, first, stop, args):
    adapter = _IdentityFinalizer()
    kwargs = {"Nens": args.nens}
    current = local
    for cycle in range(first, stop):
        forecast = LocalMemberEnsemble(
            current.layout,
            {
                member: _forecast(
                    values, member, cycle, current.layout.owned_start
                )
                for member, values in current.members.items()
            },
        )
        current = apply_ensemble_transform_local(
            adapter,
            forecast,
            _transform(args.nens, cycle),
            cycle,
            topology,
            kwargs,
            row_chunk_size=args.row_chunk_size,
        )
    return current


def _advance_reference(state: np.ndarray, first: int, stop: int) -> np.ndarray:
    current = state.copy()
    for cycle in range(first, stop):
        for member in range(current.shape[1]):
            current[:, member] = _forecast(
                current[:, member], member, cycle, 0
            )
        current = current @ _transform(current.shape[1], cycle)
    return current


def _collect_global(local, topology, number_of_members: int):
    payload = [
        (
            int(member),
            int(local.layout.owned_start),
            int(local.layout.owned_stop),
            np.asarray(values),
        )
        for member, values in local.members.items()
    ]
    gathered = topology.world.gather(payload, root=0)
    if topology.world_rank != 0:
        return None
    result = np.empty((local.layout.global_size, number_of_members))
    coverage = np.zeros_like(result, dtype=np.int8)
    for rank_payload in gathered:
        for member, start, stop, values in rank_payload:
            result[start:stop, member] = values
            coverage[start:stop, member] += 1
    if not np.all(coverage == 1):
        raise ValueError("final distributed state is not owned exactly once")
    return result


def main() -> int:
    args = _arguments()
    world = MPI.COMM_WORLD
    if not 0 < args.checkpoint_cycle < args.cycles:
        raise ValueError("checkpoint-cycle must lie strictly inside the trajectory")
    for groups in (args.source_ensemble_groups, args.target_ensemble_groups):
        if world.Get_size() % groups:
            raise ValueError("world size is incompatible with a requested process grid")

    source = create_distributed_topology(
        world, ensemble_groups=args.source_ensemble_groups
    )
    source_layout = contiguous_state_layout(
        args.global_size, source, layout_id="restart-trajectory-v1"
    )
    source_members = members_for_ensemble_slot(
        args.nens, source.ensemble_groups, source.ensemble_slot
    )
    local = LocalMemberEnsemble(
        source_layout,
        {
            member: _initial(
                member, source_layout.owned_start, source_layout.owned_stop
            )
            for member in source_members
        },
    )
    local = _advance_distributed(
        local, source, 0, args.checkpoint_cycle, args
    )

    checkpoint_root = world.bcast(
        tempfile.mkdtemp(prefix="icesee-mode3-restart-")
        if world.Get_rank() == 0
        else None,
        root=0,
    )
    saved = save_distributed_checkpoint(
        checkpoint_root,
        args.checkpoint_cycle,
        local,
        source,
        run_id="mode3-restart-parity",
        metadata={
            "Nens": args.nens,
            "next_cycle": args.checkpoint_cycle,
            "random_stream_counter": args.checkpoint_cycle * args.nens,
        },
    )

    target = create_distributed_topology(
        world, ensemble_groups=args.target_ensemble_groups
    )
    target_layout = contiguous_state_layout(
        args.global_size, target, layout_id="restart-trajectory-v1"
    )
    restarted, checkpoint = load_distributed_checkpoint(
        saved.path,
        target_layout,
        target,
        expected_run_id="mode3-restart-parity",
        number_of_members=args.nens,
    )
    restarted = _advance_distributed(
        restarted,
        target,
        int(checkpoint.metadata["next_cycle"]),
        args.cycles,
        args,
    )
    distributed = _collect_global(restarted, target, args.nens)

    passed = None
    if world.Get_rank() == 0:
        initial = np.column_stack(
            [_initial(member, 0, args.global_size) for member in range(args.nens)]
        )
        reference = _advance_reference(initial, 0, args.cycles)
        maximum_error = float(
            np.max(np.abs(distributed - reference), initial=0.0)
        )
        passed = bool(
            np.allclose(distributed, reference, atol=args.atol, rtol=args.rtol)
        )
        print("Mode-3 interrupted/restarted trajectory parity")
        print(
            f"  source grid:        {args.source_ensemble_groups} x "
            f"{world.Get_size() // args.source_ensemble_groups}"
        )
        print(
            f"  restart grid:       {args.target_ensemble_groups} x "
            f"{world.Get_size() // args.target_ensemble_groups}"
        )
        print(f"  checkpoint cycle:   {args.checkpoint_cycle}")
        print(f"  completed cycles:   {args.cycles}")
        print(f"  max absolute error: {maximum_error:.12g}")
        print(f"  result:             {'PASS' if passed else 'FAIL'}")
    passed = world.bcast(passed, root=0)
    world.Barrier()
    if world.Get_rank() == 0:
        shutil.rmtree(checkpoint_root)
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
