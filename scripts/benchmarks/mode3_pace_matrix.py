# ==============================================================================
# @des: Staged PACE experiment matrix for Mode-3's Option-C investigation.
# Purely algebraic (uses src/parallelization/mode3_projections.py's
# validated formulas) -- prints/writes the layouts and PROJECTED metrics
# to run, never fabricated PACE results. No hardcoded account/queue/
# partition/node-count details; world_size is always expressed
# algebraically as ensemble_groups * p_model so the SAME matrix can be
# submitted against whatever resource allocation is actually available.
# @date: 2026-09-27
# ==============================================================================
from __future__ import annotations

import argparse
import json
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..")))

from ICESEE.src.parallelization.mode3_projections import (
    project_memory_preflight,
    project_operation_counts,
    project_storage_preflight,
    project_traffic,
)

# Stage A: correctness/performance baseline -- start here, do not skip ahead.
STAGE_A = [40, 80, 100]
# Stage B: medium scale -- only after Stage A succeeds.
STAGE_B = [200, 300]
# Stage C: large scale -- only after Stage B succeeds.
STAGE_C = [500, 1000]

# Representative layouts per Ne: (ensemble_groups, p_model). world_size =
# ensemble_groups * p_model. These are STARTING POINTS for PACE submission,
# not guaranteed-optimal -- Phase 17/18/19 (strong-scaling, ensemble-
# concurrency, round-scaling) exist specifically to find better ones.
LAYOUTS_PER_NE = {
    40: [(40, 1), (20, 2), (10, 4)],
    80: [(80, 1), (40, 2), (20, 4)],
    100: [(100, 1), (50, 2), (25, 4)],
    200: [(200, 1), (100, 2), (50, 4), (25, 8)],
    300: [(300, 1), (150, 2), (75, 4)],
    500: [(500, 1), (250, 2), (125, 4)],
    1000: [(1000, 1), (500, 2), (250, 4), (125, 8)],
}

ROW_CHUNK_SIZES = [4096, 16384, 65536]
STORE_BACKENDS = ["memory", "hdf5_member_major"]


def build_matrix(member_state_gb: float, nx_global: int, nt: int, n_analysis_events: int):
    matrix = []
    for stage_name, ne_values in (("A", STAGE_A), ("B", STAGE_B), ("C", STAGE_C)):
        for ne in ne_values:
            for (ensemble_groups, p_model) in LAYOUTS_PER_NE[ne]:
                for row_chunk_size in ROW_CHUNK_SIZES:
                    for backend in STORE_BACKENDS:
                        if backend == "memory" and ne >= 500:
                            # memory remains the correctness/performance
                            # reference at small/medium scale only -- not
                            # expected (and not claimed) to fit Ne>=500 at
                            # the 37GB/member stress size; still included
                            # at smaller Ne for backend comparison (Gate 20).
                            continue
                        ops = project_operation_counts(
                            ne=ne, ensemble_groups=ensemble_groups, p_model=p_model,
                            nx_global=nx_global, num_variable_blocks=5,
                            row_chunk_size=row_chunk_size, nt=nt,
                            n_analysis_events=n_analysis_events,
                        )
                        traffic = project_traffic(
                            logical_member_bytes=member_state_gb * 1e9, ne=ne,
                            nt=nt, n_analysis_events=n_analysis_events,
                        )
                        storage = project_storage_preflight(
                            ne=ne, member_state_bytes=member_state_gb * 1e9,
                            ensemble_groups=ensemble_groups, p_model=p_model,
                        )
                        memory = project_memory_preflight(
                            member_state_bytes=member_state_gb * 1e9, p_model=p_model,
                            ne=ne, ensemble_groups=ensemble_groups,
                            row_chunk_size=row_chunk_size, backend=backend,
                        )
                        matrix.append({
                            "stage": stage_name,
                            "ne": ne,
                            "ensemble_groups": ensemble_groups,
                            "p_model": p_model,
                            "world_size": ensemble_groups * p_model,
                            "row_chunk_size": row_chunk_size,
                            "store_backend": backend,
                            "rounds_per_slot": ops.rounds_per_slot,
                            "store_files": ops.store_files,
                            "bulk_ops_per_run_optimized": ops.bulk_get_count_optimized,
                            "physical_hdf5_selections_optimized": ops.physical_hdf5_selections_optimized,
                            "logical_ensemble_gb": storage.logical_ensemble_bytes / 1e9,
                            "total_traffic_gb_per_run": traffic.total_bytes_per_run / 1e9,
                            "inactive_ram_gb_per_rank": memory.inactive_ram_bytes_per_rank / 1e9,
                            "active_native_gb_per_rank": memory.active_native_state_bytes_per_rank / 1e9,
                        })
    return matrix


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--member-state-gb", type=float, default=37.0)
    parser.add_argument("--nx-global", type=int, default=54123, help="Icepack idealized_pig's own mesh node count, for a concrete reference case")
    parser.add_argument("--nt", type=int, default=25)
    parser.add_argument("--n-analysis-events", type=int, default=2)
    parser.add_argument("--stage-only", type=str, default=None, choices=["A", "B", "C"])
    args = parser.parse_args()

    matrix = build_matrix(args.member_state_gb, args.nx_global, args.nt, args.n_analysis_events)
    if args.stage_only:
        matrix = [row for row in matrix if row["stage"] == args.stage_only]
    print(json.dumps(matrix, indent=2))
    print(f"\n# {len(matrix)} configurations total", file=sys.stderr)


if __name__ == "__main__":
    main()
