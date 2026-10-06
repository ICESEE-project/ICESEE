# ==============================================================================
# @des: Operation-count, traffic, and preflight-projection formulas for
# execution-mode-3's Option-C (out-of-core member store) path. Pure
# arithmetic, no MPI/Firedrake/HDF5 -- validated against real measured
# Icepack instrumentation (see module docstring below and this session's
# derivation), then used to project Ne=40..1000 without fabricating
# results at scales this Mac cannot run.
# @date: 2026-09-27
# ==============================================================================
"""Derivation and validation (Phase 1 of the four-hour PACE-prep session).

The store-streaming analysis path (``distributed_streaming_runtime.py``)
issues, per ``run_native_store_streaming_analysis_cycle`` call:

  - ``blocks_per_cycle = ceil(owned_size_total / row_chunk_size)`` bulk
    get/put pairs -- ONLY on timesteps with a scheduled analysis event
    (``number_of_batches > 0``), since this session's Phase-1 finding
    proved the ensemble transform is exactly ``eye(Ne)`` (bit-for-bit, not
    approximately) when no observation batch contributes, making the full
    row-block pass a mathematically provable no-op on every other
    timestep -- and ``run_native_store_streaming_analysis_cycle`` was
    fixed this session to skip it there (see that module's own comment).

  - ``owned_size_total`` is the SUM across every variable block in one
    member's packed state (e.g. an application with 5 packed scalar
    fields would report 5x one field's node count for P=1) -- NOT one
    field's own size. This was
    the exact source of an early hand-derivation mismatch this session
    (54123 vs the correct 270615), corrected and validated below.

Validated against a real P=1, Nens=2 Icepack run (25 timesteps, 2
scheduled analysis events, row_chunk_size=4096 (the hard-coded function
default -- icepack's mode3_runner.py does not yet expose this as a CLI
override), owned_size_total=270615):

  blocks_per_cycle = ceil(270615 / 4096) = 67
  PRE-optimization bulk_get_count = 25 timesteps * 67 blocks = 1675
    (measured instrumentation, BEFORE this session's transform-skip fix:
    exactly 1675 -- confirms the formula and the (undesirable) pattern)
  POST-optimization bulk_get_count = 2 analysis events * 67 blocks = 134
    (12.5x fewer HDF5 operations for this specific run's schedule; the
    reduction factor is exactly nt / n_analysis_events in general)
"""

from __future__ import annotations

import math
from dataclasses import dataclass


def blocks_per_cycle(owned_size_total: int, row_chunk_size: int) -> int:
    """Number of row blocks one analysis pass touches on one rank."""

    if row_chunk_size <= 0:
        raise ValueError("row_chunk_size must be positive")
    return math.ceil(int(owned_size_total) / int(row_chunk_size))


@dataclass(frozen=True)
class OperationCountProjection:
    ne: int
    ensemble_groups: int
    p_model: int
    row_chunk_size: int
    nx_local: int
    num_variable_blocks: int
    nt: int
    n_analysis_events: int
    rounds_per_slot: int
    blocks_per_cycle: int
    bulk_get_count_optimized: int
    bulk_put_count_optimized: int
    bulk_get_count_unoptimized: int
    physical_hdf5_selections_optimized: int
    physical_hdf5_selections_unoptimized: int
    store_files: int


def project_operation_counts(
    *,
    ne: int,
    ensemble_groups: int,
    p_model: int,
    nx_global: int,
    num_variable_blocks: int,
    row_chunk_size: int,
    nt: int,
    n_analysis_events: int,
) -> OperationCountProjection:
    """Phase 1/2: project HDF5 operation counts for one (Ne, layout) case.

    Does NOT multiply blindly by global Ne -- ``rounds_per_slot`` (members
    actually resident on one rank/slot over the run) is what determines
    how many members' data one rank's bulk call spans; a rank never
    touches another ensemble group's members at all.
    """

    if ensemble_groups <= 0 or p_model <= 0:
        raise ValueError("ensemble_groups and p_model must be positive")
    rounds_per_slot = math.ceil(int(ne) / int(ensemble_groups))
    nx_local = math.ceil(int(nx_global) / int(p_model))
    owned_size_total = nx_local * int(num_variable_blocks)
    blocks = blocks_per_cycle(owned_size_total, row_chunk_size)

    bulk_get_optimized = int(n_analysis_events) * blocks
    bulk_get_unoptimized = int(nt) * blocks

    # member-major's get_ensemble_rows loops internally per member on this
    # rank/slot -- physical HDF5 selections multiply by rounds_per_slot,
    # NOT by global Ne (a rank only ever holds its own slot's members).
    physical_optimized = bulk_get_optimized * rounds_per_slot
    physical_unoptimized = bulk_get_unoptimized * rounds_per_slot

    world_size = int(ensemble_groups) * int(p_model)

    return OperationCountProjection(
        ne=int(ne), ensemble_groups=int(ensemble_groups), p_model=int(p_model),
        row_chunk_size=int(row_chunk_size), nx_local=nx_local,
        num_variable_blocks=int(num_variable_blocks), nt=int(nt),
        n_analysis_events=int(n_analysis_events), rounds_per_slot=rounds_per_slot,
        blocks_per_cycle=blocks,
        bulk_get_count_optimized=bulk_get_optimized,
        bulk_put_count_optimized=bulk_get_optimized,
        bulk_get_count_unoptimized=bulk_get_unoptimized,
        physical_hdf5_selections_optimized=physical_optimized,
        physical_hdf5_selections_unoptimized=physical_unoptimized,
        store_files=world_size,
    )


@dataclass(frozen=True)
class TrafficProjection:
    """Phase 7: complete per-cycle logical traffic, forecast vs analysis
    separated (Phase 12/13). ``logical_member_bytes`` is one member's full
    dynamic-state size (e.g. 37e9 for the stress case)."""

    logical_member_bytes: float
    ne: int
    forecast_activation_read_bytes_per_timestep: float
    forecast_deactivation_write_bytes_per_timestep: float
    analysis_read_bytes_per_event: float
    analysis_write_bytes_per_event: float
    total_forecast_bytes_per_run: float
    total_analysis_bytes_per_run: float
    total_bytes_per_run: float


def project_traffic(
    *, logical_member_bytes: float, ne: int, nt: int, n_analysis_events: int
) -> TrafficProjection:
    """Phase 7's 4-pass model, corrected for the transform-skip
    optimization: forecast activation/deactivation happens every timestep
    (unavoidable -- every member forecasts every timestep in a standard
    EnKF cycle); analysis read/write happens ONLY on scheduled analysis
    events (this session's fix), not on every timestep as the original
    un-optimized 74TB-style estimate implicitly assumed.

    L = logical_member_bytes * ne (one variable's full ensemble size).
    Per timestep: activation read L + deactivation write L = 2L.
    Per analysis event: analysis read L + analysis write L = 2L (this is
    where the earlier "74TB per analysis transformation" figure came
    from -- it was always an accurate PER-EVENT number, never a per-run
    total; this function makes that distinction explicit).
    """

    L = float(logical_member_bytes) * int(ne)
    forecast_read = L
    forecast_write = L
    analysis_read = L
    analysis_write = L
    total_forecast = int(nt) * (forecast_read + forecast_write)
    total_analysis = int(n_analysis_events) * (analysis_read + analysis_write)
    return TrafficProjection(
        logical_member_bytes=float(logical_member_bytes), ne=int(ne),
        forecast_activation_read_bytes_per_timestep=forecast_read,
        forecast_deactivation_write_bytes_per_timestep=forecast_write,
        analysis_read_bytes_per_event=analysis_read,
        analysis_write_bytes_per_event=analysis_write,
        total_forecast_bytes_per_run=total_forecast,
        total_analysis_bytes_per_run=total_analysis,
        total_bytes_per_run=total_forecast + total_analysis,
    )


@dataclass(frozen=True)
class StoragePreflight:
    ne: int
    ensemble_groups: int
    p_model: int
    world_size: int
    rounds_per_slot: int
    logical_ensemble_bytes: float
    estimated_members_per_slot: int
    estimated_inactive_store_bytes_per_rank: float
    estimated_total_temporary_store_bytes: float
    store_files: int


def project_storage_preflight(
    *, ne: int, member_state_bytes: float, ensemble_groups: int, p_model: int
) -> StoragePreflight:
    """Phase 10: backend-neutral store-capacity preflight. Deliberately a
    standalone function, NOT a change to the production ``ResourcePlan``
    class (``src/parallelization/parallel_mpi/resource_plan.py``) -- that
    class's ``ranks_per_model``/``num_model_groups`` belong to modes 0-2's
    own, separately-established topology system, and conflating it with
    Mode-3's ``model_nprocs``/``ensemble_groups`` was already flagged as a
    recurring confusion earlier in this engagement. This function's
    outputs are what a future ``Mode3StoragePlan`` would compute, kept
    here as plain, testable arithmetic in the meantime.
    """

    if ensemble_groups <= 0 or p_model <= 0:
        raise ValueError("ensemble_groups and p_model must be positive")
    rounds_per_slot = math.ceil(int(ne) / int(ensemble_groups))
    world_size = int(ensemble_groups) * int(p_model)
    per_rank_member_bytes = float(member_state_bytes) / int(p_model)
    inactive_bytes_per_rank = rounds_per_slot * per_rank_member_bytes
    return StoragePreflight(
        ne=int(ne), ensemble_groups=int(ensemble_groups), p_model=int(p_model),
        world_size=world_size, rounds_per_slot=rounds_per_slot,
        logical_ensemble_bytes=float(member_state_bytes) * int(ne),
        estimated_members_per_slot=rounds_per_slot,
        estimated_inactive_store_bytes_per_rank=inactive_bytes_per_rank,
        estimated_total_temporary_store_bytes=inactive_bytes_per_rank * world_size,
        store_files=world_size,
    )


@dataclass(frozen=True)
class MemoryPreflight:
    active_native_state_bytes_per_rank: float
    inactive_ram_bytes_per_rank: float  # 0 for a file-backed store
    analysis_workspace_bytes_per_rank: float
    ensemble_space_workspace_bytes_per_rank: float
    approximate_total_working_bytes_per_rank: float


def project_memory_preflight(
    *,
    member_state_bytes: float,
    p_model: int,
    ne: int,
    ensemble_groups: int,
    row_chunk_size: int,
    backend: str,
    workspace_constant: float = 7.1,
) -> MemoryPreflight:
    """Phase 11: distinguish LOGICAL bytes from RESIDENT bytes per rank.

    ``workspace_constant`` (Phase 6) accounts for simultaneously-live
    buffers in the row-block transform (the input block, the gathered
    cross-rank columns, the output analysis block, etc.) -- measure it,
    do not assume 1; see ``recommend_row_chunk_size`` below for how this
    session's benchmark measurement informs a starting value.
    """

    per_rank_member_bytes = float(member_state_bytes) / int(p_model)
    active_native_bytes = per_rank_member_bytes  # bounded to ONE member (Option B)

    rounds_per_slot = math.ceil(int(ne) / int(ensemble_groups))
    if str(backend).lower() == "memory":
        inactive_ram = rounds_per_slot * per_rank_member_bytes
    else:
        inactive_ram = 0.0  # file-backed: resident RAM does not scale with rounds

    analysis_workspace = float(workspace_constant) * int(row_chunk_size) * int(ne) * 8
    ensemble_space_workspace = 6 * int(ne) * int(ne) * 8  # cross/gram/rhs + eigh temporaries

    total = active_native_bytes + inactive_ram + analysis_workspace + ensemble_space_workspace
    return MemoryPreflight(
        active_native_state_bytes_per_rank=active_native_bytes,
        inactive_ram_bytes_per_rank=inactive_ram,
        analysis_workspace_bytes_per_rank=analysis_workspace,
        ensemble_space_workspace_bytes_per_rank=ensemble_space_workspace,
        approximate_total_working_bytes_per_rank=total,
    )


def recommend_row_chunk_size(
    *, ne: int, memory_budget_bytes: float, safety_fraction: float = 0.5,
    workspace_constant: float = 7.1,
) -> int:
    """Phase 6/12: advisory maximum safe row_chunk_size (B) for a given
    per-rank analysis-workspace memory budget.

    B <= (safety_fraction * memory_budget) / (workspace_constant * Ne * 8)

    Advisory only -- not wired into any production default.
    """

    if ne <= 0 or memory_budget_bytes <= 0:
        raise ValueError("ne and memory_budget_bytes must be positive")
    usable = float(safety_fraction) * float(memory_budget_bytes)
    denom = float(workspace_constant) * int(ne) * 8
    return max(1, int(usable // denom))
