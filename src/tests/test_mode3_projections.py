# Check streamed analysis operation counts, memory budgets and storage projections.
from __future__ import annotations

from ICESEE.src.parallelization.mode3_projections import (
    blocks_per_cycle,
    project_memory_preflight,
    project_operation_counts,
    project_storage_preflight,
    project_traffic,
    recommend_row_chunk_size,
)


def test_blocks_per_cycle_matches_measured_icepack_value():
    # 5 variable blocks (h,u,v,s,basal_melt_field) x 54123 nodes = 270615.
    assert blocks_per_cycle(270615, 4096) == 67


def test_operation_count_projection_matches_measured_icepack_run_exactly():
    """For 67 row blocks, compare 25 full passes (1675 reads) with two analysis
    passes (134 reads)."""

    projection = project_operation_counts(
        ne=2, ensemble_groups=1, p_model=1, nx_global=54123, num_variable_blocks=5,
        row_chunk_size=4096, nt=25, n_analysis_events=2,
    )
    assert projection.blocks_per_cycle == 67
    assert projection.bulk_get_count_unoptimized == 1675  # matches measured pre-fix value
    assert projection.bulk_get_count_optimized == 134
    assert projection.rounds_per_slot == 2  # both members share this one slot
    # member-major loops per member internally -> physical selections
    # multiply by rounds_per_slot, not by global Ne.
    assert projection.physical_hdf5_selections_optimized == 134 * 2
    assert projection.store_files == 1  # ensemble_groups=1 x p_model=1


def test_operation_count_projection_never_multiplies_by_global_ne_directly():
    """A rank only ever touches its OWN slot's members -- ensemble_groups=40
    with Ne=40 means rounds_per_slot=1 regardless of how large Ne is."""

    projection = project_operation_counts(
        ne=1000, ensemble_groups=1000, p_model=4, nx_global=1_000_000,
        num_variable_blocks=5, row_chunk_size=4096, nt=100, n_analysis_events=10,
    )
    assert projection.rounds_per_slot == 1
    assert projection.store_files == 4000


def test_traffic_projection_separates_forecast_from_analysis():
    proj = project_traffic(logical_member_bytes=37e9, ne=1000, nt=25, n_analysis_events=2)
    L = 37e9 * 1000
    assert proj.total_forecast_bytes_per_run == 25 * 2 * L
    assert proj.total_analysis_bytes_per_run == 2 * 2 * L
    # The original "74TB" figure was always a PER-EVENT number (2 * L),
    # not a per-run total -- confirm that identity explicitly.
    assert proj.analysis_read_bytes_per_event + proj.analysis_write_bytes_per_event == 2 * L


def test_storage_preflight_basic_shape():
    plan = project_storage_preflight(
        ne=40, member_state_bytes=37e9, ensemble_groups=10, p_model=4
    )
    assert plan.rounds_per_slot == 4
    assert plan.world_size == 40
    assert plan.store_files == 40
    assert plan.estimated_inactive_store_bytes_per_rank == 4 * (37e9 / 4)


def test_memory_preflight_hdf5_backend_does_not_scale_inactive_ram_with_rounds():
    mem_backend = project_memory_preflight(
        member_state_bytes=37e9, p_model=4, ne=40, ensemble_groups=10,
        row_chunk_size=4096, backend="memory",
    )
    hdf5_backend = project_memory_preflight(
        member_state_bytes=37e9, p_model=4, ne=40, ensemble_groups=10,
        row_chunk_size=4096, backend="hdf5_member_major",
    )
    assert mem_backend.inactive_ram_bytes_per_rank > 0
    assert hdf5_backend.inactive_ram_bytes_per_rank == 0.0
    # Active native state and analysis workspace are identical regardless
    # of backend -- only inactive residency differs.
    assert mem_backend.active_native_state_bytes_per_rank == hdf5_backend.active_native_state_bytes_per_rank
    assert mem_backend.analysis_workspace_bytes_per_rank == hdf5_backend.analysis_workspace_bytes_per_rank


def test_recommend_row_chunk_size_respects_budget():
    b = recommend_row_chunk_size(ne=1000, memory_budget_bytes=1e9, safety_fraction=0.5)
    # Verify the recommended chunk fits the budget with the default workspace factor.
    workspace = 7.1 * b * 1000 * 8
    assert workspace <= 0.5 * 1e9
