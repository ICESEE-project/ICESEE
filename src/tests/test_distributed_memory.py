import numpy as np
import pytest

from src.parallelization.distributed_memory import estimate_mode3_memory


def test_memory_plan_scales_with_owned_slab_not_global_member():
    one_spatial_rank = estimate_mode3_memory(
        global_rows=5_000_000_000,
        total_members=40,
        ensemble_groups=8,
        spatial_ranks=1,
        state_row_chunk_size=2048,
        observation_row_chunk_size=1024,
        dtype=np.float64,
    )
    sixteen_spatial_ranks = estimate_mode3_memory(
        global_rows=5_000_000_000,
        total_members=40,
        ensemble_groups=8,
        spatial_ranks=16,
        state_row_chunk_size=2048,
        observation_row_chunk_size=1024,
        dtype=np.float64,
    )

    assert sixteen_spatial_ranks.owned_rows == 312_500_000
    assert sixteen_spatial_ranks.local_members == 5
    assert (
        sixteen_spatial_ranks.state_transform_workspace_bytes
        == one_spatial_rank.state_transform_workspace_bytes
    )
    assert (
        sixteen_spatial_ranks.observation_workspace_bytes
        == one_spatial_rank.observation_workspace_bytes
    )
    assert (
        sixteen_spatial_ranks.owned_snapshots_bytes
        * 16
        == one_spatial_rank.owned_snapshots_bytes
    )
    assert (
        sixteen_spatial_ranks.estimated_peak_icesee_bytes
        < sixteen_spatial_ranks.complete_local_members_bytes
    )
    assert (
        sixteen_spatial_ranks.estimated_minimum_peak_bytes
        == sixteen_spatial_ranks.estimated_peak_icesee_bytes
        + sixteen_spatial_ranks.minimum_native_state_bytes
    )


def test_memory_plan_accounts_for_dtype_and_round_robin_members():
    plan = estimate_mode3_memory(
        global_rows=101,
        total_members=7,
        ensemble_groups=3,
        spatial_ranks=4,
        dtype=np.float32,
    )
    assert plan.owned_rows == 26
    assert plan.local_members == 3
    assert plan.itemsize == 4
    assert plan.owned_snapshots_bytes == 2 * 26 * 3 * 4


@pytest.mark.parametrize(
    "argument,value",
    [
        ("global_rows", 0),
        ("total_members", -1),
        ("ensemble_groups", 0),
        ("spatial_ranks", 0),
        ("state_row_chunk_size", 0),
        ("observation_row_chunk_size", 0),
    ],
)
def test_memory_plan_rejects_nonpositive_sizes(argument, value):
    kwargs = dict(
        global_rows=100,
        total_members=4,
        ensemble_groups=2,
        spatial_ranks=2,
    )
    kwargs[argument] = value
    with pytest.raises(ValueError):
        estimate_mode3_memory(**kwargs)


def test_memory_plan_rejects_more_ensemble_groups_than_members():
    with pytest.raises(ValueError, match="cannot exceed"):
        estimate_mode3_memory(
            global_rows=100,
            total_members=4,
            ensemble_groups=5,
            spatial_ranks=2,
        )
