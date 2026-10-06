# ==============================================================================
# @des: Level-1 unit tests for src/parallelization/parallel_mpi/resource_plan.py.
#       Pure Python -- no MPI launch required.
# ==============================================================================
import pytest

from ICESEE.src.parallelization.parallel_mpi.resource_plan import plan_resources


# (world_size, nens, ranks_per_model, expected_groups, expected_rounds, expected_spare)
CASES = [
    (1, 1, 1, 1, 1, 0),
    (4, 4, 1, 4, 1, 0),
    (8, 4, 1, 4, 1, 4),
    (10, 4, 1, 4, 1, 6),
    (8, 4, 2, 4, 1, 0),
    (16, 4, 4, 4, 1, 0),
    (10, 4, 2, 4, 1, 2),
    (8, 20, 2, 4, 5, 0),
    (3, 4, 1, 3, 2, 0),
    (2, 4, 2, 1, 4, 0),
    (4, 8, 1, 4, 2, 0),
    (4, 5, 2, 2, 3, 0),
    (10, 4, 3, 3, 2, 1),  # world_size not divisible by ranks_per_model
]


@pytest.mark.parametrize("world_size,nens,rpm,exp_groups,exp_rounds,exp_spare", CASES)
def test_planned_topology_matches_spec(world_size, nens, rpm, exp_groups, exp_rounds, exp_spare):
    plan = plan_resources(world_size, nens, rpm)
    assert plan.num_model_groups == exp_groups
    assert plan.num_rounds == exp_rounds
    assert plan.spare_ranks == exp_spare
    assert plan.ranks_per_model == rpm
    assert plan.active_ranks == exp_groups * rpm
    assert plan.active_ranks + plan.spare_ranks == world_size


@pytest.mark.parametrize("world_size,nens,rpm,exp_groups,exp_rounds,exp_spare", CASES)
def test_every_member_covered_exactly_once(world_size, nens, rpm, exp_groups, exp_rounds, exp_spare):
    plan = plan_resources(world_size, nens, rpm)
    covered = plan.members_covered()
    assert sorted(covered) == list(range(nens)), (
        f"expected members 0..{nens - 1} each exactly once, got {sorted(covered)}"
    )
    assert len(covered) == len(set(covered)), "duplicate member assignment in schedule"


def test_p10_nens4_rpm2_matches_worked_example():
    """Stage 4 spec point 7's worked example: 10 ranks, 4 groups of 2
    (8 ranks), 2 spare -- not 5 groups of uneven size."""
    plan = plan_resources(10, 4, 2)
    assert plan.num_model_groups == 4
    assert plan.spare_ranks == 2
    assert plan.num_rounds == 1


def test_p8_nens20_rpm2_matches_worked_example():
    """Stage 4 spec point 6's worked example: 4 simultaneous groups, 5 rounds."""
    plan = plan_resources(8, 20, 2)
    assert plan.num_model_groups == 4
    assert plan.num_rounds == 5
    assert plan.spare_ranks == 0


def test_block_contiguous_group_assignment():
    """P=16, Nens=4, rpm=4: group 0 -> ranks 0-3, group 1 -> ranks 4-7, ...
    -- not the old strided 0,4,8,12 / 1,5,9,13 pattern."""
    plan = plan_resources(16, 4, 4)
    for rank in range(16):
        assert plan.group_id(rank) == rank // 4
        assert plan.rank_in_group(rank) == rank % 4
    # explicit worked check
    assert [plan.group_id(r) for r in range(16)] == [
        0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3,
    ]


def test_spare_ranks_have_no_group():
    plan = plan_resources(10, 4, 2)
    spare = [r for r in range(10) if plan.is_spare(r)]
    assert spare == [8, 9]
    for r in spare:
        assert plan.group_id(r) is None
        assert plan.rank_in_group(r) is None


def test_legacy_default_p_le_nens_gives_ranks_per_model_1():
    """No explicit ranks_per_model, Nens >= P: must reproduce
    ranks_per_model == 1 with every member running concurrently -- the
    invariant ISSM's integration currently depends on."""
    plan = plan_resources(4, 4, None)
    assert plan.ranks_per_model == 1
    assert plan.num_model_groups == 4
    assert plan.num_rounds == 1
    assert plan.spare_ranks == 0


def test_legacy_default_p_gt_nens_uses_auto_policy():
    plan_legacy = plan_resources(10, 4, None)
    plan_auto = plan_resources(10, 4, "auto")
    assert plan_legacy.ranks_per_model == plan_auto.ranks_per_model == 2
    assert plan_legacy.num_model_groups == plan_auto.num_model_groups == 4
    assert plan_legacy.spare_ranks == plan_auto.spare_ranks == 2


def test_rounds_schedule_matches_group_plus_round_times_groups():
    """The round-major assignment (member = group + round*num_groups)
    that the forecast loop's ens_id formula depends on."""
    plan = plan_resources(8, 20, 2)
    for r in range(plan.num_rounds):
        for g in range(plan.num_model_groups):
            expected = g + r * plan.num_model_groups
            actual = plan.member_for(r, g)
            if expected < plan.nens:
                assert actual == expected
            else:
                assert actual is None


def test_insufficient_world_size_raises_clear_error():
    with pytest.raises(ValueError, match="not enough ranks"):
        plan_resources(2, 4, 4)


def test_invalid_ranks_per_model_raises():
    with pytest.raises(ValueError):
        plan_resources(8, 4, 0)
    with pytest.raises(ValueError):
        plan_resources(8, 4, -1)


def test_invalid_world_size_or_nens_raises():
    with pytest.raises(ValueError):
        plan_resources(0, 4, 1)
    with pytest.raises(ValueError):
        plan_resources(8, 0, 1)


def test_out_of_range_rank_queries_raise():
    plan = plan_resources(4, 4, 1)
    with pytest.raises(ValueError):
        plan.group_id(4)
    with pytest.raises(ValueError):
        plan.group_id(-1)


def test_summary_is_one_line_and_reports_key_fields():
    plan = plan_resources(16, 4, 4)
    line = plan.summary()
    assert "\n" not in line
    for token in ("world=16", "ensembles=4", "ranks/model=4", "groups=4", "rounds=1", "spare=0"):
        assert token in line


# --- max_model_groups (Stage 4D.2) -------------------------------------------

def test_max_model_groups_none_is_a_no_op():
    plain = plan_resources(8, 4, 1)
    capped_none = plan_resources(8, 4, 1, max_model_groups=None)
    assert plain == capped_none


def test_max_model_groups_caps_below_the_uncapped_group_count():
    # Real Stage 4D.2 case: world_size=Nens=4, ranks_per_model=1 (ISSM
    # never uses ranks_per_model>1), issm_server_count=2 -- only 2
    # persistent servers active, the other 2 world ranks become ordinary
    # spares, and each server round-robins 2 of the 4 members.
    plan = plan_resources(4, 4, 1, max_model_groups=2)
    assert plan.num_model_groups == 2
    assert plan.num_rounds == 2
    assert plan.spare_ranks == 2
    assert sorted(plan.members_covered()) == [0, 1, 2, 3]


def test_max_model_groups_above_the_uncapped_group_count_is_a_no_op():
    plain = plan_resources(4, 4, 1)
    capped = plan_resources(4, 4, 1, max_model_groups=100)
    assert plain.num_model_groups == capped.num_model_groups == 4


def test_max_model_groups_of_one_forces_every_member_through_one_server():
    plan = plan_resources(4, 4, 1, max_model_groups=1)
    assert plan.num_model_groups == 1
    assert plan.num_rounds == 4
    assert plan.spare_ranks == 3
    assert sorted(plan.members_covered()) == [0, 1, 2, 3]


def test_max_model_groups_below_one_raises():
    with pytest.raises(ValueError):
        plan_resources(4, 4, 1, max_model_groups=0)
    with pytest.raises(ValueError):
        plan_resources(4, 4, 1, max_model_groups=-1)


def test_describe_all_ranks_covers_every_world_rank_including_spares():
    plan = plan_resources(10, 4, 2)
    dump = plan.describe_all_ranks()
    lines = dump.splitlines()
    # header/summary + one row per world rank
    assert len(lines) == 2 + plan.world_size
    for rank in (8, 9):  # spare ranks for this topology
        assert any(f"{rank:4d} | spare" in line for line in lines)
    for rank in range(8):
        assert any(f"{rank:4d} | " in line and "spare" not in line for line in lines)
