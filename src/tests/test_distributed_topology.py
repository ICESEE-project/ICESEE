import pytest

from src.parallelization.distributed_topology import create_distributed_topology


class FakeComm:
    def __init__(self, size, rank):
        self.size = size
        self.rank = rank
        self.splits = []

    def Get_size(self):
        return self.size

    def Get_rank(self):
        return self.rank

    def Split(self, color, key):
        result = (color, key)
        self.splits.append(result)
        return result


def test_two_dimensional_coordinates_and_communicators():
    world = FakeComm(size=12, rank=7)

    topology = create_distributed_topology(
        world, ensemble_groups=3, spatial_ranks=4
    )

    assert topology.ensemble_slot == 1
    assert topology.spatial_rank == 3
    assert topology.spatial_comm == (1, 3)
    assert topology.ensemble_comm == (3, 1)
    assert world.splits == [(1, 3), (3, 1)]
    assert not topology.is_spatial_root
    assert not topology.is_ensemble_root


@pytest.mark.parametrize(
    ("kwargs", "expected"),
    [
        ({"ensemble_groups": 3}, (3, 4)),
        ({"spatial_ranks": 4}, (3, 4)),
    ],
)
def test_infers_one_process_grid_dimension(kwargs, expected):
    topology = create_distributed_topology(FakeComm(12, 0), **kwargs)
    assert (topology.ensemble_groups, topology.spatial_ranks) == expected
    assert topology.is_spatial_root
    assert topology.is_ensemble_root


@pytest.mark.parametrize(
    "kwargs",
    [
        {},
        {"ensemble_groups": 0},
        {"spatial_ranks": 0},
        {"ensemble_groups": 5},
        {"spatial_ranks": 5},
        {"ensemble_groups": 2, "spatial_ranks": 4},
    ],
)
def test_rejects_invalid_process_grids(kwargs):
    with pytest.raises(ValueError):
        create_distributed_topology(FakeComm(12, 0), **kwargs)
