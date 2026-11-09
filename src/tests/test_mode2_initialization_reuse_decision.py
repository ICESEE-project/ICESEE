"""Regression for rank-local startup probes racing with shard creation."""
from types import SimpleNamespace
import pytest
from ICESEE.src.run_model_da import icesee_da_full_parallel as runner


class Rank:
    def __init__(self, rank, shared):
        self.rank, self.shared = rank, shared

    def Get_rank(self):
        return self.rank

    def bcast(self, value, root):
        assert root == 0
        if self.rank == root:
            self.shared["decision"] = value
        else:
            assert value is None
        return self.shared["decision"]


def decide(comm, **overrides):
    args = dict(run_already_complete=False, reuse_allowed=True,
                state_file="initial.h5", nd=12, nens=32)
    args.update(overrides)
    return runner._initialization_reuse_decision(comm, **args)


def test_newly_created_shard_cannot_split_ranks_between_reuse_and_initialization(monkeypatch):
    shared = {}
    present = {"value": False}
    probes = []

    def probe(path, dataset, shape):
        probes.append((path, dataset, shape))
        return present["value"]

    monkeypatch.setattr(runner.tools, "h5_has_dataset_with_shape", probe)
    results = [decide(Rank(0, shared))]
    # Simulate rank zero creating a correctly shaped shard before a late rank
    # reaches the decision. The late rank must consume the broadcast, not
    # inspect the newly created (but not yet initialized) file.
    present["value"] = True
    results.extend(decide(Rank(rank, shared)) for rank in [5, 1, 7, 2, 3, 4, 6])
    assert results == [False] * 8
    assert probes == [("initial.h5", "states", (12, 32))]


def test_valid_saved_state_is_reused_on_all_ranks(monkeypatch):
    monkeypatch.setattr(runner.tools, "h5_has_dataset_with_shape", lambda *args: True)
    shared = {}
    assert all(decide(Rank(rank, shared)) for rank in range(8))


@pytest.mark.parametrize("completed,reuse,expected", [(True, True, True), (False, False, False)])
def test_completed_and_forced_fresh_paths_do_not_probe_files(monkeypatch, completed, reuse, expected):
    def forbidden(*args):
        raise AssertionError("file probe must be skipped")
    monkeypatch.setattr(runner.tools, "h5_has_dataset_with_shape", forbidden)
    shared = {}
    assert [decide(Rank(rank, shared), run_already_complete=completed,
                   reuse_allowed=reuse) for rank in range(8)] == [expected] * 8
