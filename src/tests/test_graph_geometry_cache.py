import numpy as np
from ICESEE.src.run_model_da import _error_generation as errors


def reset_cache():
    errors._GRAPH_TRANSITION_CACHE.clear()
    errors._GRAPH_TRANSITION_CACHE_BYTES=0


def draw(coords,seed):
    return errors.generate_pseudo_random_field_1d(
        N=len(coords),Lx=10.,rh=2.,method="graph",coords=coords,
        seed=seed,num_passes=4,k_neighbors=3)


def test_graph_cache_preserves_fields_and_random_state(monkeypatch):
    coords=np.column_stack((np.arange(30.),np.sin(np.arange(30.))))
    reset_cache()
    monkeypatch.setattr(errors,"_GRAPH_TRANSITION_CACHE_LIMIT",0)
    baseline=[draw(coords,seed) for seed in (12,24,41)]
    reset_cache()
    monkeypatch.setattr(errors,"_GRAPH_TRANSITION_CACHE_LIMIT",64*1024**2)
    np.random.seed(271)
    before=np.random.get_state()
    cached=[draw(coords,seed) for seed in (12,24,41)]
    after=np.random.get_state()
    for a,b in zip(baseline,cached):
        np.testing.assert_array_equal(a,b)
    np.testing.assert_array_equal(before[1],after[1])
    assert before[2:]==after[2:]
    assert len(errors._GRAPH_TRANSITION_CACHE)==1


def test_mutating_coordinates_cannot_reuse_a_stale_graph(monkeypatch):
    coords=np.column_stack((np.arange(20.),np.zeros(20)))
    reset_cache()
    first=draw(coords,12)
    coords[8]=[40.,10.]
    changed=draw(coords,12)
    assert len(errors._GRAPH_TRANSITION_CACHE)==2
    reset_cache()
    monkeypatch.setattr(errors,"_GRAPH_TRANSITION_CACHE_LIMIT",0)
    expected=draw(coords,12)
    np.testing.assert_array_equal(changed,expected)
    assert not np.array_equal(first,changed)
