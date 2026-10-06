# ==============================================================================
# @des: Level-1 unit tests for src/parallelization/parallel_mpi/model_capabilities.py.
# Pure Python -- no MPI launch required.
# ==============================================================================
import pytest

from ICESEE.src.parallelization.parallel_mpi.model_capabilities import (
    ModelCapabilities,
    get_model_capabilities,
    register_model_capabilities,
    validate_ranks_per_model_request,
    validate_stochastic_method_for_decomposition,
)


@pytest.fixture(autouse=True)
def _clean_registry():
    from ICESEE.src.parallelization.parallel_mpi import model_capabilities as module

    before = set(module._REGISTRY)
    yield
    for name in set(module._REGISTRY) - before:
        del module._REGISTRY[name]


def test_all_four_current_applications_are_registered():
    for name in ("lorenz", "issm", "icepack", "flowline"):
        caps = get_model_capabilities(name)
        assert isinstance(caps, ModelCapabilities)
        assert caps.model_name == name


def test_lorenz_issm_flowline_conservatively_do_not_support_multi_rank():
    # Honest baseline: none of these three has been run/verified with
    # ranks_per_model > 1 in modes 0-2 in this codebase. icepack is the
    # one exception (Stage 4C, fifth continuation): its multi-rank path
    # is fixed and validated, see test_icepack_multirank_request_accepted
    # and model_capabilities.py's icepack notes.
    for name in ("lorenz", "issm", "flowline"):
        assert get_model_capabilities(name).supports_multi_rank_per_model is False
    assert get_model_capabilities("icepack").supports_multi_rank_per_model is True


def test_unregistered_model_gets_conservative_default_not_a_crash():
    caps = get_model_capabilities("some-future-model-nobody-registered")
    assert caps.supports_multi_rank_per_model is False
    assert caps.requires_one_rank_per_member is False


def test_none_model_name_gets_conservative_default():
    caps = get_model_capabilities(None)
    assert caps.supports_multi_rank_per_model is False


def test_register_and_get_round_trip():
    register_model_capabilities(
        "unit-test-model",
        supports_multi_rank_per_model=True,
        supports_distributed_state=True,
        default_ranks_per_model=2,
        notes="test only",
    )
    caps = get_model_capabilities("unit-test-model")
    assert caps.supports_multi_rank_per_model is True
    assert caps.supports_distributed_state is True
    assert caps.default_ranks_per_model == 2
    assert caps.notes == "test only"


def test_reregistering_same_model_name_overwrites():
    register_model_capabilities("unit-test-model", supports_multi_rank_per_model=False)
    register_model_capabilities("unit-test-model", supports_multi_rank_per_model=True)
    assert get_model_capabilities("unit-test-model").supports_multi_rank_per_model is True


def test_register_rejects_empty_name():
    with pytest.raises(ValueError):
        register_model_capabilities("")


# --- validate_ranks_per_model_request ---------------------------------------

def test_none_request_is_always_accepted_for_every_registered_model():
    # The legacy default every current shipped configuration uses; must
    # never be rejected regardless of capability.
    for name in ("lorenz", "issm", "icepack", "flowline", "unregistered-model"):
        validate_ranks_per_model_request(name, None)  # must not raise


def test_explicit_one_is_always_accepted():
    for name in ("lorenz", "issm", "icepack", "flowline"):
        validate_ranks_per_model_request(name, 1)  # must not raise


def test_explicit_multi_rank_request_rejected_for_non_supporting_model():
    with pytest.raises(ValueError, match="supports_multi_rank_per_model"):
        validate_ranks_per_model_request("lorenz", 2)
    with pytest.raises(ValueError):
        validate_ranks_per_model_request("issm", 3)
    with pytest.raises(ValueError):
        validate_ranks_per_model_request("flowline", 2)


def test_explicit_multi_rank_request_accepted_for_supporting_model():
    register_model_capabilities("unit-test-multirank-model", supports_multi_rank_per_model=True)
    validate_ranks_per_model_request("unit-test-multirank-model", 4)  # must not raise


def test_icepack_multirank_request_accepted():
    # Stage 4C (fifth continuation): the nurge-construction decomposition
    # bug that was the sole blocker is fixed (physical x-coordinate-based
    # nudging) and validated (see model_capabilities.py's icepack notes
    # and test_icepack_physical_nudge.py) -- icepack now genuinely
    # supports ranks_per_model > 1.
    validate_ranks_per_model_request("icepack", 2)  # must not raise


def test_unregistered_model_multi_rank_request_rejected():
    with pytest.raises(ValueError):
        validate_ranks_per_model_request("some-future-model-nobody-registered", 2)


def test_malformed_request_does_not_raise_here_leaves_it_to_plan_resources():
    # Garbage input (non-numeric string) is not this function's job to
    # validate -- plan_resources() raises its own clear error for it.
    validate_ranks_per_model_request("lorenz", "not-a-number")  # must not raise


# --- validate_stochastic_method_for_decomposition ----------------------------
# Stage 4C, Part 17: ranks_per_model > 1 combined with an active stochastic
# random field (nonzero sig_Q) and random_field_method="fft" is scientifically
# ambiguous -- FFT's covariance is defined over flat state-vector array-index
# distance, not physical coordinates, so it is not guaranteed
# decomposition-invariant. This must fail loudly, but ONLY in that exact
# combination -- never for a deterministic run, never for ranks_per_model=1,
# and never for the decomposition-invariant "graph" method.

def test_fft_plus_stochastic_plus_multirank_is_rejected():
    with pytest.raises(ValueError, match="scientifically ambiguous"):
        validate_stochastic_method_for_decomposition(
            "icepack", 2,
            {"sig_Q": [10, 4, 1.75, 0.14], "random_field_method": "fft"},
        )


def test_deterministic_multirank_run_not_rejected_even_though_fft_is_default():
    # A no-noise run must never be rejected merely because 'fft' is the
    # configured default random_field_method -- the guard must only fire
    # when stochastic fields are actually active.
    validate_stochastic_method_for_decomposition(
        "icepack", 2,
        {"sig_Q": [0, 0, 0, 0], "random_field_method": "fft"},
    )  # must not raise


def test_graph_method_never_rejected_regardless_of_sig_q():
    validate_stochastic_method_for_decomposition(
        "icepack", 2,
        {"sig_Q": [10, 4, 1.75, 0.14], "random_field_method": "graph"},
    )  # must not raise


def test_ranks_per_model_one_never_rejected_regardless_of_method():
    validate_stochastic_method_for_decomposition(
        "icepack", 1,
        {"sig_Q": [10, 4, 1.75, 0.14], "random_field_method": "fft"},
    )  # must not raise


def test_non_distributed_state_model_never_rejected_regardless_of_method():
    # Regression test: a model with supports_distributed_state=False
    # (e.g. lorenz -- a 3-variable ODE with no spatial decomposition)
    # has no FFT-index-space-vs-physical-coordinate ambiguity even when
    # ranks_per_model > 1 for it (e.g. Lorenz96's real oversubscribed
    # P>Nens resource-planning cases,
    # test_lorenz96_mode2_p_nens_matrix.py) -- an earlier version of this
    # guard fired for lorenz too and broke that already-passing real-MPI
    # regression suite.
    for name in ("lorenz", "issm", "flowline"):
        validate_stochastic_method_for_decomposition(
            name, 4,
            {"sig_Q": [10, 4, 1.75, 0.14], "random_field_method": "fft"},
        )  # must not raise


def test_numpy_array_sig_q_does_not_crash_the_guard():
    # Regression test: an earlier version of this guard used
    # `icesee_kwargs.get("sig_Q", []) or []`, which raises
    # "ValueError: The truth value of an array with more than one
    # element is ambiguous" when sig_Q is a numpy array (as
    # icesee_kwargs actually carries it in real runs, e.g. Lorenz96's
    # P=8/10/16,Nens=4 real-MPI matrix) rather than a plain list.
    import numpy as np

    with pytest.raises(ValueError, match="scientifically ambiguous"):
        validate_stochastic_method_for_decomposition(
            "icepack", 2,
            {"sig_Q": np.array([10.0, 4.0, 1.75, 0.14]), "random_field_method": "fft"},
        )
    validate_stochastic_method_for_decomposition(
        "icepack", 2,
        {"sig_Q": np.array([0.0, 0.0, 0.0, 0.0]), "random_field_method": "fft"},
    )  # must not raise


def test_escape_hatch_bypasses_rejection():
    validate_stochastic_method_for_decomposition(
        "icepack", 2,
        {
            "sig_Q": [10, 4, 1.75, 0.14],
            "random_field_method": "fft",
            "acknowledge_index_space_stochastic_multirank": True,
        },
    )  # must not raise
