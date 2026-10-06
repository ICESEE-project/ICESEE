import pytest

from ICESEE.src.parallelization.distributed_mode3_registry import (
    Mode3Registration,
    get_execution_mode_3,
    register_execution_mode_3,
    registered_execution_mode_3_models,
)
from ICESEE.src.run_model_da.icesee_da_distributed import (
    icesee_model_data_assimilation_distributed,
)


@pytest.fixture(autouse=True)
def _clean_registry():
    """Ensure each test starts and ends with no test-only registrations."""
    before = set(registered_execution_mode_3_models())
    yield
    from ICESEE.src.parallelization import distributed_mode3_registry as module

    for name in set(module._REGISTRY) - before:
        del module._REGISTRY[name]


def test_register_and_get_round_trips():
    def runner(**kwargs):
        return kwargs

    register_execution_mode_3("unit-test-model", runner, notes="test only")
    registration = get_execution_mode_3("unit-test-model")

    assert isinstance(registration, Mode3Registration)
    assert registration.model_name == "unit-test-model"
    assert registration.runner is runner
    assert registration.notes == "test only"
    assert "unit-test-model" in registered_execution_mode_3_models()


def test_get_unregistered_model_returns_none():
    assert get_execution_mode_3("does-not-exist") is None
    assert get_execution_mode_3(None) is None
    assert get_execution_mode_3("") is None


def test_reregistering_same_model_name_overwrites():
    register_execution_mode_3("unit-test-model", lambda **kw: "first")
    register_execution_mode_3("unit-test-model", lambda **kw: "second")

    registration = get_execution_mode_3("unit-test-model")
    assert registration.runner() == "second"


def test_register_rejects_empty_name():
    with pytest.raises(ValueError):
        register_execution_mode_3("", lambda **kw: None)


def test_register_rejects_non_callable_runner():
    with pytest.raises(TypeError):
        register_execution_mode_3("unit-test-model", "not-callable")


def test_registered_execution_mode_3_models_is_sorted():
    register_execution_mode_3("zeta-model", lambda **kw: None)
    register_execution_mode_3("alpha-model", lambda **kw: None)

    names = registered_execution_mode_3_models()
    assert list(names) == sorted(names)
    assert "alpha-model" in names and "zeta-model" in names


def test_dispatch_raises_named_error_for_unregistered_model():
    with pytest.raises(NotImplementedError, match="never-registered-model"):
        icesee_model_data_assimilation_distributed(
            execution_mode=3, model_name="never-registered-model"
        )


def test_dispatch_invokes_registered_runner_with_kwargs():
    received = {}

    def runner(**kwargs):
        received.update(kwargs)
        return "ran"

    register_execution_mode_3("unit-test-model", runner)

    result = icesee_model_data_assimilation_distributed(
        execution_mode=3, model_name="unit-test-model", Nens=5
    )

    assert result == "ran"
    assert received["model_name"] == "unit-test-model"
    assert received["Nens"] == 5
    assert received["execution_mode"] == 3
