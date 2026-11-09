from types import SimpleNamespace
import pytest
from ICESEE.src.parallelization import distributed_mode3_registry as registry
from ICESEE.src.parallelization import native_application as runtime


def test_examples_of_same_model_do_not_share_a_runner(monkeypatch):
    monkeypatch.setattr(registry, "_REGISTRY", {})
    pig = lambda **kwargs: "pig"
    synthetic = lambda **kwargs: "synthetic"
    registry.register_execution_mode_3("icepack", pig, application="idealized_pig")
    assert registry.get_execution_mode_3("icepack", application="synthetic_ice_stream") is None
    registry.register_execution_mode_3("icepack", synthetic, application="synthetic_ice_stream")
    assert registry.get_execution_mode_3("icepack", application="idealized_pig").runner is pig
    assert registry.get_execution_mode_3("icepack", application="synthetic_ice_stream").runner is synthetic
    assert registry.registered_execution_mode_3_models() == ("icepack",)


def test_binding_reuses_shared_runtime_and_preserves_options(monkeypatch):
    monkeypatch.setattr(registry, "_REGISTRY", {})
    adapter = SimpleNamespace(initialize_native_member=lambda: None, forecast_native_member=lambda: None)
    setup = lambda options: options
    monkeypatch.setattr(runtime, "run_native_application", lambda **options: options)
    registry.register_native_application("sample", application="example", adapter=adapter, setup=setup)
    bound = registry.get_execution_mode_3("sample", application="example").runner(Nens=8, execution_mode=3)
    assert bound["adapter"] is adapter and bound["setup"] is setup
    assert bound["Nens"] == 8 and bound["execution_mode"] == 3


@pytest.mark.parametrize("option", ["localization_flag", "inversion_flag", "inference_plugin_enabled"])
def test_unsupported_scientific_policy_is_never_silently_ignored(option):
    with pytest.raises(NotImplementedError, match=option):
        runtime.run_native_application(adapter=None, setup=lambda x:x, **{option:True})


def test_missing_adapter_operations_fail_at_registration():
    with pytest.raises(TypeError, match="initialize_native_member"):
        registry.register_native_application("sample", application="example", adapter=object(), setup=lambda x:x)


def test_resume_without_validated_binding_is_rejected():
    with pytest.raises(NotImplementedError,match="restart binding"):
        runtime.run_native_application(adapter=None,setup=lambda x:x,
                                       resume_from_checkpoint=True)
