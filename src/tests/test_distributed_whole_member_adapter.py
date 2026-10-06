"""Unit tests for the generic whole-member-only mode-3 native adapter.

These tests are pure Python and do not require MATLAB, ISSM, or any other
external solver: they validate the wiring contract itself using simple
in-memory callables, the same way ``test_icepack_distributed_fields.py``
validates the Firedrake bridge with fake Firedrake objects.
"""

from __future__ import annotations

import types

import numpy as np
import pytest

from src.parallelization.distributed_fields import DistributedFieldRegistry
from src.parallelization.distributed_native_runtime import NativeDistributedMember
from src.parallelization.distributed_whole_member_adapter import (
    WholeMemberField,
    WholeMemberNativeAdapter,
    WholeMemberState,
    build_whole_member_field_registry,
)


def _topology(spatial_ranks=1):
    return types.SimpleNamespace(spatial_ranks=spatial_ranks)


# --- WholeMemberField -------------------------------------------------------


def test_whole_member_field_owns_the_entire_variable():
    field = WholeMemberField(name="Thickness", values=np.array([1.0, 2.0, 3.0]))
    assert field.global_size == 3
    assert field.owned_start == 0
    assert field.owned_stop == 3
    np.testing.assert_array_equal(field.read_owned(), [1.0, 2.0, 3.0])


def test_whole_member_field_write_owned_replaces_values():
    field = WholeMemberField(name="Vx", values=np.zeros(4))
    field.write_owned(np.array([1.0, 2.0, 3.0, 4.0]))
    np.testing.assert_array_equal(field.read_owned(), [1.0, 2.0, 3.0, 4.0])


def test_whole_member_field_write_owned_rejects_size_mismatch():
    field = WholeMemberField(name="Vx", values=np.zeros(4))
    with pytest.raises(ValueError):
        field.write_owned(np.zeros(3))


def test_whole_member_field_rejects_empty_name():
    with pytest.raises(ValueError):
        WholeMemberField(name="", values=np.zeros(2))


def test_whole_member_field_synchronize_ghosts_is_a_noop():
    field = WholeMemberField(name="Vx", values=np.zeros(2))
    assert field.synchronize_ghosts() is None


# --- build_whole_member_field_registry --------------------------------------


def test_build_registry_uses_vec_inputs_order_for_blocks():
    variables = {"b": np.zeros(2), "a": np.zeros(3)}
    registry = build_whole_member_field_registry(variables, vec_inputs=["a", "b"])
    assert isinstance(registry, DistributedFieldRegistry)
    assert registry.names == ("a", "b")
    assert registry.layout.global_size == 5


def test_build_registry_defaults_to_dict_key_order():
    variables = {"a": np.zeros(2), "b": np.zeros(3)}
    registry = build_whole_member_field_registry(variables)
    assert registry.names == ("a", "b")


def test_build_registry_rejects_missing_vec_inputs_entry():
    variables = {"a": np.zeros(2)}
    with pytest.raises(KeyError):
        build_whole_member_field_registry(variables, vec_inputs=["a", "missing"])


# --- WholeMemberState ---------------------------------------------------


def test_whole_member_state_packs_owned_in_vec_inputs_order():
    state = WholeMemberState(
        variables={"Vy": np.array([3.0, 4.0]), "Vx": np.array([1.0, 2.0])},
        vec_inputs=["Vx", "Vy"],
    )
    np.testing.assert_array_equal(state.registry.pack_owned(), [1.0, 2.0, 3.0, 4.0])


def test_whole_member_state_update_from_forecast_overwrites_named_fields():
    state = WholeMemberState(variables={"Vx": np.zeros(2)}, vec_inputs=["Vx"])
    state.update_from_forecast({"Vx": np.array([5.0, 6.0])})
    np.testing.assert_array_equal(state.registry.field("Vx").read_owned(), [5.0, 6.0])


def test_whole_member_state_update_from_forecast_rejects_unknown_field():
    state = WholeMemberState(variables={"Vx": np.zeros(2)}, vec_inputs=["Vx"])
    with pytest.raises(KeyError):
        state.update_from_forecast({"unknown": np.zeros(2)})


# --- WholeMemberNativeAdapter --------------------------------------------


def _make_adapter(**overrides):
    def initialize_member(member_id, *, topology, icesee_kwargs):
        return WholeMemberState(
            variables={"Vx": np.array([float(member_id), float(member_id)])},
            vec_inputs=["Vx"],
            model_context={"ens_id": member_id},
        )

    def forecast_member(state, timestep, *, topology, icesee_kwargs):
        return {"Vx": state.registry.field("Vx").read_owned() + 1.0}

    kwargs = dict(
        initialize_member=initialize_member,
        forecast_member=forecast_member,
    )
    kwargs.update(overrides)
    return WholeMemberNativeAdapter(**kwargs)


def test_adapter_requires_callables():
    with pytest.raises(TypeError):
        WholeMemberNativeAdapter(initialize_member=None, forecast_member=lambda *a, **k: None)


def test_adapter_initialize_native_member_builds_member():
    adapter = _make_adapter()
    member = adapter.initialize_native_member(
        3, topology=_topology(), icesee_kwargs={}
    )
    assert isinstance(member, NativeDistributedMember)
    assert member.member_id == 3
    np.testing.assert_array_equal(member.fields.pack_owned(), [3.0, 3.0])


def test_adapter_initialize_native_member_rejects_multi_spatial_rank_topology():
    adapter = _make_adapter()
    with pytest.raises(ValueError):
        adapter.initialize_native_member(
            0, topology=_topology(spatial_ranks=2), icesee_kwargs={}
        )


def test_adapter_forecast_native_member_applies_returned_dict():
    adapter = _make_adapter()
    member = adapter.initialize_native_member(1, topology=_topology(), icesee_kwargs={})
    adapter.forecast_native_member(member, 0, topology=_topology(), icesee_kwargs={})
    np.testing.assert_array_equal(member.fields.pack_owned(), [2.0, 2.0])


def test_adapter_forecast_native_member_allows_in_place_none_return():
    def forecast_member(state, timestep, *, topology, icesee_kwargs):
        state.registry.field("Vx").write_owned(np.array([9.0, 9.0]))
        return None

    adapter = _make_adapter(forecast_member=forecast_member)
    member = adapter.initialize_native_member(0, topology=_topology(), icesee_kwargs={})
    adapter.forecast_native_member(member, 0, topology=_topology(), icesee_kwargs={})
    np.testing.assert_array_equal(member.fields.pack_owned(), [9.0, 9.0])


def test_adapter_observe_native_member_defaults_to_registry_lookup():
    adapter = _make_adapter()
    member = adapter.initialize_native_member(2, topology=_topology(), icesee_kwargs={})
    values = adapter.observe_native_member(
        member, np.array([0, 1]), topology=_topology(), icesee_kwargs={}
    )
    np.testing.assert_array_equal(values, [2.0, 2.0])


def test_adapter_observe_native_member_uses_custom_callback():
    def observe_member(state, rows, *, topology, icesee_kwargs):
        return np.full(rows.shape, 42.0)

    adapter = _make_adapter(observe_member=observe_member)
    member = adapter.initialize_native_member(0, topology=_topology(), icesee_kwargs={})
    values = adapter.observe_native_member(
        member, np.array([0]), topology=_topology(), icesee_kwargs={}
    )
    np.testing.assert_array_equal(values, [42.0])


def test_adapter_observe_native_member_rejects_shape_mismatch():
    def observe_member(state, rows, *, topology, icesee_kwargs):
        return np.zeros(rows.size + 1)

    adapter = _make_adapter(observe_member=observe_member)
    member = adapter.initialize_native_member(0, topology=_topology(), icesee_kwargs={})
    with pytest.raises(ValueError):
        adapter.observe_native_member(
            member, np.array([0]), topology=_topology(), icesee_kwargs={}
        )


def test_adapter_finalize_native_analysis_calls_finalizer():
    calls = []

    def finalize_analysis(state, forecast_owned, timestep, *, topology, icesee_kwargs):
        calls.append((forecast_owned.tolist(), timestep))

    adapter = _make_adapter(finalize_analysis=finalize_analysis)
    member = adapter.initialize_native_member(0, topology=_topology(), icesee_kwargs={})
    adapter.finalize_native_analysis(
        member, np.array([1.0, 2.0]), 5, topology=_topology(), icesee_kwargs={}
    )
    assert calls == [([1.0, 2.0], 5)]


def test_adapter_finalize_native_analysis_rejects_non_none_return():
    def finalize_analysis(*args, **kwargs):
        return "not none"

    adapter = _make_adapter(finalize_analysis=finalize_analysis)
    member = adapter.initialize_native_member(0, topology=_topology(), icesee_kwargs={})
    with pytest.raises(TypeError):
        adapter.finalize_native_analysis(
            member, np.zeros(2), 0, topology=_topology(), icesee_kwargs={}
        )


def test_adapter_without_inverse_member_has_no_inversion_capability():
    adapter = _make_adapter()
    assert not hasattr(adapter, "inverse_native_member")


def test_adapter_inverse_member_applies_returned_dict():
    def inverse_member(state, timestep, *, topology, icesee_kwargs):
        return {"Vx": np.array([7.0, 7.0])}

    adapter = _make_adapter(inverse_member=inverse_member)
    member = adapter.initialize_native_member(0, topology=_topology(), icesee_kwargs={})
    adapter.inverse_native_member(member, 0, topology=_topology(), icesee_kwargs={})
    np.testing.assert_array_equal(member.fields.pack_owned(), [7.0, 7.0])


def test_adapter_inverse_member_allows_in_place_none_return():
    def inverse_member(state, timestep, *, topology, icesee_kwargs):
        state.registry.field("Vx").write_owned(np.array([8.0, 8.0]))
        return None

    adapter = _make_adapter(inverse_member=inverse_member)
    member = adapter.initialize_native_member(0, topology=_topology(), icesee_kwargs={})
    adapter.inverse_native_member(member, 0, topology=_topology(), icesee_kwargs={})
    np.testing.assert_array_equal(member.fields.pack_owned(), [8.0, 8.0])


def test_adapter_without_restore_checkpoint_has_no_restart_capability():
    adapter = _make_adapter()
    assert not hasattr(adapter, "restore_native_checkpoint")


def test_adapter_restore_checkpoint_calls_callback():
    calls = []

    def restore_checkpoint(state, checkpoint, *, topology, icesee_kwargs):
        calls.append(checkpoint)

    adapter = _make_adapter(restore_checkpoint=restore_checkpoint)
    member = adapter.initialize_native_member(0, topology=_topology(), icesee_kwargs={})
    adapter.restore_native_checkpoint(
        member, {"marker": 1}, topology=_topology(), icesee_kwargs={}
    )
    assert calls == [{"marker": 1}]


def test_adapter_restore_checkpoint_rejects_non_none_return():
    def restore_checkpoint(*args, **kwargs):
        return "not none"

    adapter = _make_adapter(restore_checkpoint=restore_checkpoint)
    member = adapter.initialize_native_member(0, topology=_topology(), icesee_kwargs={})
    with pytest.raises(TypeError):
        adapter.restore_native_checkpoint(
            member, {}, topology=_topology(), icesee_kwargs={}
        )
