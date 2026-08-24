from __future__ import annotations

import sys
import types

import numpy as np
import pytest

from applications.icepack_model.icepack_utils._distributed_fields import (
    FiredrakeOwnedField,
    FiredrakeScalarOwnership,
    IcepackForecastFields,
    IcepackNativeAdapter,
    IcepackNativeState,
    build_icepack_field_registry,
    scalar_ownership,
)
from src.parallelization.distributed_native_runtime import (
    validate_native_distributed_adapter,
)


class _Vec:
    def __init__(self, size, ownership):
        self._size = size
        self._ownership = ownership

    def getSize(self):
        return self._size

    def getOwnershipRange(self):
        return self._ownership


class _Space:
    def __init__(self, size, ownership):
        self.dof_dset = types.SimpleNamespace(layout_vec=_Vec(size, ownership))


class _Dat:
    def __init__(self, values):
        self.values = np.asarray(values, dtype=float)
        self.begin_count = 0
        self.end_count = 0

    @property
    def data_ro(self):
        return self.values

    @property
    def data_wo(self):
        return self.values

    def global_to_local_begin(self, access):
        assert access == "READ"
        self.begin_count += 1

    def global_to_local_end(self, access):
        assert access == "READ"
        self.end_count += 1


class _Function:
    def __init__(self, values, *, size, ownership):
        self.dat = _Dat(values)
        self._space = _Space(size, ownership)

    def function_space(self):
        return self._space


def test_scalar_ownership_uses_petsc_layout():
    function = _Function([1, 2, 3], size=11, ownership=(4, 7))
    assert scalar_ownership(function) == FiredrakeScalarOwnership(11, 4, 7)


def test_registry_packs_native_scalar_and_vector_components(monkeypatch):
    h = _Function([1, 2], size=6, ownership=(2, 4))
    u = _Function([[10, 20], [11, 21]], size=12, ownership=(4, 8))
    s = _Function([30, 31], size=6, ownership=(2, 4))
    registry = build_icepack_field_registry(
        thickness=h, velocity=u, surface=s, layout_id="pig-native-v1"
    )
    np.testing.assert_array_equal(
        registry.pack_owned(), [1, 2, 10, 11, 20, 21, 30, 31]
    )
    assert registry.layout.block_names == ("h", "u", "v", "s")
    assert registry.layout.global_size == 24

    monkeypatch.setitem(sys.modules, "pyop2", types.SimpleNamespace(
        op2=types.SimpleNamespace(READ="READ")
    ))
    registry.unpack_owned(np.arange(8, dtype=float))
    np.testing.assert_array_equal(h.dat.values, [0, 1])
    np.testing.assert_array_equal(u.dat.values, [[2, 4], [3, 5]])
    np.testing.assert_array_equal(s.dat.values, [6, 7])
    assert u.dat.begin_count == 1
    assert u.dat.end_count == 1


def test_field_rejects_component_or_ownership_mismatch():
    ownership = FiredrakeScalarOwnership(9, 3, 5)
    scalar = _Function([1, 2], size=9, ownership=(3, 5))
    with pytest.raises(ValueError, match="invalid vector component"):
        FiredrakeOwnedField("u", scalar, ownership, component=0)

    vector = _Function([[1, 2]], size=4, ownership=(1, 3))
    with pytest.raises(ValueError, match="expected 2"):
        FiredrakeOwnedField("u", vector, ownership, component=0)


def _native_state(*, offset=0.0, include_melt=True):
    kwargs = dict(size=6, ownership=(2, 4))
    return IcepackNativeState(
        thickness=_Function(np.array([1, 2]) + offset, **kwargs),
        velocity=_Function(
            np.array([[10, 20], [11, 21]]) + offset,
            size=12,
            ownership=(4, 8),
        ),
        surface=_Function(np.array([30, 31]) + offset, **kwargs),
        basal_melt=(
            _Function(np.array([40, 41]) + offset, **kwargs)
            if include_melt else None
        ),
        model_context={"solver": object()},
        layout_id="pig-native-test-v1",
    )


def _forecast_fields(*, offset=0.0, include_melt=True):
    kwargs = dict(size=6, ownership=(2, 4))
    return IcepackForecastFields(
        thickness=_Function(np.array([101, 102]) + offset, **kwargs),
        velocity=_Function(
            np.array([[110, 120], [111, 121]]) + offset,
            size=12,
            ownership=(4, 8),
        ),
        surface=_Function(np.array([130, 131]) + offset, **kwargs),
        basal_melt=(
            _Function(np.array([140, 141]) + offset, **kwargs)
            if include_melt else None
        ),
    )


def test_native_state_copies_forecast_into_persistent_fields(monkeypatch):
    monkeypatch.setitem(
        sys.modules,
        "pyop2",
        types.SimpleNamespace(op2=types.SimpleNamespace(READ="READ")),
    )
    state = _native_state()
    original_functions = (
        state.thickness,
        state.velocity,
        state.surface,
        state.basal_melt,
    )
    state.update_from_forecast(_forecast_fields())

    assert (
        state.thickness,
        state.velocity,
        state.surface,
        state.basal_melt,
    ) == original_functions
    np.testing.assert_array_equal(
        state.registry.pack_owned(),
        [101, 102, 110, 111, 120, 121, 130, 131, 140, 141],
    )
    assert state.velocity.dat.begin_count == 1
    assert state.velocity.dat.end_count == 1


def test_native_state_rejects_changed_local_ownership(monkeypatch):
    monkeypatch.setitem(
        sys.modules,
        "pyop2",
        types.SimpleNamespace(op2=types.SimpleNamespace(READ="READ")),
    )
    state = _native_state(include_melt=False)
    forecast = _forecast_fields(include_melt=False)
    forecast.thickness = _Function([1, 2, 3], size=7, ownership=(2, 5))
    with pytest.raises(ValueError, match="changed local ownership shape"):
        state.update_from_forecast(forecast)


def test_callback_adapter_preserves_native_member_lifecycle(monkeypatch):
    monkeypatch.setitem(
        sys.modules,
        "pyop2",
        types.SimpleNamespace(op2=types.SimpleNamespace(READ="READ")),
    )
    events = []

    def initialize(member_id, **kwargs):
        events.append(("initialize", member_id))
        return _native_state(offset=member_id)

    def forecast(state, timestep, **kwargs):
        events.append(("forecast", timestep, state.model_context["solver"]))
        return _forecast_fields(offset=timestep)

    def observe(state, rows, **kwargs):
        packed = state.registry.pack_owned()
        return packed[rows]

    def finalize(state, forecast_owned, timestep, **kwargs):
        events.append(("finalize", timestep, forecast_owned.copy()))

    adapter = IcepackNativeAdapter(
        initialize_member=initialize,
        forecast_member=forecast,
        observe_member=observe,
        finalize_analysis=finalize,
    )
    member = adapter.initialize_native_member(
        3, topology=object(), icesee_kwargs={}
    )
    persistent_thickness = member.model_context.thickness
    adapter.forecast_native_member(
        member, 5, topology=object(), icesee_kwargs={}
    )
    assert member.model_context.thickness is persistent_thickness
    np.testing.assert_array_equal(
        adapter.observe_native_member(
            member, np.array([0, 4, 9]), topology=object(), icesee_kwargs={}
        ),
        [106, 125, 146],
    )
    forecast_owned = member.pack_owned()
    adapter.finalize_native_analysis(
        member, forecast_owned, 5, topology=object(), icesee_kwargs={}
    )
    assert events[0] == ("initialize", 3)
    assert events[1][0:2] == ("forecast", 5)
    assert events[2][0:2] == ("finalize", 5)


def test_callback_adapter_defaults_to_owned_identity_observations():
    adapter = IcepackNativeAdapter(
        initialize_member=lambda member_id, **kwargs: _native_state(
            offset=member_id
        ),
        forecast_member=lambda state, timestep, **kwargs: None,
    )
    member = adapter.initialize_native_member(
        3, topology=object(), icesee_kwargs={}
    )

    # Global variable-major rows: h[2], v[2], basal_melt[3].
    np.testing.assert_array_equal(
        adapter.observe_native_member(
            member,
            np.asarray([2, 14, 27]),
            topology=object(),
            icesee_kwargs={},
        ),
        [4, 23, 44],
    )


def test_callback_adapter_exposes_only_configured_optional_capabilities(monkeypatch):
    monkeypatch.setitem(
        sys.modules,
        "pyop2",
        types.SimpleNamespace(op2=types.SimpleNamespace(READ="READ")),
    )
    base = dict(
        initialize_member=lambda member_id, **kwargs: _native_state(
            offset=member_id
        ),
        forecast_member=lambda state, timestep, **kwargs: None,
    )
    without_optional = IcepackNativeAdapter(**base)
    assert not hasattr(without_optional, "inverse_native_member")
    assert not hasattr(without_optional, "restore_native_checkpoint")
    with pytest.raises(TypeError, match="inverse_native_member"):
        validate_native_distributed_adapter(
            without_optional, require_inversion=True
        )

    events = []

    def inverse(state, timestep, **kwargs):
        state.registry.unpack_owned(state.registry.pack_owned() + 2.0)
        events.append(("inverse", timestep))

    def restore(state, checkpoint, **kwargs):
        state.model_context["checkpoint"] = checkpoint
        events.append(("restore", checkpoint))

    adapter = IcepackNativeAdapter(
        **base,
        inverse_member=inverse,
        restore_checkpoint=restore,
    )
    validate_native_distributed_adapter(adapter, require_inversion=True)
    member = adapter.initialize_native_member(
        0, topology=object(), icesee_kwargs={}
    )
    before = member.pack_owned().copy()
    adapter.inverse_native_member(
        member, 4, topology=object(), icesee_kwargs={}
    )
    np.testing.assert_array_equal(member.pack_owned(), before + 2.0)
    adapter.restore_native_checkpoint(
        member, "checkpoint-4", topology=object(), icesee_kwargs={}
    )
    assert events == [("inverse", 4), ("restore", "checkpoint-4")]
