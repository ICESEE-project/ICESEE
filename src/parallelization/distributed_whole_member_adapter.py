"""Generic mode-3 native adapter for whole-member-only model boundaries.

Some model integrations never expose per-rank spatial ownership to Python.
For example, a solver driven through an external process or RPC boundary
(one solver instance per ensemble member) may only ever produce or consume a
complete member's state at that boundary -- there is no partial slab, halo,
or ownership range for Python to read. Modifying that solver's own results
path purely to satisfy mode-3's ownership contract would violate the
model-agnostic design: a wrapper must work with what the model already
exposes, not require changing the model itself.

This module provides that adapter shape once, so any such model can reuse it
instead of a bespoke bridge. A "whole member" is represented as one
contiguous owned interval per spatial communicator of size one -- exactly the
single-vector case the design doc already names as valid
(``docs/execution-mode-3-design.md``: "A single-vector adapter may describe
one contiguous owned interval"). The adapter therefore satisfies
``NativeDistributedModelAdapter`` and gains restart/checkpoint/member-wise-
inversion lifecycle reuse from ``distributed_native_runtime.py``, but it does
NOT satisfy mode-3's memory-scaling invariant (no rank holds a complete
member): every field's owned interval is the entire variable. Applications
built on this adapter must not claim mode-3 scalability. The design doc
explicitly permits this honestly-labeled category: "Adapters that cannot
expose distributed state remain supported by modes 0--2, but cannot claim
mode-3 scalability."
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Mapping, Sequence

import numpy as np

from .distributed_fields import DistributedFieldRegistry
from .distributed_native_runtime import NativeDistributedMember


@dataclass
class WholeMemberField:
    """A ``DistributedFieldHandle`` whose entire variable is locally owned.

    ``values`` is always one-dimensional and its size never changes after
    construction (only content changes via :meth:`write_owned`), matching the
    fixed-size global/owned layout computed once at registry construction.
    """

    name: str
    values: np.ndarray

    def __post_init__(self) -> None:
        self.name = str(self.name)
        if not self.name:
            raise ValueError("whole-member field name must be nonempty")
        self.values = np.ascontiguousarray(np.asarray(self.values).ravel())

    @property
    def global_size(self) -> int:
        return int(self.values.size)

    @property
    def owned_start(self) -> int:
        return 0

    @property
    def owned_stop(self) -> int:
        return int(self.values.size)

    def read_owned(self) -> np.ndarray:
        return self.values

    def write_owned(self, values: np.ndarray) -> None:
        array = np.ascontiguousarray(np.asarray(values).ravel())
        if array.size != self.values.size:
            raise ValueError(
                f"whole-member field {self.name!r} expected {self.values.size} "
                f"values, got {array.size}"
            )
        self.values = array

    def synchronize_ghosts(self) -> None:
        """No-op: a whole-member field has no ghost/halo region."""

        return None


def build_whole_member_field_registry(
    variables: Mapping[str, np.ndarray],
    *,
    vec_inputs: Sequence[str] | None = None,
    layout_id: str = "whole-member-state-v1",
) -> DistributedFieldRegistry:
    """Wrap a ``{variable_name: array}`` member into a field registry.

    ``vec_inputs`` fixes block ordering explicitly so it matches the caller's
    existing variable-major convention (for example ``icesee_get_index``'s
    single-rank block layout) rather than relying on incidental dict
    insertion order. Defaults to ``variables``' own key order when omitted.
    """

    order = list(vec_inputs) if vec_inputs is not None else list(variables)
    missing = [name for name in order if name not in variables]
    if missing:
        raise KeyError(f"variables is missing entries for: {missing}")
    fields = [WholeMemberField(name=name, values=variables[name]) for name in order]
    return DistributedFieldRegistry(fields, layout_id=layout_id)


@dataclass
class WholeMemberState:
    """Persistent whole-member state for a whole-member-only model boundary.

    ``model_context`` holds any application-owned bookkeeping (for example a
    server handle or ensemble ID). ICESEE never serializes or gathers those
    objects; only fields registered in ``registry`` enter analysis and
    checkpoint snapshots.
    """

    variables: dict[str, np.ndarray]
    vec_inputs: Sequence[str] | None = None
    model_context: dict[str, Any] = field(default_factory=dict)
    layout_id: str = "whole-member-state-v1"
    registry: DistributedFieldRegistry = field(init=False)

    def __post_init__(self) -> None:
        self.variables = dict(self.variables)
        self.model_context = dict(self.model_context)
        self.registry = build_whole_member_field_registry(
            self.variables,
            vec_inputs=self.vec_inputs,
            layout_id=self.layout_id,
        )

    def update_from_forecast(self, result: Mapping[str, np.ndarray]) -> None:
        """Overwrite named fields in place from a callback's returned dict.

        Unknown keys are rejected so a caller cannot silently introduce a
        variable outside the layout fixed at construction time.
        """

        for name, values in result.items():
            self.registry.field(name).write_owned(values)


class WholeMemberNativeAdapter:
    """Callback adapter from a whole-member model boundary to mode-3.

    Mirrors ``IcepackNativeAdapter``'s callback-composition shape (see
    ``applications/icepack_model/icepack_utils/_distributed_fields.py``), but
    the persistent native object is a plain ``{variable_name: numpy array}``
    member instead of application mesh objects. Applications supply callables
    that already work unmodified in modes 0-2 -- no change to the model or
    its own state boundary is required. A forecast/inversion callback may
    mutate the persistent state directly and return ``None``, or return a
    ``{variable_name: array}`` dict to overwrite named fields.
    """

    def __init__(
        self,
        *,
        initialize_member: Callable[..., WholeMemberState],
        forecast_member: Callable[..., Mapping[str, np.ndarray] | None],
        observe_member: Callable[..., np.ndarray] | None = None,
        finalize_analysis: Callable[..., None] | None = None,
        inverse_member: Callable[..., Mapping[str, np.ndarray] | None] | None = None,
        restore_checkpoint: Callable[..., None] | None = None,
    ) -> None:
        for name, callback in (
            ("initialize_member", initialize_member),
            ("forecast_member", forecast_member),
        ):
            if not callable(callback):
                raise TypeError(f"{name} must be callable")
        if observe_member is not None and not callable(observe_member):
            raise TypeError("observe_member must be callable")
        if finalize_analysis is not None and not callable(finalize_analysis):
            raise TypeError("finalize_analysis must be callable")
        if inverse_member is not None and not callable(inverse_member):
            raise TypeError("inverse_member must be callable")
        if restore_checkpoint is not None and not callable(restore_checkpoint):
            raise TypeError("restore_checkpoint must be callable")
        self._initialize_member = initialize_member
        self._forecast_member = forecast_member
        self._observe_member = observe_member
        self._finalize_analysis = finalize_analysis
        # Optional protocols exist only when the application supplies them, so
        # capability validation cannot confuse a placeholder with a working
        # distributed inversion or restart implementation.
        if inverse_member is not None:
            self.inverse_native_member = self._native_inversion_callback(
                inverse_member
            )
        if restore_checkpoint is not None:
            self.restore_native_checkpoint = self._native_restore_callback(
                restore_checkpoint
            )

    @staticmethod
    def _native_inversion_callback(
        callback: Callable[..., Mapping[str, np.ndarray] | None]
    ):
        def inverse_native_member(
            member: NativeDistributedMember,
            timestep: int,
            *,
            topology: Any,
            icesee_kwargs: Mapping[str, Any],
        ) -> None:
            state = WholeMemberNativeAdapter._state(member)
            result = callback(
                state,
                int(timestep),
                topology=topology,
                icesee_kwargs=icesee_kwargs,
            )
            if result is not None:
                state.update_from_forecast(result)

        return inverse_native_member

    @staticmethod
    def _native_restore_callback(callback: Callable[..., None]):
        def restore_native_checkpoint(
            member: NativeDistributedMember,
            checkpoint: Any,
            *,
            topology: Any,
            icesee_kwargs: Mapping[str, Any],
        ) -> None:
            state = WholeMemberNativeAdapter._state(member)
            result = callback(
                state,
                checkpoint,
                topology=topology,
                icesee_kwargs=icesee_kwargs,
            )
            if result is not None:
                raise TypeError("checkpoint restore callback must return None")

        return restore_native_checkpoint

    @staticmethod
    def _state(member: NativeDistributedMember) -> WholeMemberState:
        state = member.model_context
        if not isinstance(state, WholeMemberState):
            raise TypeError("native member does not contain a WholeMemberState")
        return state

    def initialize_native_member(
        self,
        member_id: int,
        *,
        topology: Any,
        icesee_kwargs: Mapping[str, Any],
    ) -> NativeDistributedMember:
        if int(topology.spatial_ranks) != 1:
            raise ValueError(
                "WholeMemberNativeAdapter requires exactly one spatial rank "
                "per ensemble member; this model exposes no finer-grained "
                "ownership than a complete member"
            )
        state = self._initialize_member(
            int(member_id), topology=topology, icesee_kwargs=icesee_kwargs
        )
        if not isinstance(state, WholeMemberState):
            raise TypeError("initializer must return WholeMemberState")
        return NativeDistributedMember(int(member_id), state.registry, state)

    def forecast_native_member(
        self,
        member: NativeDistributedMember,
        timestep: int,
        *,
        topology: Any,
        icesee_kwargs: Mapping[str, Any],
    ) -> None:
        state = self._state(member)
        result = self._forecast_member(
            state,
            int(timestep),
            topology=topology,
            icesee_kwargs=icesee_kwargs,
        )
        if result is not None:
            state.update_from_forecast(result)

    def observe_native_member(
        self,
        member: NativeDistributedMember,
        observation_rows: np.ndarray,
        *,
        topology: Any,
        icesee_kwargs: Mapping[str, Any],
    ) -> np.ndarray:
        rows = np.asarray(observation_rows, dtype=np.int64).ravel()
        state = self._state(member)
        if self._observe_member is None:
            values = state.registry.observe_owned_rows(rows)
        else:
            values = np.asarray(
                self._observe_member(
                    state,
                    rows,
                    topology=topology,
                    icesee_kwargs=icesee_kwargs,
                )
            )
        if values.ndim != 1 or values.size != rows.size:
            raise ValueError(
                "observation callback must return one value per local row"
            )
        return np.ascontiguousarray(values)

    def finalize_native_analysis(
        self,
        member: NativeDistributedMember,
        forecast_owned: np.ndarray,
        timestep: int,
        *,
        topology: Any,
        icesee_kwargs: Mapping[str, Any],
    ) -> None:
        state = self._state(member)
        if self._finalize_analysis is not None:
            result = self._finalize_analysis(
                state,
                np.asarray(forecast_owned),
                int(timestep),
                topology=topology,
                icesee_kwargs=icesee_kwargs,
            )
            if result is not None:
                raise TypeError("analysis finalizer must return None")
