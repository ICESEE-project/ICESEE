"""Firedrake field bridge for ICESEE execution mode 3.

The bridge is intentionally small: Firedrake retains its PETSc-distributed
vectors and ICESEE copies only the degrees of freedom owned by this rank.  No
complete field, member, or ensemble is gathered here.  Firedrake is imported
only when a halo exchange is requested, keeping the core execution-mode-3
contracts independent of a particular model stack.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Mapping

import numpy as np


# ICESEE.-prefixed (not "from src...."): distributed_fields.py/distributed_
# native_runtime.py each do their OWN relative "from .distributed_adapter
# import ..." internally. A relative import resolves against the importing
# module's own __package__, so importing this pair via the bare "src...."
# path here (while mode3_runner.py/distributed_native_runtime.py import
# the REST of this stack via the "ICESEE.src...." path) would bind their
# internal DistributedStateLayout/DistributedBlockStateLayout classes to a
# SEPARATE module identity from the ones distributed_native_runtime.py
# checks pool.layout against with isinstance() -- confirmed directly: this
# reproduced a real "native member returned an unsupported distributed
# layout" TypeError in initialize_native_member_pool during a real mode-3
# run, since the isinstance() checks there legitimately failed against a
# same-named-but-different class object.
from ICESEE.src.parallelization.distributed_fields import DistributedFieldRegistry
from ICESEE.src.parallelization.distributed_native_runtime import NativeDistributedMember


@dataclass(frozen=True)
class FiredrakeScalarOwnership:
    """Stable global ownership of scalar nodal rows on one MPI rank."""

    global_size: int
    owned_start: int
    owned_stop: int

    def __post_init__(self) -> None:
        size = int(self.global_size)
        start = int(self.owned_start)
        stop = int(self.owned_stop)
        if size < 0 or not 0 <= start <= stop <= size:
            raise ValueError("invalid Firedrake scalar ownership")
        object.__setattr__(self, "global_size", size)
        object.__setattr__(self, "owned_start", start)
        object.__setattr__(self, "owned_stop", stop)

    @property
    def owned_size(self) -> int:
        return self.owned_stop - self.owned_start


def scalar_ownership(function: Any) -> FiredrakeScalarOwnership:
    """Read PETSc ownership metadata from a scalar Firedrake function.

    A scalar function is used deliberately.  A vector-space layout may count
    components in its PETSc ownership range, whereas ICESEE stores velocity
    components as separate variable-major blocks.
    """

    try:
        layout_vec = function.function_space().dof_dset.layout_vec
        start, stop = layout_vec.getOwnershipRange()
        size = layout_vec.getSize()
    except (AttributeError, TypeError) as error:
        raise TypeError(
            "function does not expose Firedrake/PETSc scalar ownership metadata"
        ) from error
    ownership = FiredrakeScalarOwnership(size, start, stop)
    values = np.asarray(function.dat.data_ro)
    if values.ndim != 1 or values.size != ownership.owned_size:
        raise ValueError(
            "scalar Firedrake data shape is inconsistent with PETSc ownership"
        )
    return ownership


class FiredrakeOwnedField:
    """One scalar field, or one component of a vector Firedrake function."""

    def __init__(
        self,
        name: str,
        function: Any,
        ownership: FiredrakeScalarOwnership,
        *,
        component: int | None = None,
    ) -> None:
        self.name = str(name)
        self.function = function
        self.global_size = ownership.global_size
        self.owned_start = ownership.owned_start
        self.owned_stop = ownership.owned_stop
        self.component = None if component is None else int(component)
        self.synchronization_key = id(function.dat)
        self._owned_view(read_only=True)  # fail before entering a long run

    def _owned_view(self, *, read_only: bool) -> np.ndarray:
        attribute = "data_ro" if read_only else "data_wo"
        data = getattr(self.function.dat, attribute, None)
        if data is None and not read_only:
            data = getattr(self.function.dat, "data", None)
        if data is None:
            raise TypeError(f"Firedrake Dat does not expose {attribute}")
        array = np.asarray(data)
        if self.component is None:
            if array.ndim != 1:
                raise ValueError(
                    f"field {self.name!r} requires scalar Firedrake data"
                )
            view = array
        else:
            if array.ndim != 2 or not 0 <= self.component < array.shape[1]:
                raise ValueError(
                    f"field {self.name!r} has invalid vector component "
                    f"{self.component} for shape {array.shape}"
                )
            view = array[:, self.component]
        expected = self.owned_stop - self.owned_start
        if view.size != expected:
            raise ValueError(
                f"field {self.name!r} owns {view.size} rows; expected {expected}"
            )
        return view

    def read_owned(self) -> np.ndarray:
        return np.ascontiguousarray(self._owned_view(read_only=True))

    def write_owned(self, values: np.ndarray) -> None:
        target = self._owned_view(read_only=False)
        source = np.asarray(values)
        if source.ndim != 1 or source.size != target.size:
            raise ValueError(
                f"field {self.name!r} update must have shape ({target.size},)"
            )
        target[:] = source

    def synchronize_ghosts(self) -> None:
        """Refresh Firedrake halos after owned values have been updated."""

        try:
            from pyop2 import op2
        except ImportError as error:  # pragma: no cover - needs Firedrake runtime
            raise RuntimeError("Firedrake halo exchange requires pyop2") from error
        dat = self.function.dat
        dat.global_to_local_begin(op2.READ)
        dat.global_to_local_end(op2.READ)


def build_icepack_field_registry(
    *,
    thickness: Any,
    velocity: Any,
    surface: Any | None = None,
    basal_melt: Any | None = None,
    layout_id: str = "icepack-firedrake-state-v1",
) -> DistributedFieldRegistry:
    """Build Icepack's variable-major distributed state registry.

    Thickness supplies the scalar nodal ownership used by the two velocity
    components.  Other scalar fields retain their own PETSc ownership, which
    allows future mixed-space adapters without changing the core registry.
    """

    nodal = scalar_ownership(thickness)
    fields = [
        FiredrakeOwnedField("h", thickness, nodal),
        FiredrakeOwnedField("u", velocity, nodal, component=0),
        FiredrakeOwnedField("v", velocity, nodal, component=1),
    ]
    if surface is not None:
        fields.append(
            FiredrakeOwnedField("s", surface, scalar_ownership(surface))
        )
    if basal_melt is not None:
        fields.append(
            FiredrakeOwnedField(
                "basal_melt_field", basal_melt, scalar_ownership(basal_melt)
            )
        )
    return DistributedFieldRegistry(fields, layout_id=layout_id)


def _copy_owned_function(target: Any, source: Any, *, name: str) -> None:
    """Copy one model-native function without gathering global values.

    Firedrake's ``Dat`` arrays contain the values owned by the calling spatial
    rank.  Requiring identical local shapes prevents an application callback
    from accidentally handing ICESEE a replicated or globally gathered array.
    Halo synchronization is deliberately deferred until every field has been
    copied, so vector components sharing one ``Dat`` are exchanged only once.
    """

    target_data = getattr(target.dat, "data_wo", None)
    if target_data is None:
        target_data = getattr(target.dat, "data", None)
    source_data = getattr(source.dat, "data_ro", None)
    if target_data is None or source_data is None:
        raise TypeError(f"Icepack field {name!r} does not expose native Dat values")
    target_array = np.asarray(target_data)
    source_array = np.asarray(source_data)
    if target_array.shape != source_array.shape:
        raise ValueError(
            f"Icepack field {name!r} changed local ownership shape from "
            f"{target_array.shape} to {source_array.shape}"
        )
    target_array[...] = source_array


@dataclass
class IcepackForecastFields:
    """Native fields returned by one Icepack forecast callback."""

    thickness: Any
    velocity: Any
    surface: Any
    basal_melt: Any | None = None


@dataclass
class IcepackNativeState:
    """Persistent distributed Icepack state for one ensemble member.

    ``model_context`` holds application-owned objects such as the mesh, solver,
    bed, accumulation, and boundary data.  ICESEE never serializes or gathers
    those objects.  Only fields registered in ``registry`` enter analysis and
    checkpoint snapshots.
    """

    thickness: Any
    velocity: Any
    surface: Any
    basal_melt: Any | None = None
    model_context: dict[str, Any] = field(default_factory=dict)
    layout_id: str = "icepack-firedrake-state-v1"
    registry: DistributedFieldRegistry = field(init=False)

    def __post_init__(self) -> None:
        self.model_context = dict(self.model_context)
        self.registry = build_icepack_field_registry(
            thickness=self.thickness,
            velocity=self.velocity,
            surface=self.surface,
            basal_melt=self.basal_melt,
            layout_id=self.layout_id,
        )

    def update_from_forecast(self, result: IcepackForecastFields) -> None:
        """Copy forecast outputs into the persistent distributed functions."""

        if not isinstance(result, IcepackForecastFields):
            raise TypeError("Icepack forecast must return IcepackForecastFields")
        _copy_owned_function(self.thickness, result.thickness, name="h")
        _copy_owned_function(self.velocity, result.velocity, name="velocity")
        _copy_owned_function(self.surface, result.surface, name="s")
        if self.basal_melt is not None:
            if result.basal_melt is None:
                raise ValueError(
                    "forecast omitted basal_melt for a jointly estimated member"
                )
            _copy_owned_function(
                self.basal_melt, result.basal_melt, name="basal_melt_field"
            )
        elif result.basal_melt is not None:
            raise ValueError(
                "forecast returned basal_melt but it is absent from the member layout"
            )
        self.registry.synchronize_ghosts()


class IcepackNativeAdapter:
    """Callback adapter from an Icepack application to mode-3 native runtime.

    Applications keep their physics and observation logic in their existing
    modules.  This adapter only enforces the distributed lifecycle: initialize
    persistent native fields, forecast them without global packing, evaluate
    locally owned observations, and restore physical constraints after an
    analysis.  A forecast may mutate the persistent state and return ``None``,
    or return ``IcepackForecastFields`` when Icepack creates new Functions.
    """

    def __init__(
        self,
        *,
        initialize_member: Callable[..., IcepackNativeState],
        forecast_member: Callable[..., IcepackForecastFields | None],
        observe_member: Callable[..., np.ndarray] | None = None,
        finalize_analysis: Callable[..., None] | None = None,
        inverse_member: Callable[..., None] | None = None,
        restore_checkpoint: Callable[..., None] | None = None,
        reactivate_member: Callable[..., IcepackNativeState] | None = None,
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
        if reactivate_member is not None and not callable(reactivate_member):
            raise TypeError("reactivate_member must be callable")
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
        # Bounded-memory round scheduling (StreamingNativeDistributedMemberPool,
        # src/parallelization/distributed_streaming_runtime.py) duck-types
        # this attribute's presence to decide whether it may treat this
        # adapter as stream-capable -- only bound when the application
        # supplies a reactivation callback, exactly like the two optional
        # protocols above.
        if reactivate_member is not None:
            self.reactivate_native_member = self._native_reactivate_callback(
                reactivate_member
            )

    @staticmethod
    def _native_inversion_callback(callback: Callable[..., None]):
        def inverse_native_member(
            member: NativeDistributedMember,
            timestep: int,
            *,
            topology: Any,
            icesee_kwargs: Mapping[str, Any],
        ) -> None:
            state = IcepackNativeAdapter._state(member)
            result = callback(
                state,
                int(timestep),
                topology=topology,
                icesee_kwargs=icesee_kwargs,
            )
            if result is not None:
                raise TypeError("Icepack inversion callback must return None")
            state.registry.synchronize_ghosts()

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
            state = IcepackNativeAdapter._state(member)
            result = callback(
                state,
                checkpoint,
                topology=topology,
                icesee_kwargs=icesee_kwargs,
            )
            if result is not None:
                raise TypeError("Icepack checkpoint callback must return None")
            state.registry.synchronize_ghosts()

        return restore_native_checkpoint

    @staticmethod
    def _native_reactivate_callback(callback: Callable[..., IcepackNativeState]):
        def reactivate_native_member(
            member_id: int,
            packed_state: np.ndarray,
            *,
            topology: Any,
            icesee_kwargs: Mapping[str, Any],
        ) -> NativeDistributedMember:
            state = callback(
                int(member_id),
                np.asarray(packed_state),
                topology=topology,
                icesee_kwargs=icesee_kwargs,
            )
            if not isinstance(state, IcepackNativeState):
                raise TypeError("Icepack reactivator must return IcepackNativeState")
            return NativeDistributedMember(int(member_id), state.registry, state)

        return reactivate_native_member

    @staticmethod
    def _state(member: NativeDistributedMember) -> IcepackNativeState:
        state = member.model_context
        if not isinstance(state, IcepackNativeState):
            raise TypeError("native member does not contain an IcepackNativeState")
        return state

    def initialize_native_member(
        self,
        member_id: int,
        *,
        topology: Any,
        icesee_kwargs: Mapping[str, Any],
    ) -> NativeDistributedMember:
        state = self._initialize_member(
            int(member_id), topology=topology, icesee_kwargs=icesee_kwargs
        )
        if not isinstance(state, IcepackNativeState):
            raise TypeError("Icepack initializer must return IcepackNativeState")
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
        else:
            state.registry.synchronize_ghosts()

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
                "Icepack observation callback must return one value per local row"
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
                raise TypeError("Icepack analysis finalizer must return None")
        state.registry.synchronize_ghosts()
