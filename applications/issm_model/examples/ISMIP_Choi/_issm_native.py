"""Whole-member ISSM state, forecast and inversion callbacks.

Each Python rank controls a MATLAB server and exchanges complete member fields.
ISSM's internal FEM ranks do not expose spatial ownership to this adapter.
Optional recovery restores full MATLAB model files and RNG state as well as
Python fields; it does not claim spatially partitioned state scalability.
"""

from __future__ import annotations

from typing import Any, Mapping

import numpy as np

from ICESEE.src.parallelization.distributed_whole_member_adapter import (
    WholeMemberNativeAdapter,
    WholeMemberState,
)
from ICESEE.applications.issm_model.examples.ISMIP_Choi._issm_enkf import (
    forecast_step_single,
    initialize_ensemble,
)
from ICESEE.applications.issm_model.examples.ISMIP_Choi._issm_model import run_model_inverse


def _vec_inputs(icesee_kwargs: Mapping[str, Any]) -> list[str]:
    vec_inputs = icesee_kwargs.get("vec_inputs")
    if not vec_inputs:
        raise ValueError("icesee_kwargs['vec_inputs'] is required")
    return list(vec_inputs)


def _pack_ensemble(state: WholeMemberState, vec_inputs: list[str]) -> np.ndarray:
    """Rebuild the flat variable-major array ``run_model``/``run_model_inverse``
    expect, in the same block order ``icesee_get_index`` produces for a
    single-rank member (see ``src/utils/tools.py::icesee_get_index``,
    equal-size branch: block ``i`` occupies rows
    ``[i * hdim, (i + 1) * hdim)`` in ``vec_inputs`` order).
    """

    return np.concatenate(
        [state.registry.field(name).read_owned() for name in vec_inputs]
    )


def initialize_member(
    member_id: int,
    *,
    topology: Any,
    icesee_kwargs: Mapping[str, Any],
) -> WholeMemberState:
    """Create one persistent whole-member ISSM state.

    Delegates to the existing, unmodified ``initialize_ensemble`` (which
    already drives the member's MATLAB server and reads back its ensemble
    initialization HDF5 file exactly as modes 0-2 do).
    """

    vec_inputs = _vec_inputs(icesee_kwargs)
    updated_state = initialize_ensemble(int(member_id), **dict(icesee_kwargs))
    variables = {name: np.asarray(updated_state[name]).ravel() for name in vec_inputs}
    return WholeMemberState(
        variables=variables,
        vec_inputs=vec_inputs,
        model_context={"ens_id": int(member_id)},
    )


def forecast_member(
    state: WholeMemberState,
    timestep: int,
    *,
    topology: Any,
    icesee_kwargs: Mapping[str, Any],
) -> dict[str, np.ndarray] | None:
    """Advance one persistent member via the existing forecast_step_single.

    Builds the flat ensemble array that ``forecast_step_single``/
    ``run_model`` already expect from the persistent native fields, calls the
    unmodified function, and returns its updated variables to be written back
    into those fields.
    """

    vec_inputs = _vec_inputs(icesee_kwargs)
    ensemble = _pack_ensemble(state, vec_inputs)
    call_kwargs = dict(icesee_kwargs)
    call_kwargs.update({"k": int(timestep), "ens_id": state.model_context["ens_id"]})
    updated_state = forecast_step_single(ensemble=ensemble, **call_kwargs)
    return {
        name: np.asarray(updated_state[name]).ravel()
        for name in vec_inputs
        if name in updated_state
    }


def inverse_member(
    state: WholeMemberState,
    timestep: int,
    *,
    topology: Any,
    icesee_kwargs: Mapping[str, Any],
) -> dict[str, np.ndarray] | None:
    """Apply one member-wise friction/velocity inversion via run_model_inverse.

    ``run_model_inverse`` expects ``icesee_kwargs['km']`` to already carry the
    inversion timestep, exactly as modes 0-2 set it before calling
    ``inverse_step_single``; this wiring does not invent that value.
    """

    vec_inputs = _vec_inputs(icesee_kwargs)
    ensemble = _pack_ensemble(state, vec_inputs)
    call_kwargs = dict(icesee_kwargs)
    call_kwargs.update({"ens_id": state.model_context["ens_id"]})
    updated_state = run_model_inverse(ensemble, **call_kwargs)
    if updated_state is None:
        raise RuntimeError(
            "ISSM inversion returned no result for member "
            f"{state.model_context['ens_id']} at timestep {int(timestep)}; "
            "inspect the preceding model-adapter/ISSM error."
        )
    return {
        name: np.asarray(updated_state[name]).ravel()
        for name in vec_inputs
        if name in updated_state
    }


class _RecoverableIssmAdapter(WholeMemberNativeAdapter):
    def allocate_native_member(self, member_id, *, topology, icesee_kwargs):
        from ICESEE.src.parallelization.distributed_native_runtime import NativeDistributedMember
        names = _vec_inputs(icesee_kwargs)
        nd = int(icesee_kwargs["nd"])
        if nd % len(names):
            raise ValueError("ISSM checkpoint fields must share the native vertex count")
        state = WholeMemberState(variables={name:np.zeros(nd//len(names)) for name in names},
            vec_inputs=names,model_context={"ens_id":int(member_id)})
        return NativeDistributedMember(int(member_id),state.registry,state)

    def restore_native_checkpoint(self, member, checkpoint, *, topology, icesee_kwargs):
        from ICESEE.applications.issm_model.issm_utils.recovery import restore_backend_member
        restore_backend_member(member.member_id,checkpoint,topology,icesee_kwargs)


ISMIP_CHOI_NATIVE_ADAPTER = _RecoverableIssmAdapter(
    initialize_member=initialize_member,
    forecast_member=forecast_member,
    inverse_member=inverse_member,
)
