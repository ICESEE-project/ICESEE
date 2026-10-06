"""Mode-3 native adapter wiring for the ISSM ISMIP_Choi hybrid application.

## Why this is a whole-member adapter, not a spatially partitioned one

Investigated 2026-08-24 before writing any code here (see the ICESEE
mode-3-review memory notes for the full trail). ISSM's Python/MATLAB boundary
has no per-rank ownership concept at all today, and this file does not add
one:

- Ensemble parallelism is one ICESEE MPI rank per member (``ens_id = rank``,
  enforced by ``MatlabServer._check_conditions`` in
  ``issm_utils/matlab2python/mat2py_utils.py``); each rank drives its own
  persistent MATLAB subprocess.
- ``model_nprocs`` (``md.cluster = generic(..., 'np', nprocs)`` in
  ``run_model.m``) is ISSM's own *internal* MATLAB/PETSc parallel FEM solve.
  It is invoked and torn down entirely inside one MATLAB command and is never
  exposed to Python as a spatial ownership range.
- Every Python<->MATLAB state exchange is a complete-member HDF5 file: the
  MATLAB side gathers its internal parallel solution and does one
  ``h5create``/``h5write`` of the full ``md.mesh.numberofvertices``-sized
  array per variable (``run_model.m:1731-1754``); the Python side always
  reads the whole array back (``_issm_model.py``, ``_issm_enkf.py``). No
  ``io_gather``-style partitioned-output setting is used anywhere in this
  repository's ISSM MATLAB scripts.

Per explicit project direction, this work must not modify ISSM's own model
code (``run_model.m``, or any other MATLAB script) to manufacture ownership
metadata that does not already exist, and must not change ICESEE's core
architecture. The model-agnostic answer is therefore
``WholeMemberNativeAdapter``
(``src/parallelization/distributed_whole_member_adapter.py``): a reusable,
already-unit-tested adapter for exactly this situation -- a model whose only
unmodified Python-visible boundary is a complete ensemble member. This file
wires it to ISSM's *existing, unmodified* ``initialize_ensemble``,
``forecast_step_single``, and ``run_model_inverse`` functions.

## What this does and does not achieve

This satisfies the ``NativeDistributedModelAdapter``/
``NativeDistributedInversionAdapter`` protocols, so ISSM gains restart/
checkpoint/member-wise-inversion lifecycle reuse from
``distributed_native_runtime.py`` without touching MATLAB or core ICESEE
code. It does **not** satisfy mode-3's "no rank holds a complete member"
memory-scaling invariant: every field's owned interval is the entire
variable, exactly as it already is in modes 0-2. Per
``docs/execution-mode-3-design.md``, this is the explicitly permitted
state-only category ("Adapters that cannot expose distributed state remain
supported by modes 0--2, but cannot claim mode-3 scalability") -- ISSM
cannot claim mode-3 scalability through this adapter, and ``execution_mode:
3`` remains non-selectable regardless.
"""

from __future__ import annotations

from typing import Any, Mapping

import numpy as np

from src.parallelization.distributed_whole_member_adapter import (
    WholeMemberNativeAdapter,
    WholeMemberState,
)
from applications.issm_model.examples.ISMIP_Choi._issm_enkf import (
    forecast_step_single,
    initialize_ensemble,
)
from applications.issm_model.examples.ISMIP_Choi._issm_model import run_model_inverse


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
        return None
    return {
        name: np.asarray(updated_state[name]).ravel()
        for name in vec_inputs
        if name in updated_state
    }


ISMIP_CHOI_NATIVE_ADAPTER = WholeMemberNativeAdapter(
    initialize_member=initialize_member,
    forecast_member=forecast_member,
    inverse_member=inverse_member,
)
