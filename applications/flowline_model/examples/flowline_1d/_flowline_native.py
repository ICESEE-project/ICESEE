# ==============================================================================
# @des: Mode-3 native adapter wiring for the flowline_1d application.
# @date: 2026-08-26
# @author: Brian Kyanjo
# ==============================================================================
"""Mode-3 native adapter wiring for the flowline_1d application.

## Why this is a whole-member adapter, not a spatially partitioned one

flowline_1d solves the entire discretized ice-flow state (``h``, ``u``,
``xg``, size ``2*NX+1``) as one coupled nonlinear system per ensemble
member, via ``scipy.optimize.root`` over a JAX-autodiff Jacobian
(``_flowline_model.py::run_model``/``initialize_model``). There is no
per-rank/spatial ownership boundary anywhere in that solve -- every call
consumes and produces the complete member vector in one Python call, purely
in-process (no external server, no MPI/collective operation of any kind
inside the model itself). This is architecturally the same
"whole-member-only" category as ISSM
(``applications/issm_model/examples/ISMIP_Choi/_issm_native.py``), even
simpler since there is no persistent external process to manage per member.

Per the same project direction documented in ISSM's native adapter, this
file does not modify ``_flowline_model.py``/``_flowline_enkf.py`` to invent
finer-grained ownership that does not exist; it wires the existing,
unmodified ``initialize_ensemble``/``forecast_step_single`` functions (used
unchanged by modes 0-2) into ``WholeMemberNativeAdapter``
(``src/parallelization/distributed_whole_member_adapter.py``).

## What this does and does not achieve

Satisfies ``NativeDistributedModelAdapter``, so flowline_1d gains restart/
checkpoint lifecycle reuse from ``distributed_native_runtime.py``. It does
**not** satisfy mode-3's "no rank holds a complete member" memory-scaling
invariant -- exactly the same honestly-labeled tradeoff as ISSM's adapter;
see that module's docstring and ``docs/execution-mode-3-design.md``'s
"Adapters that cannot expose distributed state" clause. flowline_1d's state
is small per member (``2*NX+1``, e.g. 101 for the default ``N1=40``/
``N2=10`` grid) so a whole-member design is not a practical memory concern
for this application even though it cannot claim mode-3 scalability.

There is no inversion mechanism in flowline_1d (confirmed by inspection --
no ``run_model_inverse``/``inverse_step_single`` equivalent anywhere in
``applications/flowline_model``, and a grep for "inverse"/"inversion"
across that directory has zero matches), so this adapter registers no
``inverse_member`` callback, unlike ISSM's.
"""

from __future__ import annotations

from typing import Any, Mapping

import numpy as np

from src.parallelization.distributed_whole_member_adapter import (
    WholeMemberNativeAdapter,
    WholeMemberState,
)
from applications.flowline_model.examples.flowline_1d._flowline_enkf import (
    forecast_step_single,
    initialize_ensemble,
)


def _vec_inputs(icesee_kwargs: Mapping[str, Any]) -> list[str]:
    vec_inputs = icesee_kwargs.get("vec_inputs")
    if not vec_inputs:
        raise ValueError("icesee_kwargs['vec_inputs'] is required")
    return list(vec_inputs)


def _pack_ensemble(state: WholeMemberState, vec_inputs: list[str]) -> np.ndarray:
    """Rebuild the flat variable-major array ``run_model`` expects, in the
    same block order ``icesee_get_index`` produces for a single-rank member
    (see ``src/utils/tools.py::icesee_get_index``, ``var_nd``-aware branch:
    with ``comm=None``/``nranks == 1`` -- always true for this runner, since
    it never populates ``icesee_kwargs['subcomm']``/``['comm_world']`` --
    each variable's index range is a contiguous block in ``vec_inputs``
    order, exactly matching a plain ``np.concatenate`` here).
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
    """Create one persistent whole-member flowline state.

    Delegates to the existing, unmodified ``initialize_ensemble`` (which
    solves the model's initial condition and applies the same facemelt/
    nurdge-parameter setup modes 0-2 use). ``initialize_ensemble`` only
    reads ``icesee_kwargs['statevec_ens'].shape`` (never its values), so a
    lightweight ``(nd, 1)`` placeholder is enough here -- no
    ``Nens``-scaled array is ever allocated.
    """

    vec_inputs = _vec_inputs(icesee_kwargs)
    nd = int(icesee_kwargs["nd"])
    call_kwargs = dict(icesee_kwargs)
    call_kwargs["statevec_ens"] = np.zeros((nd, 1))
    updated_state = initialize_ensemble(int(member_id), **call_kwargs)
    variables = {name: np.asarray(updated_state[name]).ravel() for name in vec_inputs}
    return WholeMemberState(
        variables=variables,
        vec_inputs=vec_inputs,
        model_context={"member_id": int(member_id)},
    )


def forecast_member(
    state: WholeMemberState,
    timestep: int,
    *,
    topology: Any,
    icesee_kwargs: Mapping[str, Any],
) -> dict[str, np.ndarray] | None:
    """Advance one persistent member via the existing forecast_step_single."""

    vec_inputs = _vec_inputs(icesee_kwargs)
    ensemble = _pack_ensemble(state, vec_inputs)
    updated_state = forecast_step_single(ensemble=ensemble, **dict(icesee_kwargs))
    return {
        name: np.asarray(updated_state[name]).ravel()
        for name in vec_inputs
        if name in updated_state
    }


FLOWLINE_NATIVE_ADAPTER = WholeMemberNativeAdapter(
    initialize_member=initialize_member,
    forecast_member=forecast_member,
)
