# ==============================================================================
# @des: Lightweight state-ownership metadata for ICESEE model communicators.
# @date: 2025-08-27
# ==============================================================================
"""Distinguishes how one ensemble member's flat DA state vector relates to
the ranks of its model communicator (``subcomm``).

Two ownership kinds are supported today:

``replicated``
    Every rank in the model communicator holds an identical, complete copy
    of the state -- e.g. Lorenz-96, or any adapter that does not itself
    spatially decompose across ``subcomm``. ``local_size == global_size``
    and ``local_offset == 0`` on every rank; no rank's data needs to be
    combined with another's.

``distributed``
    Each rank owns a disjoint, contiguous slice of the global state --
    e.g. a Firedrake/PETSc field partitioned across the ranks of a mesh.
    ``global_size`` is the sum of every rank's ``local_size``, and
    ``local_offset`` is that rank's start index within the global
    flattened vector.

Conflating the two is what caused a real, reproduced deadlock/shape bug in
``_mpi_generate_true_wrong_state.py`` and ``_mpi_ensemble_intialization.py``:
both summed/concatenated every rank's ``nd`` unconditionally, which is only
correct for ``distributed`` state. A replicated model like Lorenz-96 has
the same ``nd`` on every rank of a model communicator, so summing it
inflated the state size and stacking identical per-rank copies duplicated
data instead of assembling a partition.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np

_VALID_DISTRIBUTIONS = ("replicated", "distributed")


@dataclass(frozen=True)
class StateOwnership:
    """Describes one ensemble member's state layout across its ``subcomm``.

    This describes the *flattened DA representation* only. A model whose
    native field is logically multidimensional (e.g. ``(nz, ny, nx)``) can
    still report the flattened ``local_size``/``global_size`` here;
    ``model_shape`` is an optional hook for an adapter to retain the
    logical shape alongside it, without requiring shared ICESEE runtime
    code to understand or reconstruct that shape.
    """

    distribution: str
    global_size: int
    local_size: int
    local_offset: int = 0
    model_shape: Optional[Tuple[int, ...]] = None

    def __post_init__(self):
        if self.distribution not in _VALID_DISTRIBUTIONS:
            raise ValueError(
                "StateOwnership.distribution must be one of "
                f"{_VALID_DISTRIBUTIONS}, got {self.distribution!r}."
            )


def resolve_state_ownership(icesee_kwargs, subcomm, local_size=None) -> StateOwnership:
    """Determine one ensemble member's state ownership from its ``subcomm``.

    ``local_size`` defaults to ``icesee_kwargs["nd"]`` -- the value every
    existing model adapter already reports for its own rank. What that
    value *means* is exactly the ownership question this resolves: for a
    replicated model it is the (identical) global size; for a distributed
    model it is this rank's own partition size.

    The distribution kind comes from ``icesee_kwargs["state_distribution"]``
    (``"replicated"`` or ``"distributed"``), defaulting to ``"distributed"``.
    That default is deliberate, not a guess: every ICESEE model adapter
    shipped today (Icepack, ISSM) already reports a genuine per-rank local
    partition size in ``nd`` -- summing it across a model communicator is
    exactly what the pre-existing (pre-ownership-aware) code always did,
    unconditionally, and it is arithmetically correct for those adapters.
    Preserving that as the default means this function changes nothing for
    any model that was already correct. Lorenz-96 is the one adapter whose
    ``nd`` is a replicated global size, not a partition; it opts in
    explicitly via ``state_distribution: replicated`` in its own config
    (see ``applications/lorenz_model/examples/lorenz96/params.yaml``)
    rather than changing the shared default, because Lorenz-96 never
    produced correct results at ``ranks_per_model > 1`` under the old
    unconditional-sum behavior either (it reproducibly deadlocks/raises a
    shape error) -- so there is no existing working behavior to preserve
    for it, only a genuine fix to opt into.
    """
    if local_size is None:
        local_size = int(icesee_kwargs.get("nd", icesee_kwargs["nd"]))
    else:
        local_size = int(local_size)

    distribution = str(
        icesee_kwargs.get("state_distribution", "distributed")
    ).strip().lower()
    if distribution not in _VALID_DISTRIBUTIONS:
        raise ValueError(
            "state_distribution must be one of "
            f"{_VALID_DISTRIBUTIONS}, got {distribution!r}."
        )

    sub_size = subcomm.Get_size() if subcomm is not None else 1

    if distribution == "replicated" or sub_size <= 1:
        return StateOwnership(
            distribution="replicated",
            global_size=local_size,
            local_size=local_size,
            local_offset=0,
        )

    dim_list = subcomm.allgather(local_size)
    sub_rank = subcomm.Get_rank()
    global_size = int(sum(dim_list))
    local_offset = int(sum(dim_list[:sub_rank]))
    return StateOwnership(
        distribution="distributed",
        global_size=global_size,
        local_size=local_size,
        local_offset=local_offset,
    )


def combine_member_state(ownership: StateOwnership, subcomm, data, root=0):
    """Assemble one ensemble member's global state from per-rank pieces.

    ``data`` is this rank's contribution: either a dict of ``{name: array}``
    (matching the ``generate_true_state``/``initialize_ensemble`` adapter
    convention) or a plain array. Returns the assembled global data on
    ``root`` (same container type as ``data``); the return value on other
    ranks is ``None`` and must not be used, matching the pre-existing
    ``sub_rank == 0``-only convention this replaces.

    For ``replicated`` state this does **no communication at all**: every
    rank's ``data`` is already, by definition, the complete global state,
    so ``root``'s own copy is returned unchanged -- strictly cheaper than
    the ``distributed`` path, not just simpler.

    For ``distributed`` state this gathers each rank's disjoint slice onto
    ``root`` and concatenates them along axis 0 into the assembled whole.
    ``np.concatenate(..., axis=0)`` is used rather than ``np.vstack``
    deliberately: for the 2D ``(local_rows, nt+1)`` pieces the true/nurged
    state path produces, the two are equivalent (stack as additional rows);
    for the 1D per-variable vectors the ensemble-initialization path
    produces, ``vstack`` would wrongly promote each piece to a row of a new
    2D array, whereas ``concatenate`` correctly joins them end-to-end into
    one longer 1D vector -- matching the pre-existing ``hstack`` behavior
    at that call site.
    """
    if ownership.distribution == "replicated":
        return data

    sub_rank = subcomm.Get_rank()
    if isinstance(data, dict):
        gathered = {
            key: subcomm.gather(value, root=root) for key, value in data.items()
        }
        if sub_rank != root:
            return None
        return {
            key: np.concatenate(pieces, axis=0) for key, pieces in gathered.items()
        }

    gathered = subcomm.gather(data, root=root)
    if sub_rank != root:
        return None
    return np.concatenate(gathered, axis=0)


def ensure_state_array(value):
    """Normalize one model variable's contribution to at least 1-D.

    Model adapters return one value per state variable (the
    ``generate_true_state``/``generate_nurged_state``/``initialize_ensemble``
    ``{name: value}`` convention). When a variable's local block has exactly
    one element (``hdim == 1``), some adapters return a bare NumPy scalar
    (e.g. Lorenz-96's ``initialize_ensemble`` returning ``u0b[0]``, a
    ``numpy.float64``) rather than a one-element array -- both represent the
    same one-element state *component*, just different NumPy
    representations of it. This is unrelated to state *ownership*: it is
    not a distribution concept, and a distributed model's local *partition*
    is already guaranteed to be an array by construction (sliced from a
    larger state), so this only ever matters for the scalar/0-D edge case.

    ``np.atleast_1d`` is exactly this normalization: it leaves any array
    that is already >=1-D unchanged -- including multidimensional ones, so
    this is safe to apply unconditionally before shape-derived logic
    (``.shape[0]``), concatenation, or in-place arithmetic -- and promotes a
    bare scalar or 0-D array to shape ``(1,)``.
    """
    return np.atleast_1d(np.asarray(value))
