"""Model-agnostic preflight planning for execution mode 3.

Execution mode 3 is intentionally not registered with the public runner yet.
This module resolves a model-owned native adapter and validates the distributed
contract before any expensive member is allocated.  It also applies the
analytical memory gate so applications fail early instead of exhausting a
node partway through a forecast.

The plan changes ownership and communication only.  ``global_stochastic`` and
``grouped_local_stochastic`` both retain ICESEE's existing stochastic EnKF
analysis semantics.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping

from .distributed_memory import Mode3MemoryPlan, estimate_mode3_memory
from .distributed_native_runtime import validate_native_distributed_adapter


_LOCAL_CALLBACKS = (
    "distributed_analysis_targets",
    "distributed_observation_coordinates",
)


@dataclass(frozen=True)
class NativeExecutionPlan:
    """Validated mode-3 configuration that is safe to initialize."""

    adapter: Any
    analysis_kind: str
    error_mode: str
    checkpoint_backend: str
    inversion_enabled: bool
    restart_enabled: bool
    state_row_chunk_size: int
    observation_row_chunk_size: int
    memory: Mode3MemoryPlan
    capabilities: tuple[str, ...]


def _adapter_capabilities(adapter: Any) -> tuple[str, ...]:
    optional = (
        "partition_native_observations",
        "distributed_analysis_targets",
        "distributed_observation_coordinates",
        "inverse_native_member",
        "restore_native_member_auxiliary",
    )
    return tuple(name for name in optional if callable(getattr(adapter, name, None)))


def _resolve_factory(icesee_kwargs: Mapping[str, Any]):
    direct = icesee_kwargs.get("distributed_model_adapter")
    if direct is not None:
        return direct
    factory = icesee_kwargs.get("distributed_model_adapter_factory")
    if factory is not None:
        if not callable(factory):
            raise TypeError("distributed_model_adapter_factory must be callable")
        return factory
    model_module = icesee_kwargs.get("model_module")
    module_factory = getattr(model_module, "create_distributed_model_adapter", None)
    if callable(module_factory):
        return module_factory
    raise ValueError(
        "mode 3 requires distributed_model_adapter, "
        "distributed_model_adapter_factory, or model_module."
        "create_distributed_model_adapter"
    )


def resolve_native_distributed_adapter(
    icesee_kwargs: Mapping[str, Any],
    topology: Any,
) -> tuple[Any, tuple[str, ...]]:
    """Resolve and collectively validate one model-native adapter.

    Factories receive only the topology and the shared ICESEE context.  They
    must return lightweight adapter objects; persistent model members are
    created later by :func:`initialize_native_member_pool`.
    """

    candidate = _resolve_factory(icesee_kwargs)
    if callable(candidate) and not callable(
        getattr(candidate, "initialize_native_member", None)
    ):
        candidate = candidate(topology=topology, icesee_kwargs=icesee_kwargs)

    inversion_enabled = bool(icesee_kwargs.get("inversion_flag", False))
    validate_native_distributed_adapter(
        candidate, require_inversion=inversion_enabled
    )
    capabilities = _adapter_capabilities(candidate)
    if bool(icesee_kwargs.get("local_analysis", False)):
        missing = [
            name for name in _LOCAL_CALLBACKS
            if name not in capabilities
        ]
        if missing:
            raise TypeError(
                "grouped-local mode 3 requires model localization callbacks: "
                + ", ".join(missing)
            )

    identity = (
        type(candidate).__module__,
        type(candidate).__qualname__,
        capabilities,
    )
    identities = topology.world.allgather(identity)
    if any(value != identity for value in identities):
        raise RuntimeError(
            "mode-3 ranks resolved inconsistent adapter types or capabilities"
        )
    return candidate, capabilities


def _global_state_rows(adapter: Any, icesee_kwargs: Mapping[str, Any]) -> int:
    configured = icesee_kwargs.get("mode3_global_state_rows")
    if configured is not None:
        rows = int(configured)
    else:
        callback = getattr(adapter, "distributed_global_state_size", None)
        if not callable(callback):
            raise ValueError(
                "set mode3_global_state_rows or implement "
                "distributed_global_state_size"
            )
        rows = int(callback(icesee_kwargs=icesee_kwargs))
    if rows <= 0:
        raise ValueError("mode3 global state size must be positive")
    return rows


def build_native_execution_plan(
    icesee_kwargs: Mapping[str, Any],
    topology: Any,
) -> NativeExecutionPlan:
    """Build a collective, allocation-free preflight plan for mode 3."""

    adapter, capabilities = resolve_native_distributed_adapter(
        icesee_kwargs, topology
    )
    number_of_members = int(icesee_kwargs.get("Nens", 1))
    if number_of_members < int(topology.ensemble_groups):
        raise ValueError(
            "Nens must be at least the number of mode-3 ensemble groups"
        )

    state_chunk = int(icesee_kwargs.get("mode3_state_row_chunk_size", 4096))
    observation_chunk = int(
        icesee_kwargs.get("mode3_observation_row_chunk_size", 4096)
    )
    memory = estimate_mode3_memory(
        global_rows=_global_state_rows(adapter, icesee_kwargs),
        total_members=number_of_members,
        ensemble_groups=int(topology.ensemble_groups),
        spatial_ranks=int(topology.spatial_ranks),
        state_row_chunk_size=state_chunk,
        observation_row_chunk_size=observation_chunk,
    )
    budget_gib = icesee_kwargs.get("mode3_max_rank_memory_gib")
    if budget_gib is not None:
        budget_bytes = float(budget_gib) * (1024.0 ** 3)
        if budget_bytes <= 0.0:
            raise ValueError("mode3_max_rank_memory_gib must be positive")
        if memory.estimated_minimum_peak_bytes > budget_bytes:
            required = memory.estimated_minimum_peak_bytes / (1024.0 ** 3)
            raise MemoryError(
                "mode-3 preflight exceeds the configured per-rank memory "
                f"budget: minimum estimate {required:.3f} GiB > "
                f"{float(budget_gib):.3f} GiB. This excludes model-solver "
                "workspace, so increase spatial_ranks or reduce local members."
            )

    backend = str(
        icesee_kwargs.get("mode3_checkpoint_backend", "rank_sharded_hdf5")
    ).lower()
    if backend not in {"rank_sharded_hdf5", "none"}:
        raise ValueError(
            "mode3_checkpoint_backend currently supports rank_sharded_hdf5 "
            "or none; ADIOS2 remains a planned backend"
        )

    error_mode = str(
        icesee_kwargs.get("observation_error_mode", "legacy_prior_anomalies")
    ).lower()
    if error_mode not in {
        "stochastic_r", "generated_r", "legacy_prior_anomalies"
    }:
        raise ValueError("unsupported stochastic observation-error mode")

    return NativeExecutionPlan(
        adapter=adapter,
        analysis_kind=(
            "grouped_local_stochastic"
            if bool(icesee_kwargs.get("local_analysis", False))
            else "global_stochastic"
        ),
        error_mode=error_mode,
        checkpoint_backend=backend,
        inversion_enabled=bool(icesee_kwargs.get("inversion_flag", False)),
        restart_enabled=bool(icesee_kwargs.get("restart", False)),
        state_row_chunk_size=state_chunk,
        observation_row_chunk_size=observation_chunk,
        memory=memory,
        capabilities=capabilities,
    )
