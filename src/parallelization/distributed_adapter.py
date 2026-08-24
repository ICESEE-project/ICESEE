"""Contracts for model-native distributed state in execution mode 3.

The interfaces in this module are deliberately non-selectable: execution mode
3 is not registered with the runtime yet.  They let model adapters describe a
rank-local state slab and let ICESEE validate that a spatial communicator owns
the global state exactly once before any distributed forecast or analysis is
attempted.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping, Protocol, Sequence, TypeAlias, runtime_checkable

import numpy as np


@dataclass(frozen=True)
class DistributedStateLayout:
    """Contiguous owned state slab and optional read-only ghost indices.

    The first mode-3 contract uses contiguous ownership in a stable global
    numbering.  This matches distributed vectors used by PETSc and avoids
    communicating complete index arrays during validation or checkpointing.
    Mesh-local ordering may still be arbitrary inside a model adapter; the
    adapter is responsible for mapping it to this global numbering.
    """

    global_size: int
    owned_start: int
    owned_stop: int
    ghost_indices: np.ndarray = field(
        default_factory=lambda: np.empty(0, dtype=np.int64)
    )
    layout_id: str = "state-v1"

    def __post_init__(self) -> None:
        global_size = int(self.global_size)
        owned_start = int(self.owned_start)
        owned_stop = int(self.owned_stop)
        if global_size < 0:
            raise ValueError("global_size must be nonnegative")
        if not 0 <= owned_start <= owned_stop <= global_size:
            raise ValueError(
                "owned slab must satisfy 0 <= owned_start <= owned_stop "
                "<= global_size"
            )
        if not str(self.layout_id):
            raise ValueError("layout_id must be nonempty")

        ghosts = np.asarray(self.ghost_indices, dtype=np.int64)
        if ghosts.ndim != 1:
            raise ValueError("ghost_indices must be one-dimensional")
        if ghosts.size:
            ghosts = np.unique(ghosts)
            if ghosts[0] < 0 or ghosts[-1] >= global_size:
                raise ValueError("ghost_indices must lie inside the global state")
            if np.any((ghosts >= owned_start) & (ghosts < owned_stop)):
                raise ValueError("ghost_indices cannot overlap the owned slab")

        object.__setattr__(self, "global_size", global_size)
        object.__setattr__(self, "owned_start", owned_start)
        object.__setattr__(self, "owned_stop", owned_stop)
        object.__setattr__(self, "ghost_indices", ghosts)
        object.__setattr__(self, "layout_id", str(self.layout_id))

    @property
    def owned_size(self) -> int:
        """Number of globally owned state entries on this spatial rank."""

        return self.owned_stop - self.owned_start

    @property
    def owned_slice(self) -> slice:
        """Slice selecting the owned entries in the stable global numbering."""

        return slice(self.owned_start, self.owned_stop)

    def checkpoint_metadata(self) -> dict[str, int | str]:
        """Return rank-placement metadata required for distributed restart."""

        return {
            "layout_id": self.layout_id,
            "global_size": self.global_size,
            "owned_start": self.owned_start,
            "owned_stop": self.owned_stop,
        }

    def local_to_global_rows(self, local_rows: np.ndarray) -> np.ndarray:
        """Map owned-array row numbers to stable global state rows."""

        rows = np.asarray(local_rows, dtype=np.int64).ravel()
        if rows.size and (rows.min() < 0 or rows.max() >= self.owned_size):
            raise ValueError("local rows leave the owned state slab")
        return self.owned_start + rows


@dataclass(frozen=True)
class DistributedStateBlockLayout:
    """Owned interval for one named variable in a variable-major state.

    ``owned_start`` and ``owned_stop`` use block-local numbering.  The
    ``global_offset`` retains ICESEE's existing variable-major packed-state
    ordering without requiring a rank to own the gaps between its local
    pieces of consecutive model vectors.
    """

    name: str
    global_offset: int
    global_size: int
    owned_start: int
    owned_stop: int
    ghost_indices: np.ndarray = field(
        default_factory=lambda: np.empty(0, dtype=np.int64)
    )

    def __post_init__(self) -> None:
        name = str(self.name)
        offset = int(self.global_offset)
        size = int(self.global_size)
        start = int(self.owned_start)
        stop = int(self.owned_stop)
        if not name:
            raise ValueError("state block name must be nonempty")
        if offset < 0 or size < 0:
            raise ValueError("state block offset and size must be nonnegative")
        if not 0 <= start <= stop <= size:
            raise ValueError("state block ownership lies outside the block")
        ghosts = np.asarray(self.ghost_indices, dtype=np.int64).ravel()
        if ghosts.size:
            ghosts = np.unique(ghosts)
            if ghosts[0] < 0 or ghosts[-1] >= size:
                raise ValueError("state block ghost indices leave the block")
            if np.any((ghosts >= start) & (ghosts < stop)):
                raise ValueError("state block ghosts overlap owned entries")
        object.__setattr__(self, "name", name)
        object.__setattr__(self, "global_offset", offset)
        object.__setattr__(self, "global_size", size)
        object.__setattr__(self, "owned_start", start)
        object.__setattr__(self, "owned_stop", stop)
        object.__setattr__(self, "ghost_indices", ghosts)

    @property
    def owned_size(self) -> int:
        return self.owned_stop - self.owned_start

    @property
    def global_owned_interval(self) -> tuple[int, int]:
        return (
            self.global_offset + self.owned_start,
            self.global_offset + self.owned_stop,
        )


@dataclass(frozen=True)
class DistributedBlockStateLayout:
    """Segmented ownership for model-native multi-variable state vectors.

    Local arrays concatenate only this rank's owned portion of each block in
    block order.  They therefore scale with local degrees of freedom even
    though their global rows remain compatible with the existing ICESEE
    variable-major state ordering.
    """

    blocks: tuple[DistributedStateBlockLayout, ...]
    layout_id: str = "block-state-v1"

    def __post_init__(self) -> None:
        blocks = tuple(self.blocks)
        if not blocks:
            raise ValueError("block state layout requires at least one block")
        if any(not isinstance(block, DistributedStateBlockLayout) for block in blocks):
            raise TypeError("blocks must use DistributedStateBlockLayout")
        names = [block.name for block in blocks]
        if len(set(names)) != len(names):
            raise ValueError("state block names must be unique")
        cursor = 0
        for block in blocks:
            if block.global_offset != cursor:
                raise ValueError("state blocks must form a gap-free global ordering")
            cursor += block.global_size
        layout_id = str(self.layout_id)
        if not layout_id:
            raise ValueError("layout_id must be nonempty")
        object.__setattr__(self, "blocks", blocks)
        object.__setattr__(self, "layout_id", layout_id)

    @property
    def global_size(self) -> int:
        return sum(block.global_size for block in self.blocks)

    @property
    def owned_size(self) -> int:
        return sum(block.owned_size for block in self.blocks)

    @property
    def block_names(self) -> tuple[str, ...]:
        return tuple(block.name for block in self.blocks)

    def block(self, name: str) -> DistributedStateBlockLayout:
        for block in self.blocks:
            if block.name == str(name):
                return block
        raise KeyError(f"unknown distributed state block: {name}")

    def local_slice(self, name: str) -> slice:
        cursor = 0
        for block in self.blocks:
            stop = cursor + block.owned_size
            if block.name == str(name):
                return slice(cursor, stop)
            cursor = stop
        raise KeyError(f"unknown distributed state block: {name}")

    def iter_owned_global_intervals(self):
        """Yield block name, local slice, and global interval without indices."""

        cursor = 0
        for block in self.blocks:
            local = slice(cursor, cursor + block.owned_size)
            yield block.name, local, block.global_owned_interval
            cursor = local.stop

    def checkpoint_metadata(self) -> dict[str, Any]:
        return {
            "layout_id": self.layout_id,
            "global_size": self.global_size,
            "blocks": [
                {
                    "name": block.name,
                    "global_offset": block.global_offset,
                    "global_size": block.global_size,
                    "owned_start": block.owned_start,
                    "owned_stop": block.owned_stop,
                }
                for block in self.blocks
            ],
        }

    def local_to_global_rows(self, local_rows: np.ndarray) -> np.ndarray:
        """Map rows in the compact segmented array to variable-major rows."""

        rows = np.asarray(local_rows, dtype=np.int64).ravel()
        if rows.size and (rows.min() < 0 or rows.max() >= self.owned_size):
            raise ValueError("local rows leave the owned block state")
        result = np.empty(rows.size, dtype=np.int64)
        for block in self.blocks:
            local = self.local_slice(block.name)
            selected = (rows >= local.start) & (rows < local.stop)
            result[selected] = (
                block.global_offset
                + block.owned_start
                + rows[selected]
                - local.start
            )
        return result


DistributedLayout: TypeAlias = DistributedStateLayout | DistributedBlockStateLayout


@dataclass(frozen=True)
class DistributedAnalysisTargets:
    """Locally owned state rows participating in one localized update.

    Coordinates may have any spatial dimension.  ``local_rows`` index this
    rank's owned state slab while ``global_rows`` retain the stable model-wide
    numbering used for deterministic grouping, diagnostics, and restart.
    ``observation_kinds`` lets a model associate (for example) thickness state
    rows with thickness and velocity observations without embedding model
    knowledge in the parallel runtime.
    """

    name: str
    local_rows: np.ndarray
    global_rows: np.ndarray
    coordinates: np.ndarray
    observation_kinds: tuple[str, ...]

    def __post_init__(self) -> None:
        local_rows = np.asarray(self.local_rows, dtype=np.int64).ravel()
        global_rows = np.asarray(self.global_rows, dtype=np.int64).ravel()
        coordinates = np.asarray(self.coordinates, dtype=np.float64)
        if coordinates.ndim == 1:
            coordinates = coordinates[:, None]
        if local_rows.shape != global_rows.shape:
            raise ValueError("local_rows and global_rows must have equal lengths")
        if coordinates.ndim != 2 or coordinates.shape[0] != local_rows.size:
            raise ValueError("target coordinates must have one row per state row")
        if local_rows.size and (
            np.unique(local_rows).size != local_rows.size
            or np.unique(global_rows).size != global_rows.size
        ):
            raise ValueError("analysis target rows must be unique")
        kinds = tuple(str(value) for value in self.observation_kinds)
        if not str(self.name) or not kinds or any(not value for value in kinds):
            raise ValueError("target name and observation kinds must be nonempty")
        if not np.all(np.isfinite(coordinates)):
            raise ValueError("target coordinates contain non-finite values")
        object.__setattr__(self, "name", str(self.name))
        object.__setattr__(self, "local_rows", local_rows)
        object.__setattr__(self, "global_rows", global_rows)
        object.__setattr__(self, "coordinates", coordinates)
        object.__setattr__(self, "observation_kinds", kinds)


@dataclass(frozen=True)
class DistributedObservationCoordinates:
    """Canonical locally owned observation IDs and physical coordinates."""

    kind: str
    observation_ids: np.ndarray
    coordinates: np.ndarray

    def __post_init__(self) -> None:
        observation_ids = np.asarray(self.observation_ids, dtype=np.int64).ravel()
        coordinates = np.asarray(self.coordinates, dtype=np.float64)
        if coordinates.ndim == 1:
            coordinates = coordinates[:, None]
        if coordinates.ndim != 2 or coordinates.shape[0] != observation_ids.size:
            raise ValueError(
                "observation coordinates must have one row per observation ID"
            )
        if observation_ids.size and np.unique(observation_ids).size != observation_ids.size:
            raise ValueError("locally owned observation IDs must be unique")
        if not str(self.kind) or not np.all(np.isfinite(coordinates)):
            raise ValueError("observation kind and coordinates must be valid")
        object.__setattr__(self, "kind", str(self.kind))
        object.__setattr__(self, "observation_ids", observation_ids)
        object.__setattr__(self, "coordinates", coordinates)


@dataclass(frozen=True)
class DistributedObservationCollection:
    """Joint local observation metadata for one stochastic target update.

    Mode 1 constructs one analysis from every active observation kind that may
    influence a target block. Mode 3 must preserve that ordering instead of
    applying one filter update per kind. This collection concatenates only the
    current spatial rank's metadata and never forms a global observation table.
    IDs must therefore be canonical and unique across all kinds.
    """

    observations: tuple[DistributedObservationCoordinates, ...]

    def __post_init__(self) -> None:
        observations = tuple(self.observations)
        if not observations:
            raise ValueError("observation collection must be nonempty")
        if any(
            not isinstance(item, DistributedObservationCoordinates)
            for item in observations
        ):
            raise TypeError(
                "observation collection entries must use "
                "DistributedObservationCoordinates"
            )
        kinds = tuple(item.kind for item in observations)
        if len(set(kinds)) != len(kinds):
            raise ValueError("observation kinds must be unique in a collection")
        dimensions = {int(item.coordinates.shape[1]) for item in observations}
        if len(dimensions) != 1:
            raise ValueError("observation coordinate dimensions differ")
        ids = np.concatenate([item.observation_ids for item in observations])
        if ids.size and np.unique(ids).size != ids.size:
            raise ValueError(
                "canonical observation IDs must be unique across kinds"
            )
        object.__setattr__(self, "observations", observations)

    @property
    def kinds(self) -> tuple[str, ...]:
        return tuple(item.kind for item in self.observations)

    @property
    def observation_ids(self) -> np.ndarray:
        return np.concatenate(
            [item.observation_ids for item in self.observations]
        )

    @property
    def coordinates(self) -> np.ndarray:
        return np.concatenate(
            [item.coordinates for item in self.observations], axis=0
        )


def validate_distributed_analysis_coordinates(
    targets: Sequence[DistributedAnalysisTargets],
    observations: Sequence[DistributedObservationCoordinates],
    layout: DistributedLayout,
) -> None:
    """Validate model-provided localization metadata without global fields."""

    dimensions: set[int] = set()
    target_names: set[str] = set()
    observation_kinds: set[str] = set()
    for target in targets:
        if not isinstance(target, DistributedAnalysisTargets):
            raise TypeError("analysis targets must use DistributedAnalysisTargets")
        if target.name in target_names:
            raise ValueError(f"duplicate analysis target name: {target.name}")
        target_names.add(target.name)
        dimensions.add(int(target.coordinates.shape[1]))
        if target.local_rows.size and (
            target.local_rows.min() < 0
            or target.local_rows.max() >= layout.owned_size
        ):
            raise ValueError(f"analysis target {target.name!r} leaves owned state slab")
        expected = layout.local_to_global_rows(target.local_rows)
        if not np.array_equal(target.global_rows, expected):
            raise ValueError(
                f"analysis target {target.name!r} has inconsistent global rows"
            )
    all_observation_ids: list[np.ndarray] = []
    for observation in observations:
        if not isinstance(observation, DistributedObservationCoordinates):
            raise TypeError(
                "observation metadata must use DistributedObservationCoordinates"
            )
        if observation.kind in observation_kinds:
            raise ValueError(f"duplicate observation kind: {observation.kind}")
        observation_kinds.add(observation.kind)
        all_observation_ids.append(observation.observation_ids)
        dimensions.add(int(observation.coordinates.shape[1]))
    if all_observation_ids:
        combined_ids = np.concatenate(all_observation_ids)
        if combined_ids.size and np.unique(combined_ids).size != combined_ids.size:
            raise ValueError(
                "canonical observation IDs must be unique across kinds"
            )
    if len(dimensions) > 1:
        raise ValueError("target and observation coordinate dimensions differ")
    missing = sorted({kind for target in targets for kind in target.observation_kinds}
                     - observation_kinds)
    if missing:
        raise ValueError(f"analysis targets reference missing observations: {missing}")


def validate_spatial_partition(
    layout: DistributedStateLayout,
    spatial_comm: Any,
) -> tuple[tuple[int, int], ...]:
    """Verify exact, gap-free ownership across one member's spatial ranks.

    Only four scalar values per rank are gathered.  The check therefore remains
    cheap for very large state vectors and never reconstructs a global member.
    The returned intervals are sorted by global start and are useful when
    building collective checkpoint hyperslabs.
    """

    local_descriptor = (
        layout.owned_start,
        layout.owned_stop,
        layout.global_size,
        layout.layout_id,
    )
    descriptors = spatial_comm.allgather(local_descriptor)
    if not descriptors:
        raise ValueError("spatial communicator returned no ownership descriptors")

    global_sizes = {int(item[2]) for item in descriptors}
    layout_ids = {str(item[3]) for item in descriptors}
    if global_sizes != {layout.global_size}:
        raise ValueError("spatial ranks disagree on global state size")
    if layout_ids != {layout.layout_id}:
        raise ValueError("spatial ranks disagree on layout_id")

    intervals = sorted((int(item[0]), int(item[1])) for item in descriptors)
    cursor = 0
    for start, stop in intervals:
        if not 0 <= start <= stop <= layout.global_size:
            raise ValueError("a spatial rank reported an invalid owned slab")
        if start < cursor:
            raise ValueError("distributed state ownership overlaps")
        if start > cursor:
            raise ValueError("distributed state ownership contains a gap")
        cursor = stop
    if cursor != layout.global_size:
        raise ValueError("distributed state ownership does not cover global state")

    return tuple(intervals)


def validate_block_spatial_partition(
    layout: DistributedBlockStateLayout,
    spatial_comm: Any,
) -> dict[str, tuple[tuple[int, int], ...]]:
    """Validate exact ownership independently for every model state block.

    Communication is proportional to the number of variables times the number
    of spatial ranks, never to the number of state degrees of freedom.
    """

    descriptor = (
        layout.layout_id,
        tuple(
            (
                block.name,
                block.global_offset,
                block.global_size,
                block.owned_start,
                block.owned_stop,
            )
            for block in layout.blocks
        ),
    )
    descriptors = spatial_comm.allgather(descriptor)
    if not descriptors:
        raise ValueError("spatial communicator returned no block ownership")
    if {str(item[0]) for item in descriptors} != {layout.layout_id}:
        raise ValueError("spatial ranks disagree on block layout_id")

    reference = tuple(
        (block.name, block.global_offset, block.global_size)
        for block in layout.blocks
    )
    for _, remote_blocks in descriptors:
        remote_definition = tuple(
            (str(name), int(offset), int(size))
            for name, offset, size, _, _ in remote_blocks
        )
        if remote_definition != reference:
            raise ValueError("spatial ranks disagree on state block definitions")

    result: dict[str, tuple[tuple[int, int], ...]] = {}
    for block_index, block in enumerate(layout.blocks):
        intervals = sorted(
            (
                int(remote_blocks[block_index][3]),
                int(remote_blocks[block_index][4]),
            )
            for _, remote_blocks in descriptors
        )
        cursor = 0
        for start, stop in intervals:
            if not 0 <= start <= stop <= block.global_size:
                raise ValueError(
                    f"state block {block.name!r} has invalid ownership"
                )
            if start < cursor:
                raise ValueError(f"state block {block.name!r} ownership overlaps")
            if start > cursor:
                raise ValueError(f"state block {block.name!r} ownership has a gap")
            cursor = stop
        if cursor != block.global_size:
            raise ValueError(
                f"state block {block.name!r} ownership does not cover the block"
            )
        result[block.name] = tuple(intervals)
    return result


def contiguous_block_state_layout(
    block_sizes: Mapping[str, int] | Sequence[tuple[str, int]],
    topology: Any,
    *,
    layout_id: str = "block-state-v1",
) -> DistributedBlockStateLayout:
    """Create balanced model-native ownership for each named state variable."""

    entries = (
        tuple(block_sizes.items())
        if isinstance(block_sizes, Mapping)
        else tuple(block_sizes)
    )
    if not entries:
        raise ValueError("block_sizes must be nonempty")
    spatial_ranks = int(topology.spatial_ranks)
    spatial_rank = int(topology.spatial_rank)
    if spatial_ranks <= 0 or not 0 <= spatial_rank < spatial_ranks:
        raise ValueError("topology contains invalid spatial rank coordinates")

    blocks = []
    global_offset = 0
    for name, raw_size in entries:
        size = int(raw_size)
        if size < 0:
            raise ValueError("state block sizes must be nonnegative")
        quotient, remainder = divmod(size, spatial_ranks)
        start = spatial_rank * quotient + min(spatial_rank, remainder)
        stop = start + quotient + (spatial_rank < remainder)
        blocks.append(
            DistributedStateBlockLayout(
                name=str(name),
                global_offset=global_offset,
                global_size=size,
                owned_start=start,
                owned_stop=stop,
            )
        )
        global_offset += size
    return DistributedBlockStateLayout(tuple(blocks), layout_id=layout_id)


def contiguous_state_layout(
    global_size: int,
    topology: Any,
    *,
    ghost_indices: np.ndarray | None = None,
    layout_id: str = "state-v1",
) -> DistributedStateLayout:
    """Create a balanced contiguous slab for one spatial communicator rank."""

    global_size = int(global_size)
    spatial_ranks = int(topology.spatial_ranks)
    spatial_rank = int(topology.spatial_rank)
    if spatial_ranks <= 0 or not 0 <= spatial_rank < spatial_ranks:
        raise ValueError("topology contains invalid spatial rank coordinates")

    quotient, remainder = divmod(global_size, spatial_ranks)
    owned_start = spatial_rank * quotient + min(spatial_rank, remainder)
    owned_stop = owned_start + quotient + (spatial_rank < remainder)
    return DistributedStateLayout(
        global_size=global_size,
        owned_start=owned_start,
        owned_stop=owned_stop,
        ghost_indices=(
            np.empty(0, dtype=np.int64)
            if ghost_indices is None
            else ghost_indices
        ),
        layout_id=layout_id,
    )


@runtime_checkable
class DistributedModelAdapter(Protocol):
    """Required mode-3 callbacks for a spatially distributed model adapter."""

    def distributed_state_layout(
        self,
        *,
        topology: Any,
        icesee_kwargs: Mapping[str, Any],
    ) -> DistributedLayout:
        """Describe the stable global state and this rank's owned slab."""

    def initialize_local_member(
        self,
        member_id: int,
        *,
        layout: DistributedLayout,
        topology: Any,
        icesee_kwargs: Mapping[str, Any],
    ) -> np.ndarray:
        """Return this rank's initial owned state for one ensemble member."""

    def forecast_local_member(
        self,
        local_state: np.ndarray,
        member_id: int,
        timestep: int,
        *,
        layout: DistributedLayout,
        topology: Any,
        icesee_kwargs: Mapping[str, Any],
    ) -> np.ndarray:
        """Advance one local member slab, including model-native halo exchange."""

    def observe_local_member(
        self,
        local_state: np.ndarray,
        observation_rows: np.ndarray,
        *,
        layout: DistributedLayout,
        topology: Any,
        icesee_kwargs: Mapping[str, Any],
    ) -> np.ndarray:
        """Evaluate locally owned observation rows without gathering the state."""

    def finalize_local_analysis(
        self,
        local_forecast: np.ndarray,
        local_analysis: np.ndarray,
        member_id: int,
        timestep: int,
        *,
        layout: DistributedLayout,
        topology: Any,
        icesee_kwargs: Mapping[str, Any],
    ) -> np.ndarray:
        """Apply local constraints and geometry consistency after analysis."""


@runtime_checkable
class DistributedInversionAdapter(Protocol):
    """Additional callback required by a mode-3 hybrid inversion workflow."""

    def inverse_local_member(
        self,
        local_state: np.ndarray,
        member_id: int,
        timestep: int,
        *,
        layout: DistributedLayout,
        topology: Any,
        icesee_kwargs: Mapping[str, Any],
    ) -> np.ndarray:
        """Run inversion without assembling the complete member on rank zero."""


@runtime_checkable
class DistributedLocalizationAdapter(Protocol):
    """Optional mode-3 callbacks for model-native localized analysis metadata."""

    def distributed_analysis_targets(
        self,
        *,
        layout: DistributedLayout,
        topology: Any,
        icesee_kwargs: Mapping[str, Any],
    ) -> Sequence[DistributedAnalysisTargets]:
        """Return locally owned state-row groups and physical coordinates."""

    def distributed_observation_coordinates(
        self,
        observation_rows: np.ndarray,
        *,
        topology: Any,
        icesee_kwargs: Mapping[str, Any],
    ) -> Sequence[DistributedObservationCoordinates]:
        """Return canonical locally owned observation coordinates by kind."""


_REQUIRED_DISTRIBUTED_CALLBACKS = (
    "distributed_state_layout",
    "initialize_local_member",
    "forecast_local_member",
    "observe_local_member",
    "finalize_local_analysis",
)


def validate_distributed_adapter(adapter: Any, *, require_inversion: bool = False) -> None:
    """Fail early when a model does not implement the mode-3 callback contract."""

    missing = [
        name
        for name in _REQUIRED_DISTRIBUTED_CALLBACKS
        if not callable(getattr(adapter, name, None))
    ]
    if require_inversion and not callable(getattr(adapter, "inverse_local_member", None)):
        missing.append("inverse_local_member")
    if missing:
        raise TypeError(
            "distributed model adapter is missing required callbacks: "
            + ", ".join(missing)
        )
