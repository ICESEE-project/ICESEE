"""Rank-sharded, bounded-memory checkpoints for execution mode 3.

The first backend deliberately uses ordinary serial HDF5 files: every MPI rank
writes its own shard and no rank gathers a complete member.  A manifest is
published only after every shard is durable.  Restart is independent of rank
placement because shards are addressed by stable member IDs and global state
intervals rather than by the rank that created them.

Commit protocol (unchanged since the first v1/v2 checkpoints, so every
checkpoint ever written by this module remains readable):

1. every rank writes ``.shard_<rank>.h5.tmp`` inside
   ``.checkpoint_<step>.<run_id>.staging/``, fsyncs it and renames it to
   ``shard_<rank>.h5``;
2. world rank 0 verifies member/block coverage, writes ``manifest.json``
   (``"complete": true``) and atomically renames the staging directory to
   ``checkpoint_<step>``.

A checkpoint is *committed* only once step 2 has completed.  Dot-prefixed
staging directories and ``.tmp`` files are never committed state.
:func:`validate_distributed_checkpoint` re-checks a committed directory
against its manifest before it is used for a restart.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
import os
from pathlib import Path
import re
import shutil
from typing import Any, Iterable, Mapping

import h5py
import numpy as np

from .distributed_adapter import (
    DistributedBlockStateLayout,
    DistributedStateLayout,
    validate_block_spatial_partition,
    validate_spatial_partition,
)
from .distributed_runtime import LocalMemberEnsemble, members_for_ensemble_slot


_FORMAT = "icesee-distributed-checkpoint-v1"
_BLOCK_FORMAT = "icesee-distributed-block-checkpoint-v2"
_SUPPORTED_FORMATS = {_FORMAT, _BLOCK_FORMAT}
_COMMITTED_NAME = re.compile(r"checkpoint_(\d{8})")
_SHARD_NAME = re.compile(r"shard_\d{8}\.h5")
# Transient artifacts: an uncommitted write, or a committed checkpoint that
# retention has already un-published (renamed) and is deleting.
_ABANDONED_SUFFIXES = (".staging", ".deleting")


def _checkpoint_name(timestep: int) -> str:
    return f"checkpoint_{int(timestep):08d}"


def _fsync_path(path: Path) -> None:
    """Flush a file or directory entry to stable storage where supported."""

    try:
        descriptor = os.open(path, os.O_RDONLY)
    except OSError:
        return
    try:
        os.fsync(descriptor)
    except OSError:
        # Some filesystems do not support fsync on directories.
        pass
    finally:
        os.close(descriptor)


def _json_safe_metadata(metadata: Mapping[str, Any] | None) -> dict[str, Any]:
    if metadata is None:
        return {}
    result = dict(metadata)
    # Fail before any collective write rather than producing a checkpoint that
    # cannot later be interpreted consistently.
    json.dumps(result)
    return result


@dataclass(frozen=True)
class DistributedCheckpoint:
    """Metadata returned after a checkpoint is committed or loaded."""

    path: Path
    timestep: int
    run_id: str
    metadata: Mapping[str, Any]


def latest_complete_distributed_checkpoint(
    checkpoint_root: str | os.PathLike[str],
    *,
    expected_run_id: str | None = None,
) -> Path | None:
    """Return the newest manifest-committed checkpoint, ignoring debris.

    Staging directories and checkpoint-looking directories without a valid,
    complete manifest (or with missing shards) are intentionally ignored.  This lets a restarted job
    recover from the last durable cycle after a writer exception or scheduler
    termination without treating a partial write as model state.
    """

    return _discover_local(
        Path(checkpoint_root), expected_run_id, None, deep=False
    ).path


def save_distributed_checkpoint(
    checkpoint_root: str | os.PathLike[str],
    timestep: int,
    local_ensemble: LocalMemberEnsemble,
    topology: Any,
    *,
    run_id: str,
    metadata: Mapping[str, Any] | None = None,
    layout_fingerprint: str | None = None,
) -> DistributedCheckpoint:
    """Atomically save distributed state without reconstructing any member.

    ``metadata`` is the model-agnostic sidecar for random-stream keys,
    observation-cycle state, inflation/localization settings, and inversion
    schedule state.  Model-specific restart data should be represented by
    stable identifiers or written by a model adapter into its own shard.

    ``layout_fingerprint`` is an optional, rank-local identity of the model's
    degree-of-freedom numbering (see
    ``native_layout_fingerprint`` in ``distributed_native_runtime``).  It is
    recorded per shard so a restart can refuse to load state into a
    differently numbered mesh.
    """

    root = Path(checkpoint_root)
    timestep = int(timestep)
    if timestep < 0:
        raise ValueError("checkpoint timestep must be nonnegative")
    run_id = str(run_id)
    if not run_id:
        raise ValueError("run_id must be nonempty")
    sidecar = _json_safe_metadata(metadata)
    layout = local_ensemble.layout
    if isinstance(layout, DistributedBlockStateLayout):
        return _save_distributed_block_checkpoint(
            root,
            timestep,
            local_ensemble,
            topology,
            run_id=run_id,
            metadata=sidecar,
            layout_fingerprint=layout_fingerprint,
        )
    if not isinstance(layout, DistributedStateLayout):
        raise TypeError("unsupported distributed checkpoint layout")
    intervals = validate_spatial_partition(layout, topology.spatial_comm)

    final_dir = root / _checkpoint_name(timestep)
    staging_dir = root / f".{_checkpoint_name(timestep)}.{run_id}.staging"
    setup_error = None
    if int(topology.world_rank) == 0:
        try:
            root.mkdir(parents=True, exist_ok=True)
            if final_dir.exists():
                raise FileExistsError(f"checkpoint already exists: {final_dir}")
            if staging_dir.exists():
                shutil.rmtree(staging_dir)
            staging_dir.mkdir(parents=False)
        except Exception as error:  # broadcast before another rank can block
            setup_error = f"{type(error).__name__}: {error}"
    setup_error = topology.world.bcast(setup_error, root=0)
    if setup_error is not None:
        raise RuntimeError(f"cannot prepare distributed checkpoint: {setup_error}")
    topology.world.Barrier()

    shard_name = f"shard_{int(topology.world_rank):08d}.h5"
    shard_path = staging_dir / shard_name
    shard_tmp = staging_dir / f".{shard_name}.tmp"
    descriptor = None
    local_write_error = None
    try:
        local_dtypes = {
            np.asarray(local_ensemble.members[member_id]).dtype.str
            for member_id in local_ensemble.member_ids
        }
        if len(local_dtypes) > 1:
            raise TypeError("checkpoint members on one rank must share a dtype")
        with h5py.File(shard_tmp, "w") as handle:
            handle.attrs["format"] = _FORMAT
            handle.attrs["run_id"] = run_id
            handle.attrs["timestep"] = timestep
            handle.attrs["layout_id"] = layout.layout_id
            handle.attrs["global_size"] = layout.global_size
            handle.attrs["owned_start"] = layout.owned_start
            handle.attrs["owned_stop"] = layout.owned_stop
            for member_id in local_ensemble.member_ids:
                handle.create_dataset(
                    f"member_{member_id}",
                    data=np.asarray(local_ensemble.members[member_id]),
                    compression=None,
                )
            handle.flush()
        _fsync_path(shard_tmp)
        shard_bytes = int(shard_tmp.stat().st_size)
        os.replace(shard_tmp, shard_path)
        descriptor = {
            "file": shard_name,
            "owned_start": layout.owned_start,
            "owned_stop": layout.owned_stop,
            "member_ids": list(local_ensemble.member_ids),
            "dtype": next(iter(local_dtypes), None),
            "nbytes": shard_bytes,
        }
        if layout_fingerprint is not None:
            descriptor["layout_fingerprint"] = str(layout_fingerprint)
    except Exception as error:
        local_write_error = (
            f"rank {int(topology.world_rank)}: {type(error).__name__}: {error}"
        )
        try:
            shard_tmp.unlink(missing_ok=True)
        except OSError:
            pass

    write_results = topology.world.gather(
        (descriptor, local_write_error), root=0
    )
    commit_error = None
    if int(topology.world_rank) == 0:
        try:
            write_errors = [error for _, error in write_results if error is not None]
            if write_errors:
                raise OSError("; ".join(write_errors))
            descriptors = [descriptor for descriptor, _ in write_results]
            dtypes = {item["dtype"] for item in descriptors if item["dtype"] is not None}
            if len(dtypes) != 1:
                raise TypeError("checkpoint shards must use one state dtype")
            state_dtype = next(iter(dtypes))
            by_member: dict[int, list[tuple[int, int]]] = {}
            for item in descriptors:
                for member_id in item["member_ids"]:
                    by_member.setdefault(int(member_id), []).append(
                        (int(item["owned_start"]), int(item["owned_stop"]))
                    )
            expected_count = (
                int(sidecar["Nens"])
                if "Nens" in sidecar
                else (max(by_member) + 1 if by_member else 0)
            )
            if set(by_member) != set(range(expected_count)):
                raise ValueError("checkpoint shards do not cover every ensemble member")
            for member_id, owned in by_member.items():
                if tuple(sorted(owned)) != intervals:
                    raise ValueError(
                        f"checkpoint member {member_id} does not cover the global state"
                    )
            manifest = {
                "format": _FORMAT,
                "complete": True,
                "run_id": run_id,
                "timestep": timestep,
                "layout_id": layout.layout_id,
                "global_size": layout.global_size,
                "number_of_members": expected_count,
                "dtype": state_dtype,
                "source_topology": {
                    "ensemble_groups": int(topology.ensemble_groups),
                    "spatial_ranks": int(topology.spatial_ranks),
                },
                "metadata": sidecar,
                "shards": descriptors,
            }
            manifest_tmp = staging_dir / ".manifest.json.tmp"
            manifest_path = staging_dir / "manifest.json"
            with manifest_tmp.open("w", encoding="utf-8") as stream:
                json.dump(manifest, stream, indent=2, sort_keys=True)
                stream.flush()
                os.fsync(stream.fileno())
            os.replace(manifest_tmp, manifest_path)
            _fsync_path(staging_dir)
            os.replace(staging_dir, final_dir)
            _fsync_path(root)
        except Exception as error:  # make collective failure deterministic
            commit_error = f"{type(error).__name__}: {error}"
            shutil.rmtree(staging_dir, ignore_errors=True)
    commit_error = topology.world.bcast(commit_error, root=0)
    if commit_error is not None:
        raise RuntimeError(f"cannot commit distributed checkpoint: {commit_error}")
    return DistributedCheckpoint(final_dir, timestep, run_id, sidecar)


def _save_distributed_block_checkpoint(
    root: Path,
    timestep: int,
    local_ensemble: LocalMemberEnsemble,
    topology: Any,
    *,
    run_id: str,
    metadata: Mapping[str, Any],
    layout_fingerprint: str | None = None,
) -> DistributedCheckpoint:
    """Save segmented variable blocks without packing global state gaps."""

    layout = local_ensemble.layout
    if not isinstance(layout, DistributedBlockStateLayout):
        raise TypeError("block checkpoint requires DistributedBlockStateLayout")
    block_intervals = validate_block_spatial_partition(layout, topology.spatial_comm)
    final_dir = root / _checkpoint_name(timestep)
    staging_dir = root / f".{_checkpoint_name(timestep)}.{run_id}.staging"
    setup_error = None
    if int(topology.world_rank) == 0:
        try:
            root.mkdir(parents=True, exist_ok=True)
            if final_dir.exists():
                raise FileExistsError(f"checkpoint already exists: {final_dir}")
            if staging_dir.exists():
                shutil.rmtree(staging_dir)
            staging_dir.mkdir(parents=False)
        except Exception as error:
            setup_error = f"{type(error).__name__}: {error}"
    setup_error = topology.world.bcast(setup_error, root=0)
    if setup_error is not None:
        raise RuntimeError(f"cannot prepare distributed checkpoint: {setup_error}")
    topology.world.Barrier()

    shard_name = f"shard_{int(topology.world_rank):08d}.h5"
    shard_path = staging_dir / shard_name
    shard_tmp = staging_dir / f".{shard_name}.tmp"
    descriptor = None
    local_write_error = None
    block_descriptors = [
        {
            "name": block.name,
            "global_offset": block.global_offset,
            "global_size": block.global_size,
            "owned_start": block.owned_start,
            "owned_stop": block.owned_stop,
        }
        for block in layout.blocks
    ]
    try:
        local_dtypes = {
            np.asarray(local_ensemble.members[member_id]).dtype.str
            for member_id in local_ensemble.member_ids
        }
        if len(local_dtypes) > 1:
            raise TypeError("checkpoint members on one rank must share a dtype")
        with h5py.File(shard_tmp, "w") as handle:
            handle.attrs["format"] = _BLOCK_FORMAT
            handle.attrs["run_id"] = run_id
            handle.attrs["timestep"] = timestep
            handle.attrs["layout_id"] = layout.layout_id
            handle.attrs["global_size"] = layout.global_size
            for member_id in local_ensemble.member_ids:
                group = handle.create_group(f"member_{member_id}")
                member = np.asarray(local_ensemble.members[member_id])
                for block in layout.blocks:
                    group.create_dataset(
                        block.name,
                        data=member[layout.local_slice(block.name)],
                        compression=None,
                    )
            handle.flush()
        _fsync_path(shard_tmp)
        shard_bytes = int(shard_tmp.stat().st_size)
        os.replace(shard_tmp, shard_path)
        descriptor = {
            "file": shard_name,
            "member_ids": list(local_ensemble.member_ids),
            "blocks": block_descriptors,
            "dtype": next(iter(local_dtypes), None),
            "nbytes": shard_bytes,
        }
        if layout_fingerprint is not None:
            descriptor["layout_fingerprint"] = str(layout_fingerprint)
    except Exception as error:
        local_write_error = (
            f"rank {int(topology.world_rank)}: {type(error).__name__}: {error}"
        )
        try:
            shard_tmp.unlink(missing_ok=True)
        except OSError:
            pass

    write_results = topology.world.gather((descriptor, local_write_error), root=0)
    commit_error = None
    if int(topology.world_rank) == 0:
        try:
            write_errors = [error for _, error in write_results if error is not None]
            if write_errors:
                raise OSError("; ".join(write_errors))
            descriptors = [item for item, _ in write_results]
            dtypes = {item["dtype"] for item in descriptors if item["dtype"] is not None}
            if len(dtypes) != 1:
                raise TypeError("checkpoint shards must use one state dtype")
            state_dtype = next(iter(dtypes))
            member_ids = {
                int(member_id)
                for item in descriptors
                for member_id in item["member_ids"]
            }
            expected_count = int(metadata.get("Nens", 0))
            if expected_count <= 0:
                expected_count = max(member_ids) + 1 if member_ids else 0
            expected_members = set(range(expected_count))
            definitions = [
                {
                    "name": block.name,
                    "global_offset": block.global_offset,
                    "global_size": block.global_size,
                }
                for block in layout.blocks
            ]
            coverage: dict[tuple[int, str], list[tuple[int, int]]] = {}
            for item in descriptors:
                item_blocks = {entry["name"]: entry for entry in item["blocks"]}
                if set(item_blocks) != set(layout.block_names):
                    raise ValueError("checkpoint shard block definitions disagree")
                for definition in definitions:
                    entry = item_blocks[definition["name"]]
                    if (
                        int(entry["global_offset"]) != definition["global_offset"]
                        or int(entry["global_size"]) != definition["global_size"]
                    ):
                        raise ValueError("checkpoint shard block definitions disagree")
                for member_id in item["member_ids"]:
                    for entry in item["blocks"]:
                        coverage.setdefault((int(member_id), entry["name"]), []).append(
                            (int(entry["owned_start"]), int(entry["owned_stop"]))
                        )
            if member_ids != expected_members:
                raise ValueError("checkpoint shards do not cover every ensemble member")
            for member_id in expected_members:
                for block in layout.blocks:
                    if tuple(sorted(coverage.get((member_id, block.name), []))) != tuple(
                        block_intervals[block.name]
                    ):
                        raise ValueError(
                            f"checkpoint member {member_id} block {block.name} "
                            "does not cover the global block"
                        )
            manifest = {
                "format": _BLOCK_FORMAT,
                "complete": True,
                "run_id": run_id,
                "timestep": timestep,
                "layout_id": layout.layout_id,
                "global_size": layout.global_size,
                "number_of_members": expected_count,
                "dtype": state_dtype,
                "blocks": definitions,
                "source_topology": {
                    "ensemble_groups": int(topology.ensemble_groups),
                    "spatial_ranks": int(topology.spatial_ranks),
                },
                "metadata": dict(metadata),
                "shards": descriptors,
            }
            manifest_tmp = staging_dir / ".manifest.json.tmp"
            manifest_path = staging_dir / "manifest.json"
            with manifest_tmp.open("w", encoding="utf-8") as stream:
                json.dump(manifest, stream, indent=2, sort_keys=True)
                stream.flush()
                os.fsync(stream.fileno())
            os.replace(manifest_tmp, manifest_path)
            _fsync_path(staging_dir)
            os.replace(staging_dir, final_dir)
            _fsync_path(root)
        except Exception as error:
            commit_error = f"{type(error).__name__}: {error}"
            shutil.rmtree(staging_dir, ignore_errors=True)
    commit_error = topology.world.bcast(commit_error, root=0)
    if commit_error is not None:
        raise RuntimeError(f"cannot commit distributed checkpoint: {commit_error}")
    return DistributedCheckpoint(final_dir, timestep, run_id, metadata)


def load_distributed_checkpoint(
    checkpoint_path: str | os.PathLike[str],
    layout: DistributedStateLayout | DistributedBlockStateLayout,
    topology: Any,
    *,
    expected_run_id: str | None = None,
    number_of_members: int,
    expected_layout_fingerprint: str | None = None,
) -> tuple[LocalMemberEnsemble, DistributedCheckpoint]:
    """Load only target-local overlaps under the current MPI topology.

    When ``expected_layout_fingerprint`` is given and the checkpoint recorded
    fingerprints, shards with this rank's ownership must match it.  Older
    checkpoints without fingerprints load unchanged.
    """

    path = Path(checkpoint_path)
    with (path / "manifest.json").open("r", encoding="utf-8") as stream:
        manifest = json.load(stream)
    _check_layout_fingerprint(manifest, layout, expected_layout_fingerprint, topology)
    if isinstance(layout, DistributedBlockStateLayout):
        return _load_distributed_block_checkpoint(
            path,
            manifest,
            layout,
            topology,
            expected_run_id=expected_run_id,
            number_of_members=number_of_members,
        )
    if not isinstance(layout, DistributedStateLayout):
        raise TypeError("unsupported distributed checkpoint layout")
    if manifest.get("format") != _FORMAT or manifest.get("complete") is not True:
        raise ValueError("checkpoint manifest is incomplete or unsupported")
    if expected_run_id is not None and manifest.get("run_id") != str(expected_run_id):
        raise ValueError("checkpoint run_id does not match the requested run")
    if int(manifest.get("global_size", -1)) != layout.global_size:
        raise ValueError("checkpoint global state size is incompatible")
    if str(manifest.get("layout_id")) != layout.layout_id:
        raise ValueError("checkpoint state layout_id is incompatible")
    if int(manifest.get("number_of_members", -1)) != int(number_of_members):
        raise ValueError("checkpoint ensemble size is incompatible")

    member_ids = members_for_ensemble_slot(
        int(number_of_members), topology.ensemble_groups, topology.ensemble_slot
    )
    local_members: dict[int, np.ndarray] = {}
    target_start, target_stop = layout.owned_start, layout.owned_stop
    state_dtype = np.dtype(manifest.get("dtype", np.dtype(np.float64).str))
    for member_id in member_ids:
        local = np.empty(layout.owned_size, dtype=state_dtype)
        filled = np.zeros(layout.owned_size, dtype=bool)
        for shard in manifest["shards"]:
            if member_id not in [int(value) for value in shard["member_ids"]]:
                continue
            source_start = int(shard["owned_start"])
            source_stop = int(shard["owned_stop"])
            overlap_start = max(target_start, source_start)
            overlap_stop = min(target_stop, source_stop)
            if overlap_stop <= overlap_start:
                continue
            with h5py.File(path / shard["file"], "r") as handle:
                dataset = handle[f"member_{member_id}"]
                source_slice = slice(
                    overlap_start - source_start, overlap_stop - source_start
                )
                values = np.asarray(dataset[source_slice])
            if values.dtype != state_dtype:
                values = values.astype(state_dtype, copy=False)
            target_slice = slice(
                overlap_start - target_start, overlap_stop - target_start
            )
            if np.any(filled[target_slice]):
                raise ValueError("checkpoint state shards overlap")
            local[target_slice] = values
            filled[target_slice] = True
        if not np.all(filled):
            raise ValueError(f"checkpoint is incomplete for member {member_id}")
        local_members[member_id] = local

    checkpoint = DistributedCheckpoint(
        path=path,
        timestep=int(manifest["timestep"]),
        run_id=str(manifest["run_id"]),
        metadata=dict(manifest.get("metadata", {})),
    )
    return LocalMemberEnsemble(layout=layout, members=local_members), checkpoint


def _load_distributed_block_checkpoint(
    path: Path,
    manifest: Mapping[str, Any],
    layout: DistributedBlockStateLayout,
    topology: Any,
    *,
    expected_run_id: str | None,
    number_of_members: int,
) -> tuple[LocalMemberEnsemble, DistributedCheckpoint]:
    """Read only the member/block overlaps required by the target topology."""

    if manifest.get("format") != _BLOCK_FORMAT or manifest.get("complete") is not True:
        raise ValueError(
            "checkpoint manifest is incomplete or incompatible with block layout"
        )
    if expected_run_id is not None and manifest.get("run_id") != str(expected_run_id):
        raise ValueError("checkpoint run_id does not match the requested run")
    if int(manifest.get("global_size", -1)) != layout.global_size:
        raise ValueError("checkpoint global state size is incompatible")
    if str(manifest.get("layout_id")) != layout.layout_id:
        raise ValueError("checkpoint state layout_id is incompatible")
    if int(manifest.get("number_of_members", -1)) != int(number_of_members):
        raise ValueError("checkpoint ensemble size is incompatible")
    expected_blocks = [
        {
            "name": block.name,
            "global_offset": block.global_offset,
            "global_size": block.global_size,
        }
        for block in layout.blocks
    ]
    if manifest.get("blocks") != expected_blocks:
        raise ValueError("checkpoint state block definitions are incompatible")

    member_ids = members_for_ensemble_slot(
        int(number_of_members), topology.ensemble_groups, topology.ensemble_slot
    )
    local_members: dict[int, np.ndarray] = {}
    state_dtype = np.dtype(manifest.get("dtype", np.dtype(np.float64).str))
    for member_id in member_ids:
        local = np.empty(layout.owned_size, dtype=state_dtype)
        filled = np.zeros(layout.owned_size, dtype=bool)
        for block in layout.blocks:
            target_local = layout.local_slice(block.name)
            for shard in manifest["shards"]:
                if member_id not in [int(value) for value in shard["member_ids"]]:
                    continue
                source = next(
                    (entry for entry in shard["blocks"] if entry["name"] == block.name),
                    None,
                )
                if source is None:
                    raise ValueError("checkpoint shard is missing a state block")
                source_start = int(source["owned_start"])
                source_stop = int(source["owned_stop"])
                overlap_start = max(block.owned_start, source_start)
                overlap_stop = min(block.owned_stop, source_stop)
                if overlap_stop <= overlap_start:
                    continue
                with h5py.File(path / shard["file"], "r") as handle:
                    dataset = handle[f"member_{member_id}/{block.name}"]
                    values = np.asarray(
                        dataset[
                            overlap_start - source_start : overlap_stop - source_start
                        ]
                    )
                if values.dtype != state_dtype:
                    values = values.astype(state_dtype, copy=False)
                target_slice = slice(
                    target_local.start + overlap_start - block.owned_start,
                    target_local.start + overlap_stop - block.owned_start,
                )
                if np.any(filled[target_slice]):
                    raise ValueError("checkpoint state block shards overlap")
                local[target_slice] = values
                filled[target_slice] = True
        if not np.all(filled):
            raise ValueError(f"checkpoint is incomplete for member {member_id}")
        local_members[member_id] = local

    checkpoint = DistributedCheckpoint(
        path=path,
        timestep=int(manifest["timestep"]),
        run_id=str(manifest["run_id"]),
        metadata=dict(manifest.get("metadata", {})),
    )
    return LocalMemberEnsemble(layout=layout, members=local_members), checkpoint


def _owned_signature(layout_or_shard: Any) -> tuple[Any, ...]:
    """Owned intervals of a target layout or of a manifest shard entry."""

    if isinstance(layout_or_shard, Mapping):
        blocks = layout_or_shard.get("blocks")
        if blocks is not None:
            return tuple(
                (str(entry["name"]), int(entry["owned_start"]), int(entry["owned_stop"]))
                for entry in blocks
            )
        return (
            int(layout_or_shard["owned_start"]),
            int(layout_or_shard["owned_stop"]),
        )
    if isinstance(layout_or_shard, DistributedBlockStateLayout):
        return tuple(
            (block.name, int(block.owned_start), int(block.owned_stop))
            for block in layout_or_shard.blocks
        )
    return (int(layout_or_shard.owned_start), int(layout_or_shard.owned_stop))


def _check_layout_fingerprint(
    manifest: Mapping[str, Any],
    layout: Any,
    expected: str | None,
    topology: Any,
) -> None:
    """Refuse to load into a model whose DOF numbering differs from the writer.

    Global state rows are only meaningful under the numbering that wrote
    them.  When the writer recorded fingerprints, a shard with exactly this
    rank's ownership must carry the same fingerprint; if no shard has this
    ownership the decomposition changed and the numbering cannot be verified.
    A checkpoint written before fingerprints existed is accepted only
    without spatial decomposition on both sides, the one case in which an
    adapter that provides fingerprints has a reproducible numbering.
    """

    if expected is None:
        return
    shards = [
        shard for shard in manifest.get("shards", [])
        if shard.get("layout_fingerprint") is not None
    ]
    if not shards:
        source = int((manifest.get("source_topology") or {}).get("spatial_ranks", 1))
        target = int(getattr(topology, "spatial_ranks", 1))
        if source != 1 or target != 1:
            raise ValueError(
                "checkpoint has no layout fingerprint and was written or is being "
                f"loaded with spatial decomposition ({source} -> {target} spatial "
                "ranks); its degree-of-freedom numbering cannot be verified"
            )
        return
    mine = _owned_signature(layout)
    matching = [shard for shard in shards if _owned_signature(shard) == mine]
    if not matching:
        raise ValueError(
            "checkpoint was written with a different spatial decomposition; its "
            "degree-of-freedom numbering cannot be verified for this restart "
            "(use the original spatial rank count)"
        )
    for shard in matching:
        if str(shard["layout_fingerprint"]) != str(expected):
            raise ValueError(
                "checkpoint layout fingerprint does not match the current model "
                f"numbering (shard {shard.get('file')}); refusing to load state "
                "into a differently numbered mesh"
            )


class CheckpointValidationError(ValueError):
    """A checkpoint-looking directory is not a usable committed checkpoint."""


def checkpoint_step_from_name(name: str) -> int | None:
    """Return the step of a committed checkpoint directory name, else None.

    Only ``checkpoint_<8 digits>`` qualifies; dot-prefixed staging, deleting
    and quarantined directories never do.
    """

    match = _COMMITTED_NAME.fullmatch(str(name))
    return int(match.group(1)) if match else None


def _read_manifest(path: Path) -> dict[str, Any]:
    manifest_path = path / "manifest.json"
    if not manifest_path.is_file():
        raise CheckpointValidationError("manifest.json is missing")
    try:
        with manifest_path.open("r", encoding="utf-8") as stream:
            manifest = json.load(stream)
    except (OSError, ValueError) as error:
        raise CheckpointValidationError(f"manifest.json is unreadable: {error}") from error
    if not isinstance(manifest, dict):
        raise CheckpointValidationError("manifest.json is not an object")
    return manifest


def _require_tiling(intervals: Iterable[tuple[int, int]], size: int, what: str) -> None:
    covered = 0
    for start, stop in sorted((int(a), int(b)) for a, b in intervals if int(b) > int(a)):
        if start != covered:
            raise CheckpointValidationError(
                f"{what} rows [{covered}, {start}) are missing or overlap"
            )
        covered = stop
    if covered != int(size):
        raise CheckpointValidationError(f"{what} covers {covered} of {size} rows")


def validate_distributed_checkpoint(
    checkpoint_path: str | os.PathLike[str],
    *,
    expected_run_id: str | None = None,
    number_of_members: int | None = None,
    deep: bool = True,
) -> dict[str, Any]:
    """Verify one committed checkpoint and return its manifest.

    Checks the commit contract (name, ``complete`` manifest, supported
    format, matching step and run), that every member's state is covered
    exactly once by shard intervals, that every listed shard is a committed
    file (never a ``.tmp``) of the recorded size, and -- when ``deep`` --
    that each shard's HDF5 datasets exist with the recorded shapes.  Works
    for every checkpoint this module has ever written.  Raises
    :class:`CheckpointValidationError` describing the first problem.
    """

    path = Path(checkpoint_path)
    step = checkpoint_step_from_name(path.name)
    if step is None:
        raise CheckpointValidationError(f"{path.name!r} is not a committed checkpoint name")
    if not path.is_dir():
        raise CheckpointValidationError("checkpoint is not a directory")
    manifest = _read_manifest(path)
    checkpoint_format = manifest.get("format")
    if checkpoint_format not in _SUPPORTED_FORMATS:
        raise CheckpointValidationError(f"unsupported format {checkpoint_format!r}")
    if manifest.get("complete") is not True:
        raise CheckpointValidationError("manifest is not marked complete")
    try:
        if int(manifest["timestep"]) != step:
            raise CheckpointValidationError("manifest timestep disagrees with directory name")
        members = int(manifest["number_of_members"])
        global_size = int(manifest["global_size"])
        shards = list(manifest["shards"])
    except (KeyError, TypeError, ValueError) as error:
        if isinstance(error, CheckpointValidationError):
            raise
        raise CheckpointValidationError(f"manifest field is missing or invalid: {error}") from error
    if expected_run_id is not None and manifest.get("run_id") != str(expected_run_id):
        raise CheckpointValidationError(
            f"run_id {manifest.get('run_id')!r} differs from {str(expected_run_id)!r}"
        )
    if number_of_members is not None and members != int(number_of_members):
        raise CheckpointValidationError(
            f"checkpoint holds {members} members; expected {int(number_of_members)}"
        )
    if members <= 0 or not shards:
        raise CheckpointValidationError("checkpoint lists no members or shards")

    blocks = manifest.get("blocks") if checkpoint_format == _BLOCK_FORMAT else None
    if checkpoint_format == _BLOCK_FORMAT and not blocks:
        raise CheckpointValidationError("block checkpoint lists no state blocks")
    block_sizes = {str(entry["name"]): int(entry["global_size"]) for entry in blocks or []}

    seen_files: set[str] = set()
    coverage: dict[tuple[int, str | None], list[tuple[int, int]]] = {}
    for shard in shards:
        name = str(shard.get("file", ""))
        if not _SHARD_NAME.fullmatch(name) or name in seen_files:
            raise CheckpointValidationError(f"invalid or duplicate shard name {name!r}")
        seen_files.add(name)
        shard_path = path / name
        if not shard_path.is_file():
            raise CheckpointValidationError(f"shard {name} is missing")
        if "nbytes" in shard and shard_path.stat().st_size != int(shard["nbytes"]):
            raise CheckpointValidationError(
                f"shard {name} has {shard_path.stat().st_size} bytes; "
                f"manifest recorded {int(shard['nbytes'])}"
            )
        member_ids = [int(value) for value in shard.get("member_ids", [])]
        for member_id in member_ids:
            if not 0 <= member_id < members:
                raise CheckpointValidationError(f"shard {name} lists unknown member {member_id}")
            if blocks is None:
                coverage.setdefault((member_id, None), []).append(
                    (int(shard["owned_start"]), int(shard["owned_stop"]))
                )
            else:
                entries = {str(entry["name"]): entry for entry in shard.get("blocks", [])}
                if set(entries) != set(block_sizes):
                    raise CheckpointValidationError(f"shard {name} block list is inconsistent")
                for block_name, entry in entries.items():
                    coverage.setdefault((member_id, block_name), []).append(
                        (int(entry["owned_start"]), int(entry["owned_stop"]))
                    )
    for member_id in range(members):
        if blocks is None:
            _require_tiling(
                coverage.get((member_id, None), []), global_size, f"member {member_id}"
            )
        else:
            for block_name, size in block_sizes.items():
                _require_tiling(
                    coverage.get((member_id, block_name), []),
                    size,
                    f"member {member_id} block {block_name}",
                )

    if deep:
        for shard in shards:
            name = str(shard["file"])
            try:
                with h5py.File(path / name, "r") as handle:
                    if handle.attrs.get("format") != checkpoint_format or int(
                        handle.attrs.get("timestep", -1)
                    ) != step:
                        raise CheckpointValidationError(
                            f"shard {name} header disagrees with the manifest"
                        )
                    for member_id in shard.get("member_ids", []):
                        if blocks is None:
                            expected = {
                                f"member_{int(member_id)}": int(shard["owned_stop"])
                                - int(shard["owned_start"])
                            }
                        else:
                            expected = {
                                f"member_{int(member_id)}/{entry['name']}": int(
                                    entry["owned_stop"]
                                )
                                - int(entry["owned_start"])
                                for entry in shard["blocks"]
                            }
                        for dataset_name, length in expected.items():
                            dataset = handle.get(dataset_name)
                            if dataset is None or tuple(dataset.shape) != (length,):
                                raise CheckpointValidationError(
                                    f"shard {name} dataset {dataset_name} is missing "
                                    "or has the wrong shape"
                                )
            except CheckpointValidationError:
                raise
            except Exception as error:  # unreadable/corrupt HDF5
                raise CheckpointValidationError(
                    f"shard {name} is unreadable: {type(error).__name__}: {error}"
                ) from error
    return manifest


@dataclass(frozen=True)
class CheckpointDiscovery:
    """Result of :func:`discover_latest_checkpoint`."""

    path: Path | None
    timestep: int | None
    manifest: Mapping[str, Any] | None
    rejected: tuple[tuple[str, str], ...] = ()


def committed_checkpoint_paths(
    checkpoint_root: str | os.PathLike[str],
) -> list[tuple[int, Path]]:
    """Committed-looking checkpoint directories, newest first (unvalidated)."""

    root = Path(checkpoint_root)
    if not root.is_dir():
        return []
    found = []
    for entry in root.iterdir():
        step = checkpoint_step_from_name(entry.name)
        if step is not None and entry.is_dir():
            found.append((step, entry))
    return sorted(found, key=lambda item: item[0], reverse=True)


def _discover_local(
    root: Path,
    expected_run_id: str | None,
    number_of_members: int | None,
    deep: bool,
) -> CheckpointDiscovery:
    rejected = []
    for step, path in committed_checkpoint_paths(root):
        try:
            manifest = validate_distributed_checkpoint(
                path,
                expected_run_id=expected_run_id,
                number_of_members=number_of_members,
                deep=deep,
            )
        except CheckpointValidationError as error:
            rejected.append((path.name, str(error)))
            continue
        return CheckpointDiscovery(path, step, manifest, tuple(rejected))
    return CheckpointDiscovery(None, None, None, tuple(rejected))


def discover_latest_checkpoint(
    checkpoint_root: str | os.PathLike[str],
    *,
    expected_run_id: str | None = None,
    number_of_members: int | None = None,
    comm: Any = None,
    deep: bool = True,
) -> CheckpointDiscovery:
    """Return the highest-step checkpoint that passes full validation.

    Staging directories and ``.tmp`` files are never candidates; committed
    directories that fail :func:`validate_distributed_checkpoint` are
    skipped (and reported in ``rejected``) so a damaged newest checkpoint
    falls back to the previous valid one.  With ``comm``, rank 0 scans the
    filesystem once and broadcasts the result, so every rank restarts from
    the same checkpoint.
    """

    root = Path(checkpoint_root)
    if comm is None:
        return _discover_local(root, expected_run_id, number_of_members, deep)
    result = None
    error = None
    if comm.Get_rank() == 0:
        try:
            result = _discover_local(root, expected_run_id, number_of_members, deep)
        except Exception as exc:  # broadcast instead of deadlocking peers
            error = f"{type(exc).__name__}: {exc}"
    result, error = comm.bcast((result, error), root=0)
    if error is not None:
        raise RuntimeError(f"checkpoint discovery failed: {error}")
    return result


def remove_abandoned_checkpoint_artifacts(
    checkpoint_root: str | os.PathLike[str],
) -> list[Path]:
    """Delete uncommitted staging and half-deleted directories.

    Call only while no checkpoint writer is active (e.g. before a run
    starts).  Committed ``checkpoint_<step>`` directories are never touched:
    only dot-prefixed ``.checkpoint_*.staging``/``.deleting`` directories
    are removed, and any stray ``.tmp`` file lives inside one of them.
    """

    root = Path(checkpoint_root)
    removed: list[Path] = []
    if not root.is_dir():
        return removed
    for entry in root.iterdir():
        if (
            entry.name.startswith(".checkpoint_")
            and entry.name.endswith(_ABANDONED_SUFFIXES)
            and entry.is_dir()
        ):
            shutil.rmtree(entry, ignore_errors=True)
            removed.append(entry)
    return removed


def quarantine_checkpoint(checkpoint_path: str | os.PathLike[str]) -> Path:
    """Un-publish an invalid committed directory without deleting its data.

    The directory is renamed to a dot-prefixed ``.invalid`` name so it is no
    longer a restart candidate and does not block rewriting that step.
    """

    path = Path(checkpoint_path)
    target = path.with_name(f".{path.name}.invalid")
    counter = 1
    while target.exists():
        target = path.with_name(f".{path.name}.invalid{counter}")
        counter += 1
    os.replace(path, target)
    return target


def delete_committed_checkpoint(checkpoint_path: str | os.PathLike[str]) -> None:
    """Remove a committed checkpoint so that no partial copy stays visible.

    The atomic rename to ``.checkpoint_<step>.deleting`` un-publishes it
    first; a crash during the subsequent ``rmtree`` therefore leaves only a
    ``.deleting`` directory, which discovery ignores and
    :func:`remove_abandoned_checkpoint_artifacts` cleans up.
    """

    path = Path(checkpoint_path)
    if checkpoint_step_from_name(path.name) is None:
        raise ValueError(f"refusing to delete non-checkpoint path {path}")
    doomed = path.with_name(f".{path.name}.deleting")
    if doomed.exists():
        shutil.rmtree(doomed, ignore_errors=True)
    os.replace(path, doomed)
    shutil.rmtree(doomed, ignore_errors=True)
