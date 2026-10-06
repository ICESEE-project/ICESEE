"""Rank-sharded, bounded-memory checkpoints for execution mode 3.

The first backend deliberately uses ordinary serial HDF5 files: every MPI rank
writes its own shard and no rank gathers a complete member.  A manifest is
published only after every shard is durable.  Restart is independent of rank
placement because shards are addressed by stable member IDs and global state
intervals rather than by the rank that created them.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
import os
from pathlib import Path
import shutil
from typing import Any, Mapping

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


def _checkpoint_name(timestep: int) -> str:
    return f"checkpoint_{int(timestep):08d}"


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
    complete manifest are intentionally ignored.  This lets a restarted job
    recover from the last durable cycle after a writer exception or scheduler
    termination without treating a partial write as model state.
    """

    root = Path(checkpoint_root)
    if not root.is_dir():
        return None
    candidates: list[tuple[int, Path]] = []
    for path in root.glob("checkpoint_*"):
        manifest_path = path / "manifest.json"
        if not path.is_dir() or not manifest_path.is_file():
            continue
        try:
            with manifest_path.open("r", encoding="utf-8") as stream:
                manifest = json.load(stream)
            if manifest.get("format") not in _SUPPORTED_FORMATS or manifest.get("complete") is not True:
                continue
            if expected_run_id is not None and manifest.get("run_id") != str(
                expected_run_id
            ):
                continue
            timestep = int(manifest["timestep"])
            if path.name != _checkpoint_name(timestep):
                continue
        except (OSError, ValueError, TypeError, KeyError, json.JSONDecodeError):
            continue
        candidates.append((timestep, path))
    return max(candidates, default=(None, None), key=lambda item: item[0])[1]


def save_distributed_checkpoint(
    checkpoint_root: str | os.PathLike[str],
    timestep: int,
    local_ensemble: LocalMemberEnsemble,
    topology: Any,
    *,
    run_id: str,
    metadata: Mapping[str, Any] | None = None,
) -> DistributedCheckpoint:
    """Atomically save distributed state without reconstructing any member.

    ``metadata`` is the model-agnostic sidecar for random-stream keys,
    observation-cycle state, inflation/localization settings, and inversion
    schedule state.  Model-specific restart data should be represented by
    stable identifiers or written by a model adapter into its own shard.
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
        os.replace(shard_tmp, shard_path)
        descriptor = {
            "file": shard_name,
            "owned_start": layout.owned_start,
            "owned_stop": layout.owned_stop,
            "member_ids": list(local_ensemble.member_ids),
            "dtype": next(iter(local_dtypes), None),
        }
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
            os.replace(staging_dir, final_dir)
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
        os.replace(shard_tmp, shard_path)
        descriptor = {
            "file": shard_name,
            "member_ids": list(local_ensemble.member_ids),
            "blocks": block_descriptors,
            "dtype": next(iter(local_dtypes), None),
        }
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
            os.replace(staging_dir, final_dir)
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
) -> tuple[LocalMemberEnsemble, DistributedCheckpoint]:
    """Load only target-local overlaps under the current MPI topology."""

    path = Path(checkpoint_path)
    with (path / "manifest.json").open("r", encoding="utf-8") as stream:
        manifest = json.load(stream)
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
