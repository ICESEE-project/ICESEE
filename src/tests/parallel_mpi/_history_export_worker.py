"""Small real-MPI checkpoint-to-aggregate round trip; no model backend."""
import hashlib
import json
from pathlib import Path
import sys

import h5py
import numpy as np
from mpi4py import MPI

from ICESEE.src.parallelization.distributed_adapter import DistributedBlockStateLayout, DistributedStateBlockLayout
from ICESEE.src.parallelization.distributed_history_export import finalize_ensemble_history, prepare_reference_layout
from ICESEE.src.parallelization.distributed_runtime import LocalMemberEnsemble, members_for_ensemble_slot
from ICESEE.src.parallelization.distributed_topology import create_distributed_topology
from ICESEE.src.parallelization.mode3_checkpointing import CheckpointPolicy, Mode3CheckpointManager

world = MPI.COMM_WORLD
data = Path(sys.argv[1])
topology = create_distributed_topology(world, spatial_ranks=int(sys.argv[2]))
definitions = [('h', 0, 7), ('forcing', 7, 3)]
blocks = tuple(DistributedStateBlockLayout(name, offset, size,
               size * topology.spatial_rank // topology.spatial_ranks,
               size * (topology.spatial_rank + 1) // topology.spatial_ranks)
               for name, offset, size in definitions)
layout = DistributedBlockStateLayout(blocks, layout_id='export-mpi-layout')
manager = Mode3CheckpointManager(data / '_mode3_state_history', run_id='export-mpi',
    topology=topology, number_of_members=5, nt=2, analysis_steps=[1],
    policy=CheckpointPolicy(every=1, keep_last=0), time_grid=[0, 0.05, 0.1])
permutations = {name: np.arange(size - 1, -1, -1) for name, _, size in definitions}
prepare_reference_layout(world, manager.history_root, run_id='export-mpi',
                         layout_id=layout.layout_id, row_maps=permutations)


def ensemble(column):
    members = members_for_ensemble_slot(5, topology.ensemble_groups, topology.ensemble_slot)
    return LocalMemberEnsemble(layout, {
        member: np.concatenate([100 * np.arange(block.global_offset + block.owned_start,
                                               block.global_offset + block.owned_stop, dtype=float)
                                + 10 * member + column for block in blocks])
        for member in members})


manager.write_initial(ensemble(0))
for step in range(2):
    manager.commit_step(step, ensemble(step + 1), did_analysis=(step == 1),
                        analyses_completed=int(step == 1))
world.Barrier()
output = finalize_ensemble_history(world, data, row_chunk_size=2)
error = None
if world.Get_rank() == 0:
    try:
        reference_rows = np.concatenate([offset + permutations[name] for name, offset, _ in definitions])
        expected = (100 * reference_rows[:, None, None]
                    + 10 * np.arange(5)[None, :, None] + np.arange(3)[None, None, :])
        with h5py.File(output, 'r') as handle:
            np.testing.assert_array_equal(handle['ensemble'][:], expected)
            np.testing.assert_array_equal(handle['ensemble_mean'][:], expected.mean(axis=1))
            assert not handle['ensemble'].is_virtual and handle.attrs['complete']
        before = hashlib.sha256(output.read_bytes()).hexdigest()
        missing = data / '_mode3_state_history/steps/checkpoint_00000000'
        missing.rename(missing.with_name('.removed-step'))
    except Exception as exc:
        error = str(exc)
error = world.bcast(error, root=0)
assert error is None, error
world.Barrier()
try:
    finalize_ensemble_history(world, data, row_chunk_size=2)
except RuntimeError as exc:
    refused = 'missing history snapshot' in str(exc)
else:
    refused = False
assert all(world.allgather(refused)), 'all ranks must reject partial history'
if world.Get_rank() == 0:
    assert hashlib.sha256(output.read_bytes()).hexdigest() == before
    print(json.dumps(dict(world_size=world.Get_size(), spatial_ranks=topology.spatial_ranks,
                          exact=True, missing_history_refused=True)), flush=True)
