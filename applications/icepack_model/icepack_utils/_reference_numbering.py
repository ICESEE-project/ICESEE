"""Map a spatially distributed Icepack layout onto the serial reference numbering.

Execution mode 3 generates the truth trajectory and synthetic observations on
a single-rank (``COMM_SELF``) reference mesh, so every observation id is a
row of the *serial* scalar-space numbering.  A model group with more than one
spatial rank partitions and renumbers the same mesh: its PETSc global row
``g`` is generally a different physical node from serial row ``g``.  Routing
serial observation ids directly onto distributed rows therefore compares
each observation with the model at the wrong location.

Both meshes are loaded from the same checkpoint, so every node has bitwise
identical coordinates in both numberings.  The functions here match owned
nodes to reference nodes by coordinate and translate observation ids.
"""

from __future__ import annotations

import numpy as np


def match_owned_nodes_to_reference(reference_coords, owned_coords, *, atol=1e-6):
    """Return the reference-node index of every owned node.

    Raises ``ValueError`` if an owned node has no reference node within
    ``atol`` (same length units as the coordinates) or if two owned nodes
    match the same reference node.
    """

    from scipy.spatial import cKDTree

    reference = np.asarray(reference_coords, dtype=np.float64)
    owned = np.asarray(owned_coords, dtype=np.float64)
    if reference.ndim != 2 or owned.ndim != 2 or reference.shape[1] != owned.shape[1]:
        raise ValueError(
            f"coordinate arrays must be (n, dim) with equal dim; got "
            f"{reference.shape} and {owned.shape}"
        )
    if owned.shape[0] == 0:
        return np.empty(0, dtype=np.int64)
    distance, index = cKDTree(reference).query(owned, k=1)
    if not np.all(distance <= atol):
        worst = int(np.argmax(distance))
        raise ValueError(
            "distributed mesh node has no matching reference node: "
            f"owned node {worst} at {owned[worst].tolist()} is "
            f"{float(distance[worst]):.3e} from the nearest reference node"
        )
    index = np.asarray(index, dtype=np.int64)
    if np.unique(index).size != index.size:
        raise ValueError("several owned nodes matched the same reference node")
    return index


def reference_to_layout_permutation(owned_start, owned_reference_index, comm, node_count):
    """Collective over the spatial communicator: build ``perm`` with
    ``perm[reference_node] = distributed global row``.

    Every rank of ``comm`` returns the same complete permutation.  It is
    O(node_count) integers, the size of one scalar field.
    """

    owned_reference_index = np.asarray(owned_reference_index, dtype=np.int64)
    pieces = comm.allgather((int(owned_start), owned_reference_index))
    node_count = int(node_count)
    perm = np.full(node_count, -1, dtype=np.int64)
    for start, reference_index in pieces:
        if reference_index.size and (
            reference_index.min() < 0 or reference_index.max() >= node_count
        ):
            raise ValueError("reference node index outside the reference mesh")
        rows = start + np.arange(reference_index.size, dtype=np.int64)
        if np.any(perm[reference_index] >= 0):
            raise ValueError("reference node owned by more than one distributed row")
        perm[reference_index] = rows
    if np.any(perm < 0):
        missing = int(np.count_nonzero(perm < 0))
        raise ValueError(
            f"{missing} reference nodes are not owned by any rank of the model group"
        )
    if np.unique(perm).size != node_count:
        raise ValueError("distributed rows do not form a permutation of the reference nodes")
    return perm


def translate_observation_ids(observation_ids, perm):
    """Translate variable-major serial state rows to distributed state rows.

    ``observation_ids`` index the packed ``[var0 | var1 | ...]`` vector whose
    blocks each hold ``len(perm)`` nodes.  The block (variable) is kept; only
    the node within the block is renumbered.  Order is preserved, so
    positions into canonical observation arrays remain valid.
    """

    ids = np.asarray(observation_ids, dtype=np.int64)
    node_count = int(np.asarray(perm).size)
    block, node = np.divmod(ids, node_count)
    return block * node_count + np.asarray(perm, dtype=np.int64)[node]
