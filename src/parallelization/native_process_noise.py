"""Reuse canonical process noise at a native model boundary.

Physical state stays distributed. A model-group leader generates one canonical
noise increment and scatters its owned entries, subject to a workspace limit.
This compatibility policy has an O(global state size) stochastic workspace;
it is not a fully distributed random-field generator. No physical member is
gathered. Initial fields are untouched.

On restart, deterministic AR history is replayed from seeds and the observation
schedule. No extra mutable state is required in native checkpoints. Scratch AR
files live in a private node-local temporary directory, never in shared output.
"""
from __future__ import annotations

import os
import hashlib
import json
import tempfile
import numpy as np

from ICESEE.src.parallelization._mpi_forecast_functions import (
    add_member_process_noise, process_noise_is_due,
)


def canonical_noise_identity(options):
    """Invalidate checkpoints from another stochastic policy/configuration."""
    keys = ("nd", "nt", "dt", "num_state_vars", "total_state_param_vars",
            "sig_Q", "base_seed", "Lx", "Ly", "random_field_method",
            "length_scale", "Q_rho", "coord_key_tolerance",
            "k_neighbors", "num_passes", "blend",
            "process_noise_schedule", "obs_index")
    payload = {name:options.get(name) for name in keys}
    encoded = json.dumps(payload,sort_keys=True,
                         default=lambda value:np.asarray(value).tolist())
    return "canonical-noise-v1-" + hashlib.sha256(encoded.encode()).hexdigest()[:16]


class NativeProcessNoiseAdapter:
    """Decorate a native adapter without changing its model capabilities."""

    def __init__(self, adapter, *, topology, canonical_rows, reference_coords,
                 icesee_kwargs, workspace_budget_mb=256):
        self.adapter = adapter
        self.comm = topology.spatial_comm
        self.rows = np.asarray(canonical_rows, dtype=np.int64)
        self.nd = int(icesee_kwargs['nd'])
        problem = None
        if self.rows.ndim != 1 or np.any(self.rows < 0) or np.any(self.rows >= self.nd):
            problem = 'native noise row map is outside the canonical state'
        elif self.nd * 8 * 10 > float(workspace_budget_mb) * 1024**2:
            problem = 'canonical native noise exceeds its stochastic workspace budget'
        problems = self.comm.allgather(problem)
        if any(problems):
            raise ValueError(next(p for p in problems if p))
        self.row_maps = self.comm.gather(self.rows, root=0)
        if self.comm.rank == 0:
            joined = np.concatenate(self.row_maps)
            problem = None if np.array_equal(np.sort(joined), np.arange(self.nd)) else (
                'native noise ownership must cover every canonical row exactly once')
        problem = self.comm.bcast(problem, root=0)
        if problem:
            raise ValueError(problem)
        self.options = dict(icesee_kwargs)
        self.options['mesh_coords'] = np.asarray(reference_coords)
        self.options['sub_rank'] = 0
        dt, nt = float(self.options['dt']), int(self.options['nt'])
        tau = max(500, max(200, dt))
        alpha = 1 - dt / tau
        if alpha <= 0 or alpha > 1:
            alpha = .5
        rho = np.sqrt((1/dt)*((1-alpha)**2)*(1/(nt-2*alpha-nt*alpha**2+2*alpha**(nt+1))))
        self.options.update(alpha=alpha, rho=rho)
        self._last = {}
        self._scratch = None

    def __getattr__(self, name):
        return getattr(self.adapter, name)

    def _increment(self, member_id, timestep, kwargs):
        if self._scratch is None:
            self._scratch = tempfile.TemporaryDirectory(
                prefix='icesee-native-noise-', dir=os.environ.get('TMPDIR') or None)
        options = dict(self.options)
        for name in ('obs_index', 'number_obs_instants', 'process_noise_schedule'):
            if name in kwargs:
                options[name] = kwargs[name]
        options.update(_modelrun_datasets=self._scratch.name,
                       data_path=self._scratch.name)
        # Reconstruct AR state at the saved checkpoint using exactly the same
        # generator and update order. Replay never changes physical fields.
        if member_id not in self._last:
            for k in range(timestep):
                options['k'] = k
                if process_noise_is_due(options, k):
                    add_member_process_noise(np.zeros(self.nd), member_id, options)
        elif timestep <= self._last[member_id]:
            raise ValueError('native process noise cannot advance a member twice at one step')
        options['k'] = timestep
        increment = add_member_process_noise(np.zeros(self.nd), member_id, options)
        self._last[member_id] = timestep
        return increment

    def forecast_native_member(self, member, timestep, *, topology, icesee_kwargs):
        self.adapter.forecast_native_member(member, timestep, topology=topology,
                                            icesee_kwargs=icesee_kwargs)
        if not process_noise_is_due(icesee_kwargs, timestep):
            return
        error, pieces = None, None
        if self.comm.rank == 0:
            try:
                increment = self._increment(member.member_id, int(timestep), icesee_kwargs)
                pieces = [increment[rows] for rows in self.row_maps]
            except Exception as exc:
                error = f'{type(exc).__name__}: {exc}'
        error = self.comm.bcast(error, root=0)
        if error:
            raise RuntimeError(f'native process-noise generation failed: {error}')
        local_increment = self.comm.scatter(pieces, root=0)
        owned = member.pack_owned()
        if owned.shape != local_increment.shape:
            raise ValueError('native noise map does not match the owned field layout')
        member.unpack_owned(owned + local_increment)

    def close(self):
        if self._scratch is not None:
            self._scratch.cleanup()
            self._scratch = None
