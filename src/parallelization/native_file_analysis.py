"""Whole-member native forecasts with the shared file-backed analysis operator.

This is an explicit Mode-3 analysis policy, not dispatch to the Mode-2 runner.
It reuses localization, inference and publication mathematics while native
members keep their model-owned state and stable member identities. One spatial
rank/member is required; the spatially partitioned native policy is separate.
"""
import numpy as np
from .EnKF_parallel_io import EnKF_fully_parallel_IO
from ._mpi_ensemble_intialization import apply_initial_spread_filebacked
from .distributed_runtime import LocalMemberEnsemble
from .distributed_native_cycle import NativeCycleResult


class FileBackedNativeAnalysis:
    def initialize(self, pool, adapter, *, topology, icesee_kwargs):
        if topology.spatial_ranks != 1:
            raise NotImplementedError('Whole-member analysis policy requires spatial_ranks=1')
        kw = icesee_kwargs
        world = topology.world
        kw.update(comm_world=world, subcomm=topology.spatial_comm,
                  rank_world=world.rank, sub_rank=0, color=topology.ensemble_slot,
                  global_shape=int(kw['nd']), dim_list=[int(kw['nd'])]*world.size)
        dt = float(kw['dt']); nt = int(kw['nt'])
        tau = max(500, max(200, dt))
        alpha = 1-dt/tau
        if alpha <= 0 or alpha > 1:
            alpha = 0.5
        kw['alpha'] = alpha
        kw['rho'] = np.sqrt((1/dt)*((1-alpha)**2)/(nt-2*alpha-nt*alpha**2+2*alpha**(nt+1)))
        kw.setdefault('time_analysis_file_writing',0.0)
        kw.setdefault('time_analysis_ensemble_mean_generation',0.0)
        self.io = EnKF_fully_parallel_IO(
            'native_analysis', int(kw['nd']), int(kw['Nens']), nt,
            topology.spatial_comm, world, kw, serial_file_creation=True,
            base_path=kw['data_path'], batch_size=1)
        self.io.H_matrix(**kw)
        self.io._ensure_batch(0)
        initial = pool.snapshot_owned()
        for member_id, vector in initial.members.items():
            self.io.write_forecast(0, vector, member_id)
        world.Barrier()
        apply_initial_spread_filebacked(self.io, kw.get('initial_spread_factor'))
        self._restore(pool, 0)

    def _restore(self, pool, step):
        with self.io._read_only_state(step) as dataset:
            members = {i: np.asarray(dataset[:,i],dtype=np.float64)
                       for i in pool.member_ids}
        restored = LocalMemberEnsemble(pool.layout,members)
        pool.restore_owned(restored)
        return restored

    def __call__(self, pool, adapter, timestep, observation_batches, *,
                 number_of_batches, topology, icesee_kwargs, error_mode, **unused):
        kw = icesee_kwargs
        forecast = pool.forecast_with_adapter(timestep,adapter,topology=topology,
                                              icesee_kwargs=kw)
        self.io._ensure_batch(timestep+1)
        for member_id, vector in forecast.members.items():
            self.io.write_forecast(timestep+1,vector,member_id)
        topology.world.Barrier()
        if number_of_batches:
            self.io.compute_analysis_update(**kw)
            analyzed = self._restore(pool,timestep+1)
            rows = sum(len(b.observation_ids) for b in observation_batches)
        else:
            analyzed = forecast
            rows = 0
        return NativeCycleResult(forecast,analyzed, self.io.last_analysis_transform if number_of_batches else np.eye(int(kw['Nens'])),rows)

    def close(self):
        self.io.close()
