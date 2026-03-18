from stack_icesee_data import finalize_stack
# from mpi4py import MPI

# comm = MPI.COMM_WORLD()

# if comm.rank == 0:  # or run on a post-hook
# Option A: no-copy, instant
# out_vds = finalize_stack("_modelrun_datasets", mode="vds", dset_name="states")
# print("VDS ready:", out_vds)

# Option B: portable single file
out_h5 = finalize_stack("_modelrun_datasets", mode="h5", dset_name="states",
                        allow_missing=False, compression="gzip", compression_opts=4)
print("Materialized file:", out_h5)

