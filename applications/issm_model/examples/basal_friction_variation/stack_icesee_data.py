#!/usr/bin/env python3
import os, re, glob
import h5py
import numpy as np

FNAME_PATTERN = r'icesee_enkf_ens_(\d+)\.h5$'  # matches ..._0000.h5, ..._12.h5, etc.

def _extract_time(fname: str) -> int:
    m = re.search(FNAME_PATTERN, os.path.basename(fname))
    if not m:
        raise ValueError(f"Bad filename (no time index): {fname}")
    return int(m.group(1))

def _list_sorted_files(input_dir: str):
    files = glob.glob(os.path.join(input_dir, "icesee_enkf_ens_*.h5"))
    if not files:
        raise RuntimeError(f"No input files found in {input_dir}")
    files.sort(key=_extract_time)
    return files

def _infer_dataset_name(h5path: str, prefer="states") -> str:
    with h5py.File(h5path, "r") as f:
        if prefer in f:
            return prefer
        keys = [k for k in f.keys() if isinstance(f[k], h5py.Dataset)]
        if len(keys) == 1:
            return keys[0]
        raise RuntimeError(
            f"Could not infer dataset name in {h5path}. "
            f"Found: {keys}. Pass dset_name explicitly."
        )

# ---------- Option A: VDS (no data copy) ----------
def build_vds(input_dir: str,
              dset_name: str | None = None,
              out_file: str | None = None,
              fillvalue=np.nan) -> str:
    files = _list_sorted_files(input_dir)
    if dset_name is None:
        dset_name = _infer_dataset_name(files[0], prefer="states")
    if out_file is None:
        out_file = os.path.join(input_dir, "icesee_ensemble_data.h5")

    # Probe meta
    with h5py.File(files[0], "r") as f0:
        nd, nens = f0[dset_name].shape
        dtype = f0[dset_name].dtype

    nt = len(files)

    # Build VDS layout
    layout = h5py.VirtualLayout(shape=(nd, nens, nt), dtype=dtype)
    for t, f in enumerate(files):
        vsrc = h5py.VirtualSource(f, dset_name, shape=(nd, nens))
        layout[:, :, t] = vsrc

    os.makedirs(os.path.dirname(out_file) or ".", exist_ok=True)
    with h5py.File(out_file, "w", libver="latest") as fout:
        # dset_name='ensemble'
        fout.create_virtual_dataset(dset_name, layout, fillvalue=fillvalue)
        fout.attrs.update({
            "nd": nd, "nens": nens, "nt": nt,
            "stack_type": "VDS",
            "source_dir": os.path.abspath(input_dir),
            "dataset_name": dset_name
        })
        # Compute ensemble mean (iterate time slices lazily)
        mean_dset = fout.create_dataset(
            "ensemble_mean", shape=(nd, nt), dtype=np.float64,
            chunks=(nd, 1), fillvalue=np.nan
        )
        for t in range(nt):
            arr = fout[dset_name][:, :, t]
            mean_dset[:, t] = np.nanmean(arr, axis=1)

    return out_file

# ---------- Option B: materialized 3-D HDF5 ----------
def consolidate_h5(input_dir: str,
                   dset_name: str | None = None,
                   out_file: str | None = None,
                   compression: str = "gzip",
                   compression_opts: int = 4,
                   chunks: tuple[int,int,int] | None = None,
                   allow_missing: bool = False) -> str:
    files = _list_sorted_files(input_dir)
    if dset_name is None:
        dset_name = _infer_dataset_name(files[0], prefer="states")
    if out_file is None:
        out_file = os.path.join(input_dir, "icesee_ensemble_data.h5")

    # Probe meta
    with h5py.File(files[0], "r") as f0:
        nd, nens = f0[dset_name].shape
        dtype = f0[dset_name].dtype

    nt = len(files)
    if chunks is None:
        # Good default for time-wise appends and time-slice reads
        chunks = (nd, nens, 1)

    os.makedirs(os.path.dirname(out_file) or ".", exist_ok=True)
    with h5py.File(out_file, "w") as fout:
        # dset = fout.create_dataset(
        #     dset_name, shape=(nd, nens, nt), dtype=dtype,
        #     chunks=chunks, compression=compression,
        #     compression_opts=compression_opts,
        #     shuffle=True, fletcher32=True
        # )
        # dset_name='ensemble'
        dset = fout.create_dataset(
            'ensemble', shape=(nd, nens, nt), dtype=dtype,
            chunks=chunks, compression=compression,
            compression_opts=compression_opts,
            shuffle=True, fletcher32=True
        )
        mean_dset = fout.create_dataset(
            "ensemble_mean", shape=(nd, nt), dtype=np.float64,
            chunks=(nd, 1), compression=compression,
            compression_opts=compression_opts,
            shuffle=True, fletcher32=True
        )
        fout.attrs.update({
            "nd": nd, "nens": nens, "nt": nt,
            "stack_type": "materialized",
            "source_dir": os.path.abspath(input_dir),
            "dataset_name": dset_name
        })

        # Copy time slices, one file at a time (low memory)
        for t, fpath in enumerate(files):
            try:
                with h5py.File(fpath, "r") as fi:
                    arr = fi[dset_name][...]
            except Exception as e:
                if allow_missing:
                    arr = np.full((nd, nens), np.nan, dtype=dtype)
                else:
                    raise RuntimeError(f"Failed reading {fpath}: {e}") from e

            if arr.shape != (nd, nens):
                raise ValueError(f"Shape mismatch at {fpath}: {arr.shape} != {(nd, nens)}")

            dset[:, :, t] = arr
            mean_dset[:, t] = np.nanmean(arr, axis=1)  # (nd,) → store column

    return out_file

# ---------- Convenience entry point for your pipeline ----------
def finalize_stack(output_dir: str,
                   mode: str = "vds",
                   dset_name: str | None = "states",
                   **kwargs) -> str:
    """
    mode: 'vds' (no copy) or 'h5' (materialized).
    kwargs are passed to the underlying builder (e.g., allow_missing=True).
    """
    if mode.lower() == "vds":
        return build_vds(output_dir, dset_name=dset_name, **kwargs)
    elif mode.lower() in ("h5", "materialized"):
        return consolidate_h5(output_dir, dset_name=dset_name, **kwargs)
    else:
        raise ValueError("mode must be 'vds' or 'h5'")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="Stack ICESEE DA slices into (nd, nens, nt)")
    ap.add_argument("output_dir", help="Directory containing icesee_enkf_ens_*.h5")
    ap.add_argument("--mode", choices=["vds", "h5"], default="vds",
                    help="vds = virtual dataset (no copy), h5 = materialized 3D file")
    ap.add_argument("--dset-name", default="states",
                    help="Dataset name inside each file (default: states)")
    ap.add_argument("--out-file", default=None, help="Path of output file")
    ap.add_argument("--allow-missing", action="store_true",
                    help="Fill missing/unreadable slices with NaN for materialized mode")
    ap.add_argument("--compression", default="gzip", help="HDF5 compression (h5 mode)")
    ap.add_argument("--compression-opts", type=int, default=4, help="Compression level (h5 mode)")
    ap.add_argument("--chunk-nd", type=int, default=None)
    ap.add_argument("--chunk-nens", type=int, default=None)
    ap.add_argument("--chunk-nt", type=int, default=1)
    args = ap.parse_args()

    chunks = None
    if args.chunk_nd or args.chunk_nens or args.chunk_nt:
        chunks = (args.chunk_nd, args.chunk_nens, args.chunk_nt)

    out = finalize_stack(
        args.output_dir, mode=args.mode, dset_name=args.dset_name,
        out_file=args.out_file, compression=args.compression,
        compression_opts=args.compression_opts, chunks=chunks,
        allow_missing=args.allow_missing
    )
    print(f"[Finalize] Stacked dataset written: {out}")

