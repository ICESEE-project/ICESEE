# Deploying ICESEE on workstations and HPC systems

This guide separates the scientific source from its compiled numerical environment. A working Firedrake, PETSc, MPI and HDF5 installation can often be retained while ICESEE Python source is updated. Every run must record both identities. A container recipe revision or a Spack parent revision alone does not identify the scientific code actually imported.

The commands below were checked against repository contents on 10 November 2026. No dependency stack was rebuilt, container image pulled, or production allocation submitted during this documentation audit. Installation and import probes establish different things from numerical validation; use the [testing guide](testing-and-contributing.md) and [operations guide](operations.md) before expanding a run.

## 1 Choose and record the deployment path

| Situation | Practical route | What must be retained |
|---|---|---|
| Working lab HPC installation | Activate that exact environment and select an explicit ICESEE source checkout | Environment activation path, lock file, compiler/MPI details and scientific source hash |
| New native HPC installation | ICESEE-Spack reviewed revision, installed in a new directory | Parent revision, all submodule revisions, actual concretization and installed packages |
| Container-enabled HPC site | Approved Apptainer image plus explicit source/input/output binds | SIF hash, recipe revision, embedded source identity, host launcher and PMIx compatibility |
| Local development | Existing model-compatible environment, or a container appropriate to the model | Architecture, image digest or complete Python/native stack identity |
| Core Lorenz learning exercises | Python 3.11+ and ICESEE dependencies with functional MPI where needed | Source/import location, Python packages and test logs |

Do not copy a working environment's directory as a portable installation: absolute paths and native libraries may be embedded. Preserve its manifests and rebuild separately when migration is needed. Keep the approved installation available until a replacement has passed the same application checks.

### Audited revisions

| Component | Audited revision | Consequence |
|---|---|---|
| ICESEE scientific baseline | `3748e416d3a8120cc860119bd7224c6f9ad41ff7` | Core numerical/PIG evidence baseline |
| ICESEE teaching/report snapshot | `2fa088a1f8bba8f0e11f1dfc2f94b9f77fde1289` | Includes the new teaching files and later ISSM lifecycle corrections |
| ICESEE `main` at audit | `3d7e6e757e07b80dcca1fdd65dffb8ca8853316a` | Different source; do not infer it contains the same runtime repairs |
| ICESEE-Spack `main` | `68e37fffcd13cc60636028d051d4637e7f573d4f` | Uses Git submodules and legacy Firedrake defaults |
| ICESEE-Spack legacy tag | `v0.1.0-icepack-legacy` -> `c512f4b4a253e9342cf7366d0e226f1304c8c301` | Older installer release; inspect its own scripts and submodule/tree layout |
| ICESEE-Spack modern tag | `v0.2.0-icepack-modern` -> `d40e51b8151695116e705acae77a2b05823557d9` | Experimental dependency track, not evidence of current PIG validation |
| ICESEE-Containers `main` | `884cc3046a26b5999f6db16a98ad38a96bfdbcfd` | Multiple distinct recipes; recipe and image must both be identified |

The Spack parent at the audited revision records `ICESEE/` at `a793cd198ab39eedc668f3f67e08d56ba9bf22dc` and `spack/` at `f36409b78591f8b02e8332eb4ad78da62c02571e`. Its bundled ICESEE is therefore not this handbook's scientific baseline. The `.gitmodules` branch hint `main` does not move that recorded pin during an ordinary submodule checkout. Earlier local documentation describing these directories as subtrees predates the current layout.

## 2 Select the scientific source explicitly

Use a fresh checkout for an experiment that must remain reproducible. The commands below select the complete teaching/report snapshot, including the ISSM corrections. The earlier scientific baseline above predates the new tutorial files; use it only when deliberately reproducing that earlier source. Later documentation-only commits can add prose while retaining this executable snapshot.

```bash
git clone https://github.com/ICESEE-project/ICESEE.git ICESEE
cd ICESEE
git checkout --detach 2fa088a1f8bba8f0e11f1dfc2f94b9f77fde1289
git status --short
git rev-parse HEAD
export ICESEE_SOURCE_ROOT="$PWD"
```

A detached checkout is useful for execution because a later `git pull` cannot silently move its branch. Develop changes on a named branch. Do not use `git reset --hard`, force submodule updates, or reinstall dependencies to resolve an import-path mistake.

For an existing working environment, activate it first. Then prepend the directory containing the `ICESEE` package to `PYTHONPATH`. In the example above the checkout directory is named `ICESEE`, and its parent makes the root package importable:

```bash
export PYTHONPATH="$(dirname "$ICESEE_SOURCE_ROOT"):$ICESEE_SOURCE_ROOT:${PYTHONPATH:-}"
python - <<'CHECK'
from pathlib import Path
import importlib.util, os, sys
spec = importlib.util.find_spec("ICESEE")
assert spec is not None, "ICESEE cannot be resolved"
locations = [Path(p).resolve() for p in (spec.submodule_search_locations or [])]
expected = Path(os.environ["ICESEE_SOURCE_ROOT"]).resolve()
print("python", sys.executable)
print("ICESEE locations", locations)
assert expected in locations, (expected, locations)
CHECK
```

This path check does not import Firedrake or launch MPI. If it selects an unexpected source, correct activation and path order before running. `pip install -e` is optional for a working source-based environment and can resolve additional dependencies; retain that environment's constraints rather than allowing an unconstrained upgrade over its native stack. The repository's Python metadata declares version `0.1.9`; the earlier paper labels its archived release `v1.0.0`. Use commit hashes until maintainers reconcile release metadata.

## 3 Reuse ICESEE-Spack safely

Start a fresh Bash shell or batch environment so activation has an explicit scope:

```bash
export ICESEE_SPACK_ROOT=/absolute/path/to/approved/ICESEE-Spack
source "$ICESEE_SPACK_ROOT/scripts/activate.sh"
command -v python
command -v mpirun
command -v mpicc
python -m pip check
```

At the audited Spack revision, `activate.sh` activates `.spack-env/icesee`, selects Spack Open MPI, may activate `venv-firedrake`, prepends the bundled source, and attempts MATLAB module loading. It also contains `module purge` and a GCC 13 `libstdc++` preload. These side effects mean its output and final environment matter more than the environment present before sourcing. Restore the explicit scientific-source path afterward and rerun the import-resolution check above. Save `module list`, paths and `LD_PRELOAD` in private run provenance if relevant.

Do not rerun `scripts/install.sh --with-icepack` merely to obtain new ICESEE Python source: `build_firedrake.sh` removes and recreates an existing `venv-firedrake`. An environment rebuild is a separate operation with an installation log and application revalidation.

### A new native stack

Use a new installation directory and a site-approved build allocation with enough storage and walltime. The current installer expects GCC, G++, GFortran, Git, downloads and the site's module environment. For ISSM, arrange a compatible MATLAB installation and license; a successful Icepack installation does not establish ISSM readiness.

```bash
git clone https://github.com/ICESEE-project/ICESEE-Spack.git ICESEE-Spack-new
cd ICESEE-Spack-new
git checkout --detach 68e37fffcd13cc60636028d051d4637e7f573d4f
bash scripts/submodules.sh --with-icepack
# On the site-approved build node or allocation:
JOBS=8 bash scripts/install.sh --with-icepack
source scripts/activate.sh
```

For an ISSM build use `--with-issm` as well, after checking the site's MATLAB module and compiler configuration. The installer selects system Spack if it is already on `PATH`, otherwise the pinned `spack/` tree; record which executable it actually used. `activate.sh` subsequently uses the pinned tree. This difference deserves inspection rather than an assumption that both commands selected one Spack installation.

At audited `main`, the requested native stack includes Python 3.11, NumPy 1.26.4, Open MPI 5.0, HDF5 1.14.5 with MPI, PETSc/petsc4py 3.23.4, and a Firedrake 2025.4.2 legacy constraint set. These are requested versions, not a claim that every package in an existing environment has those versions. In particular `build_icepack.sh` installs `h5py>3.12.1` through pip; check its resulting MPI capability explicitly. Preserve `.spack-env/icesee/spack.yaml`, its generated `spack.lock` when present, `spack find -lv`, submodule status, installation logs and `python -m pip freeze`.

Existing Icepack/ISSM checkouts may be deliberately reused at their current commits. Do not set `ICESEE_SUBMODULES_FORCE=1` for routine source updates: it changes dependency revisions and requires rebuilding affected components. The Spack repository's [submodule guide](https://github.com/ICESEE-project/ICESEE-Spack/blob/68e37fffcd13cc60636028d051d4637e7f573d4f/docs/submodules.md) explains the migration from old subtree layouts and preservation of ignored build files.

## 4 Require functional MPI and HDF5 probes

The audited `scripts/test.sh` is diagnostic, not a strict acceptance gate: it catches import failures and continues, and its ICESEE probe imports uppercase `ICESEE` but then references lowercase `icesee`. Read each result. Require direct checks that terminate nonzero on failure:

```bash
python - <<'CHECK'
import sys, h5py, numpy, ICESEE
from mpi4py import MPI
print("python", sys.executable)
print("ICESEE", list(ICESEE.__path__))
print("numpy", numpy.__version__)
print("MPI", MPI.Get_library_version())
print("h5py", h5py.__version__, "MPI enabled", h5py.get_config().mpi)
assert h5py.get_config().mpi, "This deployment requires MPI-enabled h5py"
CHECK
```

Run the following under the same launcher that the application will use, on allocated compute resources. Replace `mpiexec -n 2` with the site's approved equivalent if required:

```bash
mpiexec -n 2 python -c 'from mpi4py import MPI; c=MPI.COMM_WORLD; r=c.allgather(c.rank); print(c.rank,c.size,r,flush=True); assert c.size == 2 and r == [0,1]'
```

The required result is one world of size two, not two independent singleton processes. Matching vendor strings alone cannot prove launcher compatibility. For a shared-file probe, create a new test script in a writable shared scratch directory:

```python
# mpi_hdf5_probe.py: run with exactly two ranks on an allocated resource.
from pathlib import Path
import sys
import h5py
import numpy as np
from mpi4py import MPI
comm = MPI.COMM_WORLD
assert comm.size == 2
assert h5py.get_config().mpi
path = Path(sys.argv[1])
if comm.rank == 0 and path.exists():
    comm.Abort(2)  # Never overwrite an existing artifact.
comm.Barrier()
with h5py.File(path, "w", driver="mpio", comm=comm) as handle:
    ds = handle.create_dataset("rank", (comm.size,), dtype="i8")
    ds[comm.rank] = comm.rank
comm.Barrier()
if comm.rank == 0:
    with h5py.File(path, "r") as handle:
        np.testing.assert_array_equal(handle["rank"][:], [0, 1])
    print("parallel HDF5 probe passed", path)
```

```bash
mpiexec -n 2 python mpi_hdf5_probe.py /absolute/shared/scratch/new-mpi-probe.h5
```

All ranks in the file communicator participate in file creation and metadata creation. Setting `HDF5_MPI=ON` at runtime does not convert a serial h5py wheel into an MPI build; it is a build-time setting. See the [h5py parallel-HDF5 instructions](https://docs.h5py.org/en/stable/mpi.html). Never disable file locking as a generic hang workaround.

For Icepack add `python -c 'import firedrake, icepack; from petsc4py import PETSc; print(PETSc.Sys.getVersion())'`, followed by the small application regression. A successful import does not exercise distributed point location, solver convergence or restart.

## 5 Identify container recipes and immutable images

Clone and pin the recipe repository independently of ICESEE:

```bash
git clone https://github.com/ICESEE-project/ICESEE-Containers.git
cd ICESEE-Containers
git checkout --detach 884cc3046a26b5999f6db16a98ad38a96bfdbcfd
cd spack-managed/combined-container
```

The real combined Docker filenames at this revision are `Dockerfile.matlab-runtime`, `Dockerfile.no-matlab-runtime`, and `Dockerfile.cloudruns`. Some examples in the combined README use `Dockerfile.nomatlab-runtime` or other spellings that do not exist in this directory. Use the filenames present in the checked-out source.

The combined external-MATLAB Apptainer recipe is `combined-env-external-matlab.def`; the bundled-MATLAB recipe is `combined-env-inbuilt-matlab.def`. These bootstrap published Docker tags. A recipe hash therefore does not freeze its base image. Many full builder recipes also clone moving Spack, ICESEE `main`, and Icepack branches. Rebuilding the same recipe later may produce different software. Preserve the built SIF and its SHA-256, or use immutable OCI digests and archived build manifests. Do not equate a `v1.0` image tag with ICESEE's source version.

On a permitted build host, a reviewed external-MATLAB recipe can be built with:

```bash
apptainer build combined-external-matlab.sif combined-env-external-matlab.def
sha256sum combined-external-matlab.sif
apptainer inspect combined-external-matlab.sif
```

The build may require site-supported fakeroot or a remote build service. That command was verified against the recipe path, not executed in this audit. The full external-MATLAB Docker builder still uses a MathWorks MATLAB image in its builder stage even though the final runtime uses Ubuntu 24.04; plan build access accordingly.

`Dockerfile.cloudruns` pins its base digest and records a previously qualified OCI image in comments. Its recorded tag/digest is provenance in source, not a newly verified registry pull. It uses Linux amd64 in its documented build command. Record the actual architecture and image digest used on the target site.

### Use the wrapper matching the application

Combined Apptainer images install `with-firedrake`, `with-icepack`, `with-icesee`, and `with-issm`. The Icepack/ICESEE/Firedrake side uses Open MPI, whereas the ISSM side uses a separate MPICH stack. Selecting `with-issm` is not interchangeable with `with-icepack` for a Python MPI job. The smaller standalone `icepack-lean-env.def` recipe has different activation mechanics; do not assume every image includes the combined wrappers.

Probe the actual combined image before binding new source:

```bash
export ICESEE_IMAGE=/absolute/path/to/approved-combined.sif
apptainer exec "$ICESEE_IMAGE" with-icepack python -c 'import sys,ICESEE,icepack,h5py; from mpi4py import MPI; print(sys.executable,list(ICESEE.__path__),MPI.Get_library_version(),h5py.get_config().mpi)'
```

For a current source bind, bind the checkout to `/opt/ICESEE`, because the combined wrapper prepends `/opt` to Python's package path. Binding source changes the imported code, not the installed binary dependencies:

```bash
export ICESEE_SOURCE_ROOT=/absolute/path/to/ICESEE
apptainer exec   --bind "$ICESEE_SOURCE_ROOT:/opt/ICESEE:ro"   "$ICESEE_IMAGE" with-icepack python -c 'import ICESEE; print(list(ICESEE.__path__))'
```

A read-only source bind is suitable for import probes and runs whose outputs and input caches have been redirected to writable locations. Some examples expect to write beside the script or create data in their `data/` directory. Supply those required files ahead of time and bind dedicated writable input-cache/output locations; otherwise use a disposable writable source checkout rather than making the approved source mutable. Verify the selected Python environment has the new source's required dependencies.

### MPI across container boundaries

Use the site's supported hybrid or bind MPI model, as described in [Apptainer's MPI guide](https://apptainer.org/docs/user/latest/mpi.html). Host launcher, PMI/PMIx integration, container MPI and interconnect libraries must agree. `srun --mpi=pmix` is a site-specific choice, not a universal fix. Do not load an arbitrary host MPICH module over an Open MPI Python stack merely because an old example names it.

After host/container MPI coordination is agreed, run the two-rank world probe through that exact launch path:

```bash
# Inside an allocated Slurm resource, only where PMIx is the site's supported path:
srun --mpi=pmix --ntasks=2 apptainer exec   --bind "$ICESEE_SOURCE_ROOT:/opt/ICESEE:ro"   "$ICESEE_IMAGE" with-icepack python -c 'from mpi4py import MPI; c=MPI.COMM_WORLD; r=c.allgather(c.rank); print(c.rank,c.size,r,flush=True); assert c.size == 2 and r == [0,1]'
```

Repeat the parallel-HDF5 probe with a shared output bind. Then run the small model case. A C MPI hello-world probe is useful but does not prove that the selected Python interpreter, mpi4py, HDF5 and model use the same MPI libraries.

### MATLAB and writable ISSM paths

For an external-MATLAB combined image, bind the approved host MATLAB installation to `/opt/matlab/R2024b` and set `MATLABROOT` explicitly. Supply the site's license configuration as an environment setting; do not put license server addresses or credentials in portable run scripts. Bind writable ISSM execution directories separately from the read-only application source. Check the wrapper's MATLAB startup/search path using `issmversion` before a coupled forecast. The current wrapper adds MATLAB flags itself; inspect it before combining `-r` and `-batch` options. No GNU Octave equivalence is established by these recipes.

## 6 Run the PACE scripts from the correct source

The PIG source includes `pace/resolve_environment.sh`, `pace/launch_rank.sh`, `pace/runtime_check.sbatch`, `pace/validate_mode3_24.sbatch`, and `run_idealized_1.sbatch`. The native launchers expect a retained Spack installation; they are not generic container launchers. Use a separately reviewed container batch wrapper for container execution.

The explicit roots are:

```bash
export ICESEE_SOURCE_ROOT=/absolute/path/to/approved/ICESEE
export ICESEE_SPACK_ROOT=/absolute/path/to/approved/ICESEE-Spack
cd "$ICESEE_SOURCE_ROOT/applications/icepack_model/examples/idealized_pig"
```

Before submission, inspect the `#SBATCH` account, partition, node count, ranks, memory and walltime. The checked-in values target a specific lab allocation and are not portable authorization. `runtime_check.sbatch` requests two nodes/four MPI ranks, two members with two ranks per model, and a 0.1-year run plus a completed-checkpoint resume. `validate_mode3_24.sbatch` requests one node/eight ranks, 24 members, two ranks per model and the five-year configuration. Source inputs under `data/` must be available; the resolver checks the layout but does not download the scientific dataset for you.

```bash
# Submit only after editing site account/partition in a personal wrapper as needed:
sbatch pace/runtime_check.sbatch
```

The current native process grid is `P = Pe * Ps`, with `Ps = model_nprocs`; it schedules multiple ensemble members per model group when needed. Choose ranks per node divisible by `Ps` and block placement for node-local groups. `Nens` is the number of ensemble members, not a CPU count. Budget MATLAB/ISSM subprocesses and threads separately from MPI ranks.

The PIG rank launcher places JIT caches in node-local temporary storage and makes XDG state rank-private, with job/step namespaces. The container recipes' default persistent cache helper may instead choose shared `/scratch`; review and adapt that choice for the target filesystem. Avoid shared SQLite state across MPI ranks. Preserve final and checkpoint outputs in persistent storage, not node-local caches.

## 7 Promotion and rollback

Maintain a small deployment ledger containing approved source hash, environment/image identity, input hashes, command/YAML, acceptance evidence and a pointer to the previous approved bundle. After any compiler, MPI, PETSc, Firedrake, Icepack, container or scheduler change, rerun the relevant probes and small application/restart tests before a long simulation. After a Python-only ICESEE update, still run its clock, analysis and adapter regressions: source compatibility does not prove scientific equivalence.

Use new output directories when promoting the corrected `state-columns-v2` clock. Its initial state is at column zero; forecast cycle `k` produces column `k+1`, and observations are analyzed at their matching state column. Old-clock histories cannot be relabeled or silently resumed into the new clock. Preserve them with their original source identity.

Rollback selects the previous approved source and environment and its matching checkpoint/output contract. Do not run old code over new outputs or alter checkpoint identity checks to force reuse. A scheduler exit zero, a readable plot or finite values alone is insufficient: inspect step exits, model times, complete member/ownership coverage, physical extrema and restart parity as described in the [operations guide](operations.md).

## Source and reference trail

- [ICESEE runtime source at the documented commit](https://github.com/ICESEE-project/ICESEE/tree/3748e416d3a8120cc860119bd7224c6f9ad41ff7)
- [ICESEE-Spack source at the audited parent](https://github.com/ICESEE-project/ICESEE-Spack/tree/68e37fffcd13cc60636028d051d4637e7f573d4f)
- [ICESEE-Containers recipes at the audited parent](https://github.com/ICESEE-project/ICESEE-Containers/tree/884cc3046a26b5999f6db16a98ad38a96bfdbcfd)
- [Firedrake installation](https://www.firedrakeproject.org/install.html) and [Icepack installation](https://icepack.github.io/install/)
- [Apptainer MPI](https://apptainer.org/docs/user/latest/mpi.html), [bind mounts](https://apptainer.org/docs/user/latest/bind_paths_and_mounts.html), and [h5py parallel HDF5](https://docs.h5py.org/en/stable/mpi.html)
