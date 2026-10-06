from __future__ import annotations

import os
import sys
import tempfile
import types
from pathlib import Path

import numpy as np
import pytest

# The real ``_icepack_model.py`` import chain (pulled in transitively below)
# has three standalone-script assumptions that only hold when it is imported
# from inside ``examples/idealized_pig`` with that directory on ``sys.path``:
# ``modelfunc`` is imported as a bare top-level package, and
# ``config._utility_imports`` parses ``sys.argv`` with ``argparse``, loads
# ``params.yaml`` relative to the current working directory, AND deletes/
# recreates ``icesee_kwargs['data_path']`` -- all at import time. An explicit
# ``--data_path`` override keeps that auto-clean guard off idealized_pig's
# own real (default, relative) ``_modelrun_datasets`` directory; see
# ``test_lorenz_distributed_native_adapter.py`` for the concrete damage this
# gap caused once for Lorenz96's tracked CI figure before this fix was added
# everywhere. Satisfy all three here so this test runs safely under plain
# ``pytest`` from the repo root like every other test, without special
# invocation flags.
_REPO_ROOT = Path(__file__).resolve().parents[2]
_IDEALIZED_PIG_DIR = _REPO_ROOT / "applications" / "icepack_model" / "examples" / "idealized_pig"
for _extra_path in (str(_REPO_ROOT), str(_REPO_ROOT.parent), str(_IDEALIZED_PIG_DIR)):
    if _extra_path not in sys.path:
        sys.path.insert(0, _extra_path)

_SAFE_DATA_PATH = tempfile.mkdtemp(prefix="icesee_icepack_idealized_pig_native_test_")
_ARGV_BACKUP = sys.argv[:]
_CWD_BACKUP = os.getcwd()
sys.argv = [sys.argv[0], "--data_path", _SAFE_DATA_PATH]
os.chdir(_IDEALIZED_PIG_DIR)
try:
    from ICESEE.applications.icepack_model.examples.idealized_pig import _icepack_native as native
finally:
    sys.argv = _ARGV_BACKUP
    os.chdir(_CWD_BACKUP)
# ICESEE.-prefixed (not the bare "applications...."/"src...." forms): a
# real mode-3 run reproduced "native member returned an unsupported
# distributed layout" traced to exactly this kind of split import --
# _icepack_native.py/_distributed_fields.py were fixed to import
# consistently via the ICESEE.-prefixed path, so this test must match or
# its own isinstance()/identity checks below split the same way.
from ICESEE.applications.icepack_model.icepack_utils._distributed_fields import (
    IcepackForecastFields,
    IcepackNativeState,
)
from ICESEE.src.parallelization.distributed_native_runtime import (
    validate_native_distributed_adapter,
)
from ICESEE.src.run_model_da._error_generation import coordinate_keyed_white_noise
from ICESEE.src.utils.random_streams import initialization_seed


class _FakeDat:
    def __init__(self, values):
        self.values = np.asarray(values, dtype=float)

    @property
    def data(self):
        return self.values

    @property
    def data_ro(self):
        return self.values


class _FakeVec:
    """Stand-in for the PETSc layout vector read by ``scalar_ownership``."""

    def __init__(self, size):
        self._size = size

    def getSize(self):
        return self._size

    def getOwnershipRange(self):
        return (0, self._size)


class _FakeSpace:
    def __init__(self, size):
        self.size = size
        self.dof_dset = types.SimpleNamespace(layout_vec=_FakeVec(size))


class _FakeCoordExpr:
    """Stand-in for one component (x or y) of a fake SpatialCoordinate."""

    def __init__(self, values):
        self.values = np.asarray(values, dtype=float)


class _FakeFunction:
    """Stand-in for ``firedrake.Function`` that just wraps a numpy array."""

    def __init__(self, space):
        self.space = space
        self.dat = _FakeDat(np.zeros(space.size))

    def function_space(self):
        return self.space

    def interpolate(self, expr):
        if isinstance(expr, _FakeCoordExpr):
            self.dat.values = np.array(expr.values, dtype=float)
        return self


@pytest.fixture(autouse=True)
def _patch_function(monkeypatch):
    monkeypatch.setattr(native, "Function", _FakeFunction)


def _scalar_field(size=2):
    return _FakeFunction(_FakeSpace(size))


def _vector_field(size=2, components=2):
    field = _FakeFunction(_FakeSpace(size))
    field.dat.values = np.zeros((size, components))
    return field


def _topology(spatial_comm="spatial-comm-a"):
    return types.SimpleNamespace(spatial_comm=spatial_comm)


def _base_icesee_kwargs(**overrides):
    kwargs = dict(
        initFile="init.h5",
        meshFile="mesh.exp",
        SMBFile="smb.tif",
        uThresh=300.0,
        dt=0.2,
        GLThresh=41,
        hThresh=30,
        water_to_ice=1.087241003,
        wrong_basal_melt_field=2.0,
        # basal_melt_field's participation in the packed/registered state
        # is driven by vec_inputs (2026-09-28 reconciliation -- see
        # _icepack_native._basal_melt_in_state), not a joint_estimation
        # flag. Default here matches the current icesee/features
        # convention (5 state vars, basal_melt_field included).
        vec_inputs=["h", "u", "v", "s", "basal_melt_field"],
        base_seed=42,
    )
    kwargs.update(overrides)
    return kwargs


def _no_basal_melt_kwargs(**overrides):
    """A state convention where basal_melt_field is NOT declared at all --
    the adapter must omit it from the packed state entirely."""

    overrides.setdefault("vec_inputs", ["h", "u", "v", "s"])
    return _base_icesee_kwargs(**overrides)


def _install_fake_mesh_pipeline(monkeypatch, *, build_calls):
    """Patch initializeMesh/initializeRun/icepack solver construction with
    lightweight fakes and record how many times the mesh is (re)built."""

    Q = _FakeSpace(4)
    V = _FakeSpace(8)
    h0 = _FakeFunction(Q)
    h0.dat.values[:] = [10.0, 11.0, 12.0, 13.0]
    s0 = _FakeFunction(Q)
    s0.dat.values[:] = [1.0, 2.0, 3.0, 4.0]
    u0 = _FakeFunction(V)
    bed = _FakeFunction(Q)
    floating0 = _FakeFunction(Q)
    A0 = _FakeFunction(Q)
    beta0 = _FakeFunction(Q)
    smb = _FakeFunction(Q)

    # Fake "locally owned" physical coordinates for Q's 4 fake DOFs -- an
    # arbitrary but fixed (x, y) per row, standing in for whatever subset
    # of the real mesh's coordinates this rank happens to own.
    coord_x = np.array([500.0, 1500.0, 2500.0, 3500.0])
    coord_y = np.array([-100.0, -200.0, -300.0, -400.0])
    coords = np.stack((coord_x, coord_y), axis=1)

    def fake_initialize_mesh(**mesh_kwargs):
        assert "comm" in mesh_kwargs
        return ("mesh-sentinel", {}, Q, V)

    def fake_spatial_coordinate(mesh):
        assert mesh == "mesh-sentinel"
        return _FakeCoordExpr(coord_x), _FakeCoordExpr(coord_y)

    def fake_initialize_run(mesh_kwargs, forward_solver, mesh, Q_, V_):
        build_calls.append(mesh_kwargs["comm"])
        return (
            None,  # h (unused)
            h0,
            None,  # s (unused)
            s0,
            u0,
            bed,
            None,  # zF (unused)
            None,  # grounded0 (unused)
            floating0,
            A0,
            beta0,
            smb,
        )

    class _FakeIceStream:
        def __init__(self, friction=None, viscosity=None):
            self.friction = friction
            self.viscosity = viscosity

    class _FakeFlowSolver:
        def __init__(self, model, **opts):
            self.model = model
            self.opts = opts

    fake_icepack = types.SimpleNamespace(
        models=types.SimpleNamespace(IceStream=_FakeIceStream),
        solvers=types.SimpleNamespace(FlowSolver=_FakeFlowSolver),
    )

    monkeypatch.setattr(native, "initializeMesh", fake_initialize_mesh)
    monkeypatch.setattr(native, "initializeRun", fake_initialize_run)
    monkeypatch.setattr(native, "icepack", fake_icepack)
    monkeypatch.setattr(native, "SpatialCoordinate", fake_spatial_coordinate)
    return dict(
        Q=Q, V=V, h0=h0, s0=s0, u0=u0, bed=bed, floating0=floating0, A0=A0,
        beta0=beta0, smb=smb, coords=coords,
    )


def test_shared_context_is_built_once_per_topology_and_scoped_to_spatial_comm(monkeypatch):
    native._SHARED_CONTEXT_CACHE.clear()
    build_calls = []
    _install_fake_mesh_pipeline(monkeypatch, build_calls=build_calls)
    icesee_kwargs = _base_icesee_kwargs()

    topo_a = _topology("comm-a")
    ctx1 = native._shared_context(topo_a, icesee_kwargs)
    ctx2 = native._shared_context(topo_a, icesee_kwargs)
    assert ctx1 is ctx2
    assert build_calls == ["comm-a"]

    topo_b = _topology("comm-b")
    ctx3 = native._shared_context(topo_b, icesee_kwargs)
    assert ctx3 is not ctx1
    assert build_calls == ["comm-a", "comm-b"]


@pytest.mark.parametrize(
    "step,dt,expected_scenario",
    [
        (0, 0.2, "control"),
        (10, 1.0, "warm"),  # 6 <= 10 < 15
        (16, 1.0, "control"),  # 15 <= 16 < 18
        (19, 1.0, "warm"),  # 18 <= 19 < 20
        (22, 1.0, "control"),  # 20 <= 22 < 25
        (26, 1.0, "warm"),  # 25 <= 26 < 27
        (29, 1.0, "control"),  # 27 <= 29 < 31
        (35, 1.0, "warm"),  # 31 <= 35 < 40
        (45, 1.0, "control"),  # 40 <= 45 < 48
        (49, 1.0, "warm"),  # 48 <= 49 < 50
        (55, 1.0, "control"),  # 50 <= 55 < 59
        (62, 1.0, "warm"),  # 59 <= 62 < 64
        (67, 1.0, "control"),  # 64 <= 67 < 69
        (72, 1.0, "warm"),  # 69 <= 72 < 76
        (80, 1.0, "control"),  # 76 <= 80
    ],
)
def test_basal_melt_schedule_matches_reference_scenario_timeline(
    monkeypatch, step, dt, expected_scenario
):
    seen = {}

    def fake_basal_melt_rate(icesee_kwargs, step_, floating, Q, s, h, scenario, experiment=None):
        seen["scenario"] = scenario
        seen["experiment"] = experiment
        return "melt-field", 0.0

    monkeypatch.setattr(native, "BasalMeltRate", fake_basal_melt_rate)
    icesee_kwargs = _base_icesee_kwargs(dt=dt)
    native._basal_melt_for_step(icesee_kwargs, step, "floating", "Q", "s", "h")
    assert seen["scenario"] == expected_scenario
    # Mode-3 always forecasts the ensemble/forecast trajectory (never the
    # true trajectory, which is generated once, execution-mode-independently,
    # by _icepack_enkf.generate_true_state).
    assert seen["experiment"] == native.EXPERIMENT_WRONG


def test_initialize_member_uses_deterministic_coordinate_keyed_perturbation(monkeypatch):
    """Same base seed + same logical member + same physical coordinate =>
    same perturbation, independent of decomposition/ownership -- the
    scientific invariant this fix establishes (replacing the old
    array-position-keyed, decomposition-dependent draw)."""
    native._SHARED_CONTEXT_CACHE.clear()
    ctx_objs = _install_fake_mesh_pipeline(monkeypatch, build_calls=[])

    icepack_calls = []

    def fake_icepack(solver, h, u, smb, basal_melt_field, bed, dt, h0, run_kwargs):
        icepack_calls.append(
            dict(h=h, u=u, smb=smb, basal_melt_field=basal_melt_field,
                 bed=bed, dt=dt, h0=h0, run_kwargs=run_kwargs)
        )
        return (_scalar_field(), _vector_field(), _scalar_field(), "floating-next", "grounded-next")

    def fake_basal_melt_rate(icesee_kwargs, step, floating, Q, s, h, scenario, experiment=None):
        assert scenario == "control"
        assert s is ctx_objs["s0"]
        assert h is ctx_objs["h0"]
        return "melt-field", 0.0

    monkeypatch.setattr(native, "Icepack", fake_icepack)
    monkeypatch.setattr(native, "BasalMeltRate", fake_basal_melt_rate)

    icesee_kwargs = _no_basal_melt_kwargs()
    topo = _topology()

    # Reproduce the expected coordinate-keyed perturbation draw
    # independently: (base_seed, member_id) -> seed (same derivation
    # modes 1/2 use), then coordinate_keyed_white_noise over this rank's
    # OWNED physical coordinates only -- NOT np.random.normal over local
    # array size, which was the pre-fix (decomposition-dependent) draw.
    seed = initialization_seed(icesee_kwargs["base_seed"], member_id=7, variable_index=0)
    expected_noise = coordinate_keyed_white_noise(ctx_objs["coords"], seed)
    expected_perturbation = 15 * expected_noise

    state = native.initialize_member(7, topology=topo, icesee_kwargs=icesee_kwargs)

    assert isinstance(state, IcepackNativeState)
    np.testing.assert_allclose(
        icepack_calls[0]["h"].dat.data_ro,
        ctx_objs["h0"].dat.data_ro + expected_perturbation,
    )
    assert icepack_calls[0]["u"] is ctx_objs["u0"]
    assert icepack_calls[0]["smb"] is ctx_objs["smb"]
    assert icepack_calls[0]["bed"] is ctx_objs["bed"]
    assert icepack_calls[0]["h0"] is ctx_objs["h0"]
    assert icepack_calls[0]["basal_melt_field"] == "melt-field"
    assert state.basal_melt is None  # basal_melt_field not in vec_inputs


def _run_initialize_member(monkeypatch, member_id, *, icesee_kwargs=None):
    """Shared boilerplate for the perturbation-focused tests below: install
    the fake mesh/Icepack pipeline, run one member through
    ``initialize_member``, and return (ctx_objs, perturbed_h_array).

    Captures the actual "h" (perturbed thickness Function) Icepack() is
    called with -- NOT ctx.h0, which is a separate, never-mutated
    baseline reference the real code passes alongside it."""

    ctx_objs = _install_fake_mesh_pipeline(monkeypatch, build_calls=[])
    icepack_calls = []

    def fake_icepack(solver, h, u, smb, basal_melt_field, bed, dt, h0, run_kwargs):
        icepack_calls.append(h)
        # u/s must share h's OWN ownership size (Q=4 here) --
        # build_icepack_field_registry derives "nodal" ownership from
        # thickness (h) and validates u/v/s against it.
        return (
            h, _vector_field(size=h.dat.data_ro.size), _scalar_field(size=h.dat.data_ro.size),
            "floating-next", "grounded-next",
        )

    def fake_basal_melt_rate(icesee_kwargs, step, floating, Q, s, h, scenario, experiment=None):
        return "melt-field", 0.0

    monkeypatch.setattr(native, "Icepack", fake_icepack)
    monkeypatch.setattr(native, "BasalMeltRate", fake_basal_melt_rate)

    kwargs = icesee_kwargs if icesee_kwargs is not None else _no_basal_melt_kwargs()
    topo = _topology()
    native.initialize_member(member_id, topology=topo, icesee_kwargs=kwargs)
    return ctx_objs, np.array(icepack_calls[0].dat.data_ro)


def test_initialize_member_perturbation_uses_only_locally_owned_coordinates(monkeypatch):
    """The coordinate array passed to coordinate_keyed_white_noise must be
    exactly this rank's OWNED coordinates (same size as the local h0) --
    never a larger/global-sized array -- confirming the fix generates only
    locally owned values and never materializes a global Nx-sized
    temporary."""
    native._SHARED_CONTEXT_CACHE.clear()
    ctx_objs = _install_fake_mesh_pipeline(monkeypatch, build_calls=[])

    captured = {}

    def fake_coordinate_keyed_white_noise(coords, seed, coord_tolerance=1e-6):
        captured["coords"] = np.asarray(coords)
        captured["seed"] = seed
        return np.zeros(np.asarray(coords).shape[0])

    monkeypatch.setattr(native, "coordinate_keyed_white_noise", fake_coordinate_keyed_white_noise)
    monkeypatch.setattr(
        native, "Icepack",
        lambda solver, h, u, smb, basal_melt_field, bed, dt, h0, run_kwargs: (
            h, _vector_field(size=h.dat.data_ro.size), _scalar_field(size=h.dat.data_ro.size),
            "floating-next", "grounded-next",
        ),
    )
    monkeypatch.setattr(
        native, "BasalMeltRate",
        lambda icesee_kwargs, step, floating, Q, s, h, scenario, experiment=None: ("melt-field", 0.0),
    )

    icesee_kwargs = _no_basal_melt_kwargs()
    topo = _topology()
    native.initialize_member(7, topology=topo, icesee_kwargs=icesee_kwargs)

    assert "coords" in captured
    np.testing.assert_array_equal(captured["coords"], ctx_objs["coords"])
    assert captured["coords"].shape[0] == ctx_objs["h0"].dat.data_ro.size


def test_initialize_member_perturbation_reproducible_for_same_member(monkeypatch):
    """Calling initialize_member twice for the SAME member_id (fresh
    shared-context build each time) must produce the identical perturbed
    thickness -- reproducibility of the member-keyed stream."""
    native._SHARED_CONTEXT_CACHE.clear()
    _ctx1, h_first = _run_initialize_member(monkeypatch, 7)
    native._SHARED_CONTEXT_CACHE.clear()
    _ctx2, h_second = _run_initialize_member(monkeypatch, 7)
    np.testing.assert_array_equal(h_first, h_second)


def test_initialize_member_perturbation_differs_across_members(monkeypatch):
    """Different member_ids must receive distinct perturbation
    realizations at the same physical coordinates."""
    native._SHARED_CONTEXT_CACHE.clear()
    _ctx1, h_member7 = _run_initialize_member(monkeypatch, 7)
    native._SHARED_CONTEXT_CACHE.clear()
    _ctx2, h_member8 = _run_initialize_member(monkeypatch, 8)
    assert not np.allclose(h_member7, h_member8)


def test_initialize_member_perturbation_amplitude_is_exactly_15(monkeypatch):
    """h = h0 + 15 * N(0,1): the amplitude constant must still apply
    exactly, unchanged by the coordinate-keying fix."""
    native._SHARED_CONTEXT_CACHE.clear()
    ctx_objs = _install_fake_mesh_pipeline(monkeypatch, build_calls=[])

    known_noise = np.array([0.1, -0.2, 0.3, -0.4])

    def fake_coordinate_keyed_white_noise(coords, seed, coord_tolerance=1e-6):
        return known_noise

    monkeypatch.setattr(native, "coordinate_keyed_white_noise", fake_coordinate_keyed_white_noise)

    icepack_calls = []

    def fake_icepack(solver, h, u, smb, basal_melt_field, bed, dt, h0, run_kwargs):
        icepack_calls.append(h)
        return (
            h, _vector_field(size=h.dat.data_ro.size), _scalar_field(size=h.dat.data_ro.size),
            "floating-next", "grounded-next",
        )

    monkeypatch.setattr(native, "Icepack", fake_icepack)
    monkeypatch.setattr(
        native, "BasalMeltRate",
        lambda icesee_kwargs, step, floating, Q, s, h, scenario, experiment=None: ("melt-field", 0.0),
    )

    icesee_kwargs = _no_basal_melt_kwargs()
    topo = _topology()
    native.initialize_member(7, topology=topo, icesee_kwargs=icesee_kwargs)

    np.testing.assert_array_equal(
        icepack_calls[0].dat.data_ro,
        ctx_objs["h0"].dat.data_ro + 15 * known_noise,
    )


def test_initialize_member_basal_melt_present_when_declared_in_state(monkeypatch):
    """2026-09-28 reconciliation (second pass): basal_melt_field's presence
    in the packed state is driven by vec_inputs, not joint_estimation.
    With joint_estimation left unset/False (the current Idealized PIG
    config), the joint-estimation-specific `+ wrong_basal_melt_field`
    nudge does not fire, so initialize_member returns whatever
    BasalMeltRate produced, unmodified -- see
    test_initialize_member_nudges_basal_melt_when_joint_estimation for the
    joint_estimation=True case, where the nudge correctly still fires."""
    native._SHARED_CONTEXT_CACHE.clear()
    ctx_objs = _install_fake_mesh_pipeline(monkeypatch, build_calls=[])

    def fake_icepack(solver, h, u, smb, basal_melt_field, bed, dt, h0, run_kwargs):
        return (_scalar_field(), _vector_field(), _scalar_field(), "floating-next", "grounded-next")

    # initialize_member now passes BasalMeltRate's own return value straight
    # through as state.basal_melt (no more wrapping via a nudge step), so
    # the fake must itself satisfy build_icepack_field_registry's
    # scalar_ownership() introspection -- a real _FakeFunction, not a bare
    # stand-in class.
    melt_field = _scalar_field(size=4)
    melt_field.dat.values[:] = [1.0, 2.0, 3.0, 4.0]

    seen = {}

    def fake_basal_melt_rate(icesee_kwargs, step, floating, Q, s, h, scenario, experiment=None):
        seen["experiment"] = experiment
        return melt_field, 0.0

    monkeypatch.setattr(native, "Icepack", fake_icepack)
    monkeypatch.setattr(native, "BasalMeltRate", fake_basal_melt_rate)

    icesee_kwargs = _base_icesee_kwargs(wrong_basal_melt_field=2.0)  # vec_inputs includes basal_melt_field
    topo = _topology()
    state = native.initialize_member(0, topology=topo, icesee_kwargs=icesee_kwargs)

    assert state.basal_melt is not None
    # Unmodified -- joint_estimation is unset/False here, so the
    # joint-estimation-specific nudge does not fire.
    np.testing.assert_allclose(state.basal_melt.dat.data_ro, [1.0, 2.0, 3.0, 4.0])
    assert seen["experiment"] == native.EXPERIMENT_WRONG


def test_initialize_member_nudges_basal_melt_when_joint_estimation(monkeypatch):
    """2026-09-28 reconciliation (second pass): joint_estimation means "are
    we performing joint state-parameter estimation," not "does
    basal_melt_field exist." The `+ wrong_basal_melt_field` nudge is
    genuinely joint-estimation-algorithm-specific (it initializes the
    ensemble's WRONG prior belief about a jointly-estimated parameter) and
    must still fire when joint_estimation=True, regardless of whether
    basal_melt_field's existence is (separately, correctly) driven by
    vec_inputs."""
    native._SHARED_CONTEXT_CACHE.clear()
    ctx_objs = _install_fake_mesh_pipeline(monkeypatch, build_calls=[])

    def fake_icepack(solver, h, u, smb, basal_melt_field, bed, dt, h0, run_kwargs):
        return (_scalar_field(), _vector_field(), _scalar_field(), "floating-next", "grounded-next")

    melt_field = _scalar_field(size=4)
    melt_field.dat.values[:] = [1.0, 2.0, 3.0, 4.0]

    def fake_basal_melt_rate(icesee_kwargs, step, floating, Q, s, h, scenario, experiment=None):
        return melt_field, 0.0

    monkeypatch.setattr(native, "Icepack", fake_icepack)
    monkeypatch.setattr(native, "BasalMeltRate", fake_basal_melt_rate)

    icesee_kwargs = _base_icesee_kwargs(joint_estimation=True, wrong_basal_melt_field=2.0)
    topo = _topology()
    state = native.initialize_member(0, topology=topo, icesee_kwargs=icesee_kwargs)

    assert state.basal_melt is not None
    # Nudged: [1,2,3,4] + 2.0
    np.testing.assert_allclose(state.basal_melt.dat.data_ro, [3.0, 4.0, 5.0, 6.0])


def test_forecast_member_advances_persistent_state_without_ensemble_vector(monkeypatch):
    native._SHARED_CONTEXT_CACHE.clear()
    ctx_objs = _install_fake_mesh_pipeline(monkeypatch, build_calls=[])

    icepack_calls = []

    def fake_icepack(solver, h, u, smb, basal_melt_field, bed, dt, h0, run_kwargs):
        icepack_calls.append(dict(h=h, u=u, h0=h0, basal_melt_field=basal_melt_field))
        return ("h-next", "u-next", "s-next", "floating-next", "grounded-next")

    # 2026-09-28 reconciliation: forecast_member recomputes floating fresh
    # from the member's CURRENT surface each step (mf.flotationHeight/
    # mf.flotationMask), never the frozen ctx.floating0 an earlier version
    # of this adapter incorrectly kept reusing.
    persistent_surface = object()

    def fake_flotation_height(bed, Q):
        assert bed is ctx_objs["bed"]
        return "zF-sentinel"

    def fake_flotation_mask(surface, zF, Q):
        assert surface is persistent_surface
        assert zF == "zF-sentinel"
        return "floating-recomputed", "grounded-recomputed"

    def fake_basal_melt_rate(icesee_kwargs, step, floating, Q, s, h, scenario, experiment=None):
        assert floating == "floating-recomputed"
        assert experiment == native.EXPERIMENT_WRONG
        return "melt-field", 0.0

    monkeypatch.setattr(native, "Icepack", fake_icepack)
    monkeypatch.setattr(native, "BasalMeltRate", fake_basal_melt_rate)
    monkeypatch.setattr(native.mf, "flotationHeight", fake_flotation_height)
    monkeypatch.setattr(native.mf, "flotationMask", fake_flotation_mask)

    icesee_kwargs = _base_icesee_kwargs()
    topo = _topology()

    persistent_thickness = object()
    persistent_velocity = object()
    state = types.SimpleNamespace(
        thickness=persistent_thickness,
        velocity=persistent_velocity,
        surface=persistent_surface,
        basal_melt=None,
    )

    result = native.forecast_member(state, 5, topology=topo, icesee_kwargs=icesee_kwargs)

    assert isinstance(result, IcepackForecastFields)
    assert icepack_calls[0]["h"] is persistent_thickness
    assert icepack_calls[0]["u"] is persistent_velocity
    assert icepack_calls[0]["h0"] is ctx_objs["h0"]  # inflow BC is always the shared reference
    assert result.thickness == "h-next"
    assert result.velocity == "u-next"
    assert result.surface == "s-next"
    assert result.basal_melt is None  # state had no basal_melt block


def test_forecast_member_carries_basal_melt_only_when_member_has_it(monkeypatch):
    native._SHARED_CONTEXT_CACHE.clear()
    _install_fake_mesh_pipeline(monkeypatch, build_calls=[])

    def fake_icepack(solver, h, u, smb, basal_melt_field, bed, dt, h0, run_kwargs):
        return ("h-next", "u-next", "s-next", "floating-next", "grounded-next")

    def fake_basal_melt_rate(icesee_kwargs, step, floating, Q, s, h, scenario, experiment=None):
        return "melt-field", 0.0

    monkeypatch.setattr(native, "Icepack", fake_icepack)
    monkeypatch.setattr(native, "BasalMeltRate", fake_basal_melt_rate)
    monkeypatch.setattr(native.mf, "flotationHeight", lambda bed, Q: "zF-sentinel")
    monkeypatch.setattr(native.mf, "flotationMask", lambda surface, zF, Q: ("floating-recomputed", "grounded-recomputed"))

    icesee_kwargs = _base_icesee_kwargs()
    topo = _topology()
    state = types.SimpleNamespace(
        thickness=object(), velocity=object(), surface=object(), basal_melt=object()
    )
    result = native.forecast_member(state, 1, topology=topo, icesee_kwargs=icesee_kwargs)
    assert result.basal_melt == "melt-field"


def test_native_adapter_satisfies_required_distributed_protocol():
    validate_native_distributed_adapter(native.IDEALIZED_PIG_NATIVE_ADAPTER)
