import importlib.util
from pathlib import Path

import numpy as np


SCRIPT = (
    Path(__file__).parents[2]
    / "applications"
    / "issm_model"
    / "examples"
    / "ISMIP_Choi"
    / "generate_bed_kringing.py"
)


def _load_script():
    spec = importlib.util.spec_from_file_location("generate_bed_kringing", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_background_is_invariant_to_mesh_vertex_order():
    module = _load_script()
    rng = np.random.default_rng(17)
    x, y = np.meshgrid(np.linspace(0, 40_000, 9), np.linspace(0, 20_000, 5))
    x = x.ravel()
    y = y.ravel()
    bed = 200.0 * np.sin(x / 12_000.0) - 80.0 * np.cos(y / 7_000.0)

    expected = module.smooth_background_2d(
        bed, x, y, corr_len_m=8_000.0, truncation_sigma=3.0
    )
    permutation = rng.permutation(bed.size)
    shuffled = module.smooth_background_2d(
        bed[permutation],
        x[permutation],
        y[permutation],
        corr_len_m=8_000.0,
        truncation_sigma=3.0,
    )
    restored = np.empty_like(shuffled)
    restored[permutation] = shuffled

    np.testing.assert_allclose(restored, expected, rtol=0.0, atol=1.0e-12)


def test_background_preserves_a_constant_bed():
    module = _load_script()
    x = np.array([0.0, 1_000.0, 0.0, 1_000.0])
    y = np.array([0.0, 0.0, 1_000.0, 1_000.0])
    bed = np.full(4, -350.0)

    result = module.smooth_background_2d(
        bed, x, y, corr_len_m=2_000.0, truncation_sigma=3.0
    )
    np.testing.assert_allclose(result, bed)
