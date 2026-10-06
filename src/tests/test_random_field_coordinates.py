import numpy as np
import pytest

from ICESEE.src.run_model_da._error_generation import (
    generate_enkf_field,
    generate_pseudo_random_field_1d,
    generate_pseudo_random_field_2D,
)
from ICESEE.src.utils.localization import (
    prepare_random_field_coordinates,
    register_coord_provider,
)


def test_fft_does_not_fetch_registered_coordinates():
    model_name = "coordinate_test_fft"

    def provider(_):
        raise AssertionError("FFT must not request mesh coordinates")

    register_coord_provider(model_name, provider)
    kwargs = {"model_name": model_name, "random_field_method": "fft"}

    assert prepare_random_field_coordinates(kwargs, expected_nodes=4) is None
    assert "mesh_coords" not in kwargs


def test_graph_packs_registered_physical_coordinates():
    model_name = "coordinate_test_graph"
    expected = np.array(
        [[0.0, 0.0], [2.0, 0.0], [0.5, 3.0], [5.0, 4.0]]
    )
    register_coord_provider(model_name, lambda _: expected.copy())
    kwargs = {"model_name": model_name, "random_field_method": "graph"}

    result = prepare_random_field_coordinates(kwargs, expected_nodes=4)

    np.testing.assert_allclose(result, expected)
    np.testing.assert_allclose(kwargs["mesh_coords"], expected)


def test_graph_rejects_coordinate_count_mismatch():
    model_name = "coordinate_test_bad_size"
    register_coord_provider(model_name, lambda _: np.zeros((3, 2)))

    with pytest.raises(ValueError, match="coordinate count"):
        prepare_random_field_coordinates(
            {"model_name": model_name, "random_field_method": "graph"},
            expected_nodes=4,
        )


def test_graph_field_uses_one_physical_mesh_for_each_variable():
    coords = np.array(
        [[0.0, 0.0], [1.0, 0.0], [0.0, 2.0], [3.0, 2.0], [5.0, 4.0]]
    )

    field = generate_enkf_field(
        random_field_method="graph",
        mesh_coords=coords,
        noise_dim=coords.shape[0],
        num_vars=3,
        rh=2.0,
        num_passes=2,
    )

    assert field.shape == (15,)
    for block in field.reshape(3, 5):
        assert np.isfinite(block).all()
        np.testing.assert_allclose(np.mean(block), 0.0, atol=1.0e-12)
        np.testing.assert_allclose(np.std(block), 1.0, atol=1.0e-12)


def test_fft_length_scale_yaml_alias_reaches_rh():
    common = {
        "random_field_method": "fft",
        "noise_dim": 128,
        "num_vars": 1,
        "Lx_dim": 1000.0,
    }

    np.random.seed(23)
    from_yaml = generate_enkf_field(length_scale=[75.0], **common)
    np.random.seed(23)
    from_api = generate_enkf_field(rh=[75.0], **common)

    np.testing.assert_allclose(from_yaml, from_api)


def test_fft_explicit_scalar_reproduces_legacy_default_field():
    """Lock the pre-alias Icepack fallback without reintroducing the bug."""
    domain_scale = np.sqrt(5000.0 * 1200.0)
    common = {
        "random_field_method": "fft",
        "noise_dim": 128,
        "num_vars": 4,
        "Lx_dim": domain_scale,
    }

    np.random.seed(31)
    legacy_fallback = generate_enkf_field(**common)
    np.random.seed(31)
    explicit_effective_scale = generate_enkf_field(
        length_scale=domain_scale / 10.0,
        **common,
    )

    np.testing.assert_allclose(explicit_effective_scale, legacy_fallback)


def test_fft_scalar_and_per_variable_scales_have_distinct_semantics():
    """A repeated list must not be mistaken for the legacy scalar path."""
    common = {
        "random_field_method": "fft",
        "noise_dim": 128,
        "num_vars": 4,
        "Lx_dim": 1000.0,
    }

    np.random.seed(37)
    combined = generate_enkf_field(length_scale=100.0, **common)
    np.random.seed(37)
    independent_blocks = generate_enkf_field(
        length_scale=[100.0] * 4,
        **common,
    )

    assert not np.allclose(combined, independent_blocks)


def test_small_fft_field_honors_member_seed():
    """The low-dimensional fast path must not create an unseeded RNG."""
    common = {
        "random_field_method": "fft",
        "noise_dim": 3,
        "num_vars": 1,
        "Lx_dim": 10.0,
        "length_scale": 1.0,
    }

    first = generate_enkf_field(seed=8123, **common)
    repeated = generate_enkf_field(seed=8123, **common)
    other_member = generate_enkf_field(seed=8124, **common)

    np.testing.assert_array_equal(first, repeated)
    assert not np.array_equal(first, other_member)


def test_1d_generator_accepts_multidimensional_mesh_coordinates():
    coords = np.array(
        [[0.0, 0.0], [1.0, 0.0], [0.0, 2.0], [4.0, 2.0]]
    )
    field = generate_pseudo_random_field_1d(
        N=4, coords=coords, method="auto", num_passes=2
    )

    assert field.shape == (4,)
    assert np.isfinite(field).all()


def test_graph_field_changes_when_physical_mesh_geometry_changes():
    line_coords = np.column_stack((np.arange(8, dtype=float), np.zeros(8)))
    interleaved_coords = np.array(
        [
            [0.0, 0.0],
            [100.0, 0.0],
            [1.0, 0.0],
            [101.0, 0.0],
            [2.0, 0.0],
            [102.0, 0.0],
            [3.0, 0.0],
            [103.0, 0.0],
        ]
    )

    np.random.seed(17)
    line_field = generate_pseudo_random_field_1d(
        N=8,
        coords=line_coords,
        method="graph",
        k_neighbors=2,
        num_passes=3,
    )
    np.random.seed(17)
    interleaved_field = generate_pseudo_random_field_1d(
        N=8,
        coords=interleaved_coords,
        method="graph",
        k_neighbors=2,
        num_passes=3,
    )

    assert not np.allclose(line_field, interleaved_field)


def test_2d_generator_uses_explicit_mesh_coordinates():
    coords = np.array(
        [
            [0.0, 0.0],
            [2.0, 0.0],
            [5.0, 0.0],
            [0.5, 3.0],
            [2.5, 3.0],
            [6.0, 3.0],
        ]
    )
    field = generate_pseudo_random_field_2D(
        nx=3,
        ny=2,
        coords=coords,
        random_field_method="graph",
        num_passes=2,
    )

    assert field.shape == (2, 3)
    assert np.isfinite(field).all()


@pytest.mark.parametrize(
    "generator, kwargs",
    [
        (generate_pseudo_random_field_1d, {"N": 4}),
        (generate_pseudo_random_field_2D, {"nx": 2, "ny": 2}),
    ],
)
def test_low_level_generators_reject_unknown_method(generator, kwargs):
    with pytest.raises(ValueError, match="Unsupported random-field method"):
        generator(random_field_method="not-a-method", **kwargs)
